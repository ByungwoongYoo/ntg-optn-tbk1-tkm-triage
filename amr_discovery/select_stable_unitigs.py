#!/usr/bin/env python3
"""Select a prespecified, direction-stable unitig screen for independent validation.

Unitigs are selected only from the discovery partition. Reverse-complement duplicates are
collapsed. Selection is a screening step, not a resistance or novelty claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pyseer", nargs="+", required=True)
    p.add_argument("--out", default="unitig_candidates")
    p.add_argument("--alpha-screen", type=float, default=0.10)
    p.add_argument("--p-screen", type=float, default=1e-5)
    p.add_argument("--max-candidates", type=int, default=5000)
    p.add_argument("--min-length", type=int, default=31)
    p.add_argument("--max-length", type=int, default=5000)
    return p.parse_args()


def bh(s: pd.Series) -> pd.Series:
    p = pd.to_numeric(s, errors="coerce").to_numpy(float); out = np.full(len(p), np.nan); ok = np.isfinite(p); v = p[ok]
    if len(v):
        order = np.argsort(v); r = v[order]; q = np.minimum.accumulate((r * len(r) / np.arange(1, len(r) + 1))[::-1])[::-1]; q = np.clip(q, 0, 1); out[np.flatnonzero(ok)[order]] = q
    return pd.Series(out, index=s.index)


def canonical(seq: str) -> tuple[str, str]:
    s = seq.upper().replace("U", "T"); rc = str(Seq(s).reverse_complement())
    return (s, "+") if s <= rc else (rc, "-")


def read_one(path: str, idx: int) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", dtype={"variant": str})
    need = {"variant", "beta", "lrt-pvalue"}
    if not need.issubset(x.columns): raise RuntimeError(f"{path} lacks {sorted(need - set(x.columns))}")
    x["beta"] = pd.to_numeric(x.beta, errors="coerce"); x["lrt-pvalue"] = pd.to_numeric(x["lrt-pvalue"], errors="coerce"); x["q_bh"] = bh(x["lrt-pvalue"])
    x = x[["variant", "beta", "lrt-pvalue", "q_bh"]].rename(columns={"beta": f"beta_{idx}", "lrt-pvalue": f"p_{idx}", "q_bh": f"q_{idx}"})
    return x.drop_duplicates("variant", keep="first")


def main() -> None:
    a = args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    merged = None
    for i, p in enumerate(a.pyseer):
        x = read_one(p, i); merged = x if merged is None else merged.merge(x, on="variant", how="inner", validate="one_to_one")
    assert merged is not None
    bcols = [c for c in merged if c.startswith("beta_")]; pcols = [c for c in merged if c.startswith("p_")]; qcols = [c for c in merged if c.startswith("q_")]
    merged["beta_min"] = merged[bcols].min(axis=1); merged["p_max"] = merged[pcols].max(axis=1); merged["q_max"] = merged[qcols].max(axis=1)
    merged["direction_stable_positive"] = merged.beta_min > 0
    merged["discovery_screen"] = merged.direction_stable_positive & ((merged.q_max <= a.alpha_screen) | (merged.p_max <= a.p_screen))
    merged["sequence_length"] = merged.variant.astype(str).str.len(); merged = merged[(merged.sequence_length >= a.min_length) & (merged.sequence_length <= a.max_length)].copy()
    can = merged.variant.astype(str).map(canonical); merged["canonical_sequence"] = [x[0] for x in can]; merged["original_orientation"] = [x[1] for x in can]
    merged = merged.sort_values(["discovery_screen", "p_max", "q_max", "sequence_length"], ascending=[False, True, True, True])
    best = merged.drop_duplicates("canonical_sequence", keep="first").copy(); screened = best[best.discovery_screen].head(a.max_candidates).copy(); screened["candidate_id"] = [f"UG{n:06d}" for n in range(1, len(screened) + 1)]
    merged.to_csv(out / "ALL_DISCOVERY_UNITIG_RESULTS.csv", index=False); best.to_csv(out / "CANONICAL_DISCOVERY_UNITIGS.csv", index=False); screened.to_csv(out / "SELECTED_DISCOVERY_UNITIGS.csv", index=False)
    with open(out / "selected_unitigs.fasta", "w") as fh:
        for _, r in screened.iterrows(): fh.write(f">{r.candidate_id}\n{r.canonical_sequence}\n")
    with open(out / "selected_unitigs.tsv", "w") as fh:
        fh.write("candidate_id\tsequence\n");
        for _, r in screened.iterrows(): fh.write(f"{r.candidate_id}\t{r.canonical_sequence}\n")
    summary = {"n_models": len(a.pyseer), "n_intersection": int(len(merged)), "n_direction_stable_positive": int(merged.direction_stable_positive.sum()), "n_discovery_screen_before_rc_collapse": int(merged.discovery_screen.sum()), "n_selected_after_rc_collapse": int(len(screened)), "alpha_screen": a.alpha_screen, "p_screen": a.p_screen, "max_candidates": a.max_candidates, "status": "CANDIDATES_SELECTED_FOR_INDEPENDENT_VALIDATION" if len(screened) else "NO_DISCOVERY_UNITIG_PASSED_SCREEN", "boundary": "Selection used the discovery partition only. Candidates require independent held-out replication, population-structure sensitivity, genomic-context annotation, known-mechanism exclusion, and external review."}
    (out / "UNITIG_SELECTION_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]; (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__": main()
