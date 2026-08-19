#!/usr/bin/env python3
"""Evaluate discovery-selected unitigs in a source-held-out validation partition.

Final statistical promotion requires independent one-sided validation with multiplicity control,
positive population-structure-adjusted whole-cohort effects across sensitivity models, and
cross-source replication. It remains an association screen, not a causal or novelty claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--selection", required=True); p.add_argument("--manifest", required=True)
    p.add_argument("--validation-rtab", required=True); p.add_argument("--all-rtab", required=True)
    p.add_argument("--whole-pyseer", nargs="+", required=True); p.add_argument("--out", default="unitig_replication")
    p.add_argument("--alpha", type=float, default=0.05); p.add_argument("--min-validation-present", type=int, default=5)
    p.add_argument("--min-source-groups", type=int, default=3); p.add_argument("--min-st-groups", type=int, default=2)
    return p.parse_args()


def bh(s: pd.Series) -> pd.Series:
    p = pd.to_numeric(s, errors="coerce").to_numpy(float); out = np.full(len(p), np.nan); ok = np.isfinite(p); v = p[ok]
    if len(v):
        o = np.argsort(v); r = v[o]; q = np.minimum.accumulate((r * len(r) / np.arange(1, len(r) + 1))[::-1])[::-1]; q = np.clip(q, 0, 1); out[np.flatnonzero(ok)[o]] = q
    return pd.Series(out, index=s.index)


def read_rtab(path: str) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0); x.index = x.index.astype(str).str.upper(); x.columns = x.columns.astype(str); return x.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def contingency(row: pd.Series, meta: pd.DataFrame) -> tuple[int, int, int, int]:
    z = meta[["assembly_ID", "phenotype"]].copy(); z["x"] = z.assembly_ID.map(row).fillna(0).astype(int); r = z.phenotype.astype(str) == "R"
    return int((r & (z.x == 1)).sum()), int((~r & (z.x == 1)).sum()), int((r & (z.x == 0)).sum()), int((~r & (z.x == 0)).sum())


def or_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0: aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
    l = math.log((aa * dd) / (bb * cc)); se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd); return math.exp(l), math.exp(l - 1.96 * se), math.exp(l + 1.96 * se)


def meta_random(row: pd.Series, meta: pd.DataFrame, group_col: str) -> dict[str, Any]:
    effects = []
    for g, sub in meta.groupby(group_col, dropna=False):
        a, b, c, d = contingency(row, sub)
        if (a + c) == 0 or (b + d) == 0 or (a + b) == 0 or (c + d) == 0: continue
        aa, bb, cc, dd = map(float, (a, b, c, d))
        if min(aa, bb, cc, dd) == 0: aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
        lor = math.log((aa * dd) / (bb * cc)); var = 1 / aa + 1 / bb + 1 / cc + 1 / dd; effects.append((str(g), lor, var, a, b, c, d))
    if not effects: return {"n_groups": 0}
    y = np.array([e[1] for e in effects]); w = 1 / np.array([e[2] for e in effects]); mu = float(np.sum(w * y) / np.sum(w)); q = float(np.sum(w * (y - mu) ** 2)); df = len(y) - 1; cval = float(np.sum(w) - np.sum(w ** 2) / np.sum(w)); tau2 = max(0.0, (q - df) / cval) if df > 0 and cval > 0 else 0.0
    wr = 1 / (np.array([e[2] for e in effects]) + tau2); mur = float(np.sum(wr * y) / np.sum(wr)); se = math.sqrt(1 / float(np.sum(wr))); p = float(norm.sf(mur / se)) if se > 0 else 1.0; i2 = max(0.0, (q - df) / q * 100) if q > 0 and df > 0 else 0.0
    return {"n_groups": len(effects), "random_or": math.exp(mur), "ci_low": math.exp(mur - 1.96 * se), "ci_high": math.exp(mur + 1.96 * se), "one_sided_p": p, "I2": i2, "positive_groups": int(sum(e[1] > 0 for e in effects)), "details": [{"group": e[0], "log_or": e[1], "variance": e[2], "a": e[3], "b": e[4], "c": e[5], "d": e[6]} for e in effects]}


def read_pyseer(path: str, idx: int) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", dtype={"variant": str}); x.variant = x.variant.astype(str).str.upper(); x.beta = pd.to_numeric(x.beta, errors="coerce"); x["lrt-pvalue"] = pd.to_numeric(x["lrt-pvalue"], errors="coerce"); x["q_bh"] = bh(x["lrt-pvalue"])
    return x[["variant", "beta", "lrt-pvalue", "q_bh"]].rename(columns={"beta": f"whole_beta_{idx}", "lrt-pvalue": f"whole_p_{idx}", "q_bh": f"whole_q_{idx}"}).drop_duplicates("variant")


def main() -> None:
    a = args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sel = pd.read_csv(a.selection, dtype={"canonical_sequence": str, "candidate_id": str}); sel.canonical_sequence = sel.canonical_sequence.str.upper(); meta = pd.read_csv(a.manifest, dtype={"assembly_ID": str}); val_meta = meta[meta.split == "validation"].copy(); vr = read_rtab(a.validation_rtab); ar = read_rtab(a.all_rtab)
    whole = None
    for i, p in enumerate(a.whole_pyseer):
        x = read_pyseer(p, i); whole = x if whole is None else whole.merge(x, on="variant", how="outer")
    assert whole is not None
    bcols = [c for c in whole if c.startswith("whole_beta_")]; pcols = [c for c in whole if c.startswith("whole_p_")]; qcols = [c for c in whole if c.startswith("whole_q_")]
    whole["whole_models_present"] = whole[bcols].notna().sum(axis=1); whole["whole_beta_min"] = whole[bcols].min(axis=1); whole["whole_p_max"] = whole[pcols].max(axis=1); whole["whole_q_max"] = whole[qcols].max(axis=1)
    rows = []; details = {}; nmodels = len(a.whole_pyseer)
    for _, s in sel.iterrows():
        seq = s.canonical_sequence
        if seq not in vr.index or seq not in ar.index: continue
        vrow = vr.loc[seq]; arow = ar.loc[seq]; va = contingency(vrow, val_meta); ao = contingency(arow, meta); vor, vlo, vhi = or_ci(*va); aor, alo, ahi = or_ci(*ao)
        vp = float(fisher_exact([[va[0], va[1]], [va[2], va[3]]], alternative="greater").pvalue); ap = float(fisher_exact([[ao[0], ao[1]], [ao[2], ao[3]]], alternative="greater").pvalue)
        src = meta_random(arow, meta, "source_group"); st_col = "Kleborate_ST" if "Kleborate_ST" in meta.columns else None; stm = meta_random(arow, meta, st_col) if st_col else {"n_groups": 0}; country = meta_random(arow, meta, "ISO_country_code") if "ISO_country_code" in meta.columns else {"n_groups": 0}; details[str(s.candidate_id)] = {"source": src, "ST": stm, "country": country}
        d = s.to_dict(); d.update({"variant": seq, "validation_R_present": va[0], "validation_S_present": va[1], "validation_R_absent": va[2], "validation_S_absent": va[3], "validation_or": vor, "validation_ci_low": vlo, "validation_ci_high": vhi, "validation_p_one_sided": vp, "all_R_present": ao[0], "all_S_present": ao[1], "all_R_absent": ao[2], "all_S_absent": ao[3], "all_unadjusted_or": aor, "all_unadjusted_ci_low": alo, "all_unadjusted_ci_high": ahi, "all_p_one_sided": ap, "source_n": src.get("n_groups", 0), "source_or": src.get("random_or"), "source_ci_low": src.get("ci_low"), "source_ci_high": src.get("ci_high"), "source_p": src.get("one_sided_p"), "source_I2": src.get("I2"), "source_positive_groups": src.get("positive_groups", 0), "st_n": stm.get("n_groups", 0), "st_or": stm.get("random_or"), "st_ci_low": stm.get("ci_low"), "st_ci_high": stm.get("ci_high"), "country_n": country.get("n_groups", 0), "country_or": country.get("random_or"), "country_ci_low": country.get("ci_low"), "country_ci_high": country.get("ci_high")}); rows.append(d)
    res = pd.DataFrame(rows)
    if res.empty: raise RuntimeError("No selected unitig was present in both validation and all-cohort matrices")
    res["validation_q_bh"] = bh(res.validation_p_one_sided); res = res.merge(whole, on="variant", how="left")
    res["validation_replication"] = (res.validation_R_present + res.validation_S_present >= a.min_validation_present) & (res.validation_or > 1) & (res.validation_ci_low > 1) & (res.validation_q_bh <= a.alpha)
    res["whole_adjusted_stable"] = (res.whole_models_present == nmodels) & (res.whole_beta_min > 0) & (res.whole_q_max <= a.alpha)
    res["source_replication"] = (res.source_n >= a.min_source_groups) & (pd.to_numeric(res.source_ci_low, errors="coerce") > 1) & (pd.to_numeric(res.source_p, errors="coerce") <= a.alpha)
    res["st_replication"] = (res.st_n >= a.min_st_groups) & (pd.to_numeric(res.st_ci_low, errors="coerce") > 1)
    res["strict_statistical_replication"] = res.validation_replication & res.whole_adjusted_stable & res.source_replication & res.st_replication
    res["evidence_score"] = res.validation_replication.astype(int) * 4 + res.whole_adjusted_stable.astype(int) * 3 + res.source_replication.astype(int) * 4 + res.st_replication.astype(int) * 2 + (pd.to_numeric(res.country_ci_low, errors="coerce") > 1).fillna(False).astype(int)
    res = res.sort_values(["strict_statistical_replication", "evidence_score", "validation_q_bh", "whole_q_max"], ascending=[False, False, True, True]); res.to_csv(out / "ALL_UNITIG_REPLICATION_EVIDENCE.csv", index=False); strict = res[res.strict_statistical_replication].copy(); strict.to_csv(out / "STRICT_STATISTICALLY_REPLICATED_UNITIGS.csv", index=False); (out / "UNITIG_META_ANALYSIS_DETAILS.json").write_text(json.dumps(details, indent=2, ensure_ascii=False, default=str) + "\n")
    summary = {"n_selected": int(len(sel)), "n_evaluated": int(len(res)), "n_validation_replicated": int(res.validation_replication.sum()), "n_whole_adjusted_stable": int(res.whole_adjusted_stable.sum()), "n_source_replicated": int(res.source_replication.sum()), "n_strict_statistically_replicated": int(res.strict_statistical_replication.sum()), "strict_candidate_ids": strict.candidate_id.astype(str).tolist(), "status": "STATISTICALLY_REPLICATED_UNITIGS_REQUIRE_CONTEXT_AND_NOVELTY_AUDIT" if len(strict) else "NO_UNITIG_SURVIVED_COMPLETE_REPLICATION_GATE", "boundary": "Statistical replication does not establish a new resistance determinant. Strict unitigs must still be mapped to genomic context, compared with known colistin mechanisms and databases, checked for linkage/lineage confounding, and validated biologically."}
    (out / "UNITIG_REPLICATION_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    cols = [c for c in ["candidate_id", "sequence_length", "p_max", "q_max", "validation_or", "validation_ci_low", "validation_ci_high", "validation_q_bh", "whole_beta_min", "whole_q_max", "source_n", "source_or", "source_ci_low", "source_ci_high", "source_I2", "st_n", "st_or", "st_ci_low", "st_ci_high", "strict_statistical_replication"] if c in res.columns]
    report = ["# Whole-genome unitig replication audit", "", f"- Discovery-selected canonical unitigs: **{len(sel):,}**", f"- Evaluated in held-out data: **{len(res):,}**", f"- Held-out validation replicates after BH correction: **{int(res.validation_replication.sum()):,}**", f"- Positive in all population-structure sensitivity models: **{int(res.whole_adjusted_stable.sum()):,}**", f"- Cross-source random-effects replicates: **{int(res.source_replication.sum()):,}**", f"- Complete statistical gate: **{int(res.strict_statistical_replication.sum()):,}**", "", "## Claim boundary", "", summary["boundary"], "", "## Highest-ranked unitigs", "", res.head(30)[cols].to_markdown(index=False), ""]
    (out / "UNITIG_REPLICATION_REPORT.md").write_text("\n".join(report)); hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]; (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__": main()
