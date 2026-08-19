#!/usr/bin/env python3
"""Evaluate source-held-out replication and cross-stratum robustness of AMR markers.

The script combines pyseer discovery/whole-cohort results with a prespecified source-group
holdout. It reports strict replicated candidates, but never equates association with causality or
novelty. Literature/database novelty review and biological validation remain separate gates.
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True); p.add_argument("--rtab", required=True); p.add_argument("--feature-meta", required=True)
    p.add_argument("--discovery", nargs="+", required=True); p.add_argument("--whole", nargs="+", required=True)
    p.add_argument("--out", default="gwas_evaluation"); p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5); p.add_argument("--min-informative-groups", type=int, default=3)
    return p.parse_args()


def read_pyseer(path: str) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", dtype={"variant": str})
    for c in ["af", "filter-pvalue", "lrt-pvalue", "beta", "beta-std-err"]:
        if c in x.columns: x[c] = pd.to_numeric(x[c], errors="coerce")
    x["source_file"] = Path(path).name; return x


def bh(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(float); out = np.full(len(p), np.nan); ok = np.isfinite(p); vals = p[ok]
    if len(vals):
        order = np.argsort(vals); ranked = vals[order]; q = ranked * len(ranked) / np.arange(1, len(ranked) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; q = np.clip(q, 0, 1); out[np.flatnonzero(ok)[order]] = q
    return pd.Series(out, index=pvals.index)


def contingency(feature: pd.Series, meta: pd.DataFrame) -> tuple[int, int, int, int]:
    z = meta[["assembly_ID", "phenotype"]].copy(); z["x"] = z.assembly_ID.map(feature).fillna(0).astype(int); r = z.phenotype.astype(str) == "R"
    return int((r & (z.x == 1)).sum()), int((~r & (z.x == 1)).sum()), int((r & (z.x == 0)).sum()), int((~r & (z.x == 0)).sum())


def or_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float, float]:
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0: aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
    lor = math.log((aa * dd) / (bb * cc)); se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return math.exp(lor), math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se), se


def fixed_random_meta(rows: list[tuple[str, int, int, int, int]]) -> dict[str, Any]:
    effects = []
    for group, a, b, c, d in rows:
        if (a + c) == 0 or (b + d) == 0 or (a + b) == 0 or (c + d) == 0: continue
        aa, bb, cc, dd = map(float, (a, b, c, d))
        if min(aa, bb, cc, dd) == 0: aa += 0.5; bb += 0.5; cc += 0.5; dd += 0.5
        lor = math.log((aa * dd) / (bb * cc)); var = 1 / aa + 1 / bb + 1 / cc + 1 / dd; effects.append((group, lor, var, a, b, c, d))
    if not effects: return {"n_groups": 0}
    w = np.array([1 / e[2] for e in effects]); y = np.array([e[1] for e in effects]); mu = float(np.sum(w * y) / np.sum(w)); q = float(np.sum(w * (y - mu) ** 2)); df = len(effects) - 1; cval = float(np.sum(w) - np.sum(w ** 2) / np.sum(w)); tau2 = max(0.0, (q - df) / cval) if df > 0 and cval > 0 else 0.0
    wr = 1 / (np.array([e[2] for e in effects]) + tau2); mur = float(np.sum(wr * y) / np.sum(wr)); ser = math.sqrt(1 / float(np.sum(wr))); z = mur / ser if ser > 0 else 0.0; p = float(2 * norm.sf(abs(z))); i2 = max(0.0, (q - df) / q * 100) if q > 0 and df > 0 else 0.0
    return {"n_groups": len(effects), "random_or": math.exp(mur), "random_ci_low": math.exp(mur - 1.96 * ser), "random_ci_high": math.exp(mur + 1.96 * ser), "random_p": p, "tau2": tau2, "I2": i2, "positive_groups": int(sum(e[1] > 0 for e in effects)), "group_details": [{"group": e[0], "log_or": e[1], "variance": e[2], "a_R_present": e[3], "b_S_present": e[4], "c_R_absent": e[5], "d_S_absent": e[6]} for e in effects]}


def meta_by(feature: pd.Series, meta: pd.DataFrame, group_col: str) -> dict[str, Any]:
    return fixed_random_meta([(str(g), *contingency(feature, sub)) for g, sub in meta.groupby(group_col, dropna=False)])


def main() -> None:
    a = parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(a.manifest, dtype={"assembly_ID": str}); rtab = pd.read_csv(a.rtab, sep="\t", index_col=0); rtab.columns = rtab.columns.astype(str); rtab = rtab.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    fmap = pd.read_csv(a.feature_meta, dtype={"feature": str}).drop_duplicates("feature")
    disc_frames = []
    for i, pth in enumerate(a.discovery):
        x = read_pyseer(pth); x["model_index"] = i; x["q_bh"] = bh(x["lrt-pvalue"]); disc_frames.append(x)
    whole_frames = []
    for i, pth in enumerate(a.whole):
        x = read_pyseer(pth); x["model_index"] = i; x["q_bh"] = bh(x["lrt-pvalue"]); whole_frames.append(x)
    disc = pd.concat(disc_frames, ignore_index=True); whole = pd.concat(whole_frames, ignore_index=True); disc.to_csv(out / "discovery_pyseer_combined.tsv", sep="\t", index=False); whole.to_csv(out / "whole_pyseer_combined.tsv", sep="\t", index=False)
    n_disc_models = disc.model_index.nunique(); n_whole_models = whole.model_index.nunique()
    dagg = disc.groupby("variant").agg(discovery_models=("model_index", "nunique"), discovery_beta_min=("beta", "min"), discovery_beta_max=("beta", "max"), discovery_p_max=("lrt-pvalue", "max"), discovery_p_min=("lrt-pvalue", "min"), discovery_q_max=("q_bh", "max"), discovery_q_min=("q_bh", "min")).reset_index()
    wagg = whole.groupby("variant").agg(whole_models=("model_index", "nunique"), whole_beta_min=("beta", "min"), whole_beta_max=("beta", "max"), whole_p_max=("lrt-pvalue", "max"), whole_p_min=("lrt-pvalue", "min"), whole_q_max=("q_bh", "max"), whole_q_min=("q_bh", "min")).reset_index()
    candidates = dagg.merge(wagg, on="variant", how="outer").merge(fmap, left_on="variant", right_on="feature", how="left")
    val = manifest[manifest.split == "validation"].copy(); allm = manifest.copy(); results = []; meta_details = {}
    for _, row in candidates.iterrows():
        feat = str(row.variant)
        if feat not in rtab.index: continue
        series = rtab.loc[feat]; va_a, va_b, va_c, va_d = contingency(series, val); va_or, va_lo, va_hi, _ = or_ci(va_a, va_b, va_c, va_d)
        try: va_p = float(fisher_exact([[va_a, va_b], [va_c, va_d]], alternative="two-sided").pvalue)
        except Exception: va_p = float("nan")
        al_a, al_b, al_c, al_d = contingency(series, allm); al_or, al_lo, al_hi, _ = or_ci(al_a, al_b, al_c, al_d)
        try: al_p = float(fisher_exact([[al_a, al_b], [al_c, al_d]], alternative="two-sided").pvalue)
        except Exception: al_p = float("nan")
        source_meta = meta_by(series, allm, "source_group") if "source_group" in allm.columns else {"n_groups": 0}; country_meta = meta_by(series, allm, "ISO_country_code") if "ISO_country_code" in allm.columns else {"n_groups": 0}; st_col = "Kleborate_ST" if "Kleborate_ST" in allm.columns else None; st_meta = meta_by(series, allm, st_col) if st_col else {"n_groups": 0}
        meta_details[feat] = {"source_group": source_meta, "country": country_meta, "ST": st_meta}; d = row.to_dict(); d.update({"validation_R_present": va_a, "validation_S_present": va_b, "validation_R_absent": va_c, "validation_S_absent": va_d, "validation_or": va_or, "validation_ci_low": va_lo, "validation_ci_high": va_hi, "validation_fisher_p": va_p, "all_R_present": al_a, "all_S_present": al_b, "all_R_absent": al_c, "all_S_absent": al_d, "all_unadjusted_or": al_or, "all_unadjusted_ci_low": al_lo, "all_unadjusted_ci_high": al_hi, "all_fisher_p": al_p, "source_meta_n": source_meta.get("n_groups", 0), "source_meta_or": source_meta.get("random_or"), "source_meta_ci_low": source_meta.get("random_ci_low"), "source_meta_ci_high": source_meta.get("random_ci_high"), "source_meta_p": source_meta.get("random_p"), "source_meta_I2": source_meta.get("I2"), "source_meta_positive_groups": source_meta.get("positive_groups", 0), "country_meta_n": country_meta.get("n_groups", 0), "country_meta_or": country_meta.get("random_or"), "country_meta_ci_low": country_meta.get("random_ci_low"), "country_meta_ci_high": country_meta.get("random_ci_high"), "ST_meta_n": st_meta.get("n_groups", 0), "ST_meta_or": st_meta.get("random_or"), "ST_meta_ci_low": st_meta.get("random_ci_low"), "ST_meta_ci_high": st_meta.get("random_ci_high")}); results.append(d)
    res = pd.DataFrame(results)
    if res.empty: raise RuntimeError("No pyseer variants overlapped the supplied Rtab matrix")
    raw_known = res.get("known_mechanism_screen", pd.Series(False, index=res.index)).fillna(False); known = raw_known.map(lambda x: str(x).strip().lower() in {"true", "1", "yes"})
    res["discovery_stable"] = (res.discovery_models == n_disc_models) & (res.discovery_beta_min > 0) & (res.discovery_q_max <= a.alpha)
    res["whole_stable"] = (res.whole_models == n_whole_models) & (res.whole_beta_min > 0) & (res.whole_p_max <= a.alpha)
    res["validation_replication"] = (res.validation_R_present + res.validation_S_present >= a.min_validation_present) & (res.validation_or > 1) & (res.validation_ci_low > 1) & (res.validation_fisher_p <= a.alpha)
    res["source_replication"] = (res.source_meta_n >= a.min_informative_groups) & (pd.to_numeric(res.source_meta_ci_low, errors="coerce") > 1) & (pd.to_numeric(res.source_meta_p, errors="coerce") <= a.alpha)
    res["strict_replicated_marker"] = res.discovery_stable & res.whole_stable & res.validation_replication & res.source_replication & ~known; res["known_mechanism_feature"] = known
    res["evidence_score"] = res.discovery_stable.astype(int) * 3 + res.whole_stable.astype(int) * 2 + res.validation_replication.astype(int) * 3 + res.source_replication.astype(int) * 3 + (pd.to_numeric(res.source_meta_ci_low, errors="coerce") > 1).fillna(False).astype(int) - known.astype(int) * 2
    res = res.sort_values(["strict_replicated_marker", "evidence_score", "discovery_p_max", "validation_fisher_p"], ascending=[False, False, True, True]); res.to_csv(out / "ALL_MARKER_EVIDENCE.csv", index=False); strict = res[res.strict_replicated_marker].copy(); strict.to_csv(out / "STRICT_REPLICATED_MARKERS.csv", index=False); (out / "MARKER_META_ANALYSIS_DETAILS.json").write_text(json.dumps(meta_details, indent=2, ensure_ascii=False, default=str) + "\n")
    summary = {"n_features_in_matrix": int(len(rtab)), "n_discovery_models": int(n_disc_models), "n_whole_models": int(n_whole_models), "n_discovery_stable": int(res.discovery_stable.sum()), "n_validation_replicated": int(res.validation_replication.sum()), "n_source_replicated": int(res.source_replication.sum()), "n_strict_replicated_markers": int(res.strict_replicated_marker.sum()), "strict_markers": strict.variant.astype(str).tolist(), "status": "CANDIDATES_REQUIRE_NOVELTY_AND_BIOLOGICAL_VALIDATION" if len(strict) else "NO_STRICT_REPLICATED_MARKER", "boundary": "A strict marker is an association that survived source-held-out replication and population-structure sensitivity. It is not proof of causality, a new mechanism, clinical validity, or treatment relevance."}
    (out / "GWAS_FINAL_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    top = res.head(30); cols = [c for c in ["variant", "gene", "type", "evidence_score", "discovery_p_max", "validation_or", "validation_ci_low", "validation_ci_high", "validation_fisher_p", "source_meta_n", "source_meta_or", "source_meta_ci_low", "source_meta_ci_high", "source_meta_I2", "known_mechanism_feature", "strict_replicated_marker"] if c in top.columns]
    report = ["# K. pneumoniae–colistin targeted marker evaluation", "", f"- Features evaluated: **{len(res):,}**", f"- Stable in all discovery models: **{int(res.discovery_stable.sum()):,}**", f"- Replicated in the source-held-out validation set: **{int(res.validation_replication.sum()):,}**", f"- Replicated across independent source groups: **{int(res.source_replication.sum()):,}**", f"- Strict replicated non-screen-known markers: **{int(res.strict_replicated_marker.sum()):,}**", "", "## Claim boundary", "", summary["boundary"], "", "## Top evidence-ranked features", "", top[cols].to_markdown(index=False), "", "## Final result", "", "One or more markers met the statistical replication gate. They advance only to database/literature novelty audit and biological interpretation; they are not yet novel resistance determinants." if len(strict) else "No feature met the complete prespecified replication gate. The public-data analysis therefore does not support a new genomic-marker claim at this stage.", ""]
    (out / "GWAS_FINAL_REPORT.md").write_text("\n".join(report)); hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]; (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
