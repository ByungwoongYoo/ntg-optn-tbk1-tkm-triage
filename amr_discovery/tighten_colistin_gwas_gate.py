#!/usr/bin/env python3
"""Apply the final conservative gate to an evaluated colistin GWAS table.

This postprocessor is intentionally stricter than the exploratory evaluator. It requires
BH-corrected held-out replication, BH-corrected whole-cohort stability, random-effects
replication across source groups and countries, and exclusion of pre-screen-known
mechanism features. Passing remains an association, not causality or novelty.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--evidence", required=True)
    p.add_argument("--meta-details", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5)
    p.add_argument("--min-groups", type=int, default=3)
    return p.parse_args()


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if len(vals):
        order = np.argsort(vals)
        ranked = vals[order]
        q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out[np.flatnonzero(ok)[order]] = np.clip(q, 0, 1)
    return pd.Series(out, index=values.index)


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(a.evidence, dtype={"variant": str})
    details = json.loads(Path(a.meta_details).read_text())

    x["validation_q_bh"] = bh(x["validation_fisher_p"])
    for field in ["n_groups", "random_or", "random_ci_low", "random_ci_high", "random_p", "I2", "positive_groups"]:
        x[f"country_{field}"] = x.variant.map(lambda v: (details.get(str(v), {}).get("country") or {}).get(field))

    n_disc_models = int(pd.to_numeric(x.discovery_models, errors="coerce").max())
    n_whole_models = int(pd.to_numeric(x.whole_models, errors="coerce").max())
    known = x.get("known_mechanism_feature", False)
    if not isinstance(known, pd.Series):
        known = pd.Series(False, index=x.index)
    known = known.fillna(False).map(lambda v: str(v).strip().lower() in {"true", "1", "yes"})

    x["discovery_stable_final"] = (
        (x.discovery_models == n_disc_models)
        & (x.discovery_beta_min > 0)
        & (x.discovery_q_max <= a.alpha)
    )
    x["whole_stable_final"] = (
        (x.whole_models == n_whole_models)
        & (x.whole_beta_min > 0)
        & (x.whole_q_max <= a.alpha)
    )
    x["validation_replication_final"] = (
        (x.validation_R_present + x.validation_S_present >= a.min_validation_present)
        & (x.validation_or > 1)
        & (x.validation_ci_low > 1)
        & (x.validation_q_bh <= a.alpha)
    )
    x["source_replication_final"] = (
        (x.source_meta_n >= a.min_groups)
        & (pd.to_numeric(x.source_meta_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(x.source_meta_p, errors="coerce") <= a.alpha)
    )
    x["country_replication_final"] = (
        (pd.to_numeric(x.country_n_groups, errors="coerce") >= a.min_groups)
        & (pd.to_numeric(x.country_random_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(x.country_random_p, errors="coerce") <= a.alpha)
    )
    x["known_mechanism_feature_final"] = known
    x["strict_final_marker"] = (
        x.discovery_stable_final
        & x.whole_stable_final
        & x.validation_replication_final
        & x.source_replication_final
        & x.country_replication_final
        & ~known
    )
    x["final_evidence_score"] = (
        3 * x.discovery_stable_final.astype(int)
        + 2 * x.whole_stable_final.astype(int)
        + 3 * x.validation_replication_final.astype(int)
        + 3 * x.source_replication_final.astype(int)
        + 2 * x.country_replication_final.astype(int)
        - 2 * known.astype(int)
    )
    x = x.sort_values(
        ["strict_final_marker", "final_evidence_score", "discovery_q_max", "validation_q_bh"],
        ascending=[False, False, True, True],
    )
    x.to_csv(out / "ALL_MARKERS_FINAL_GATE.csv", index=False)
    strict = x[x.strict_final_marker].copy()
    strict.to_csv(out / "STRICT_FINAL_MARKERS.csv", index=False)

    summary = {
        "n_candidates_evaluated": int(len(x)),
        "n_discovery_stable": int(x.discovery_stable_final.sum()),
        "n_whole_stable": int(x.whole_stable_final.sum()),
        "n_validation_replicated_after_bh": int(x.validation_replication_final.sum()),
        "n_source_replicated": int(x.source_replication_final.sum()),
        "n_country_replicated": int(x.country_replication_final.sum()),
        "n_strict_final_markers": int(x.strict_final_marker.sum()),
        "strict_markers": strict.variant.astype(str).tolist(),
        "status": "CANDIDATES_REQUIRE_MECHANISM_NOVELTY_AND_BIOLOGICAL_VALIDATION" if len(strict) else "NO_STRICT_REPLICATED_MARKER",
        "boundary": "A surviving marker is a multiplicity-corrected, population-structure-aware, source-held-out association replicated across source groups and countries. It is not proof of causality, novelty, clinical validity, or treatment relevance.",
    }
    (out / "STRICT_GATE_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    cols = [c for c in [
        "variant", "gene", "type", "final_evidence_score", "discovery_q_max", "whole_q_max",
        "validation_or", "validation_ci_low", "validation_ci_high", "validation_q_bh",
        "source_meta_n", "source_meta_or", "source_meta_ci_low", "source_meta_ci_high", "source_meta_p",
        "country_n_groups", "country_random_or", "country_random_ci_low", "country_random_ci_high", "country_random_p",
        "known_mechanism_feature_final", "strict_final_marker"
    ] if c in x.columns]
    report = [
        "# Final conservative K. pneumoniae–colistin marker gate", "",
        f"- Candidate features evaluated: **{len(x):,}**",
        f"- Discovery-stable: **{int(x.discovery_stable_final.sum()):,}**",
        f"- Whole-cohort stable: **{int(x.whole_stable_final.sum()):,}**",
        f"- Held-out replication after BH correction: **{int(x.validation_replication_final.sum()):,}**",
        f"- Cross-source replication: **{int(x.source_replication_final.sum()):,}**",
        f"- Cross-country replication: **{int(x.country_replication_final.sum()):,}**",
        f"- Strict non-known survivors: **{int(x.strict_final_marker.sum()):,}", "",
        "## Claim boundary", "", summary["boundary"], "",
        "## Top evidence-ranked features", "", x.head(30)[cols].to_markdown(index=False), "",
        "## Final result", "",
        "At least one association advances to independent known-mechanism, sequence-context, database/literature, and biological-validation audits." if len(strict) else "No feature survives the complete final statistical gate; no new genomic-marker claim is supported by this analysis.", "",
    ]
    (out / "STRICT_GATE_REPORT.md").write_text("\n".join(report))
    hashes = [
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}"
        for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
