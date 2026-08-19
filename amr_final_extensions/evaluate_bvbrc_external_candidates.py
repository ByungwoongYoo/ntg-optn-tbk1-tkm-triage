#!/usr/bin/env python3
"""Evaluate frozen AMR Portal candidates in a laboratory-only BV-BRC cohort.

The candidate list originates entirely from the earlier AMR Portal discovery workflow.
The BV-BRC cohort is BioProject-disjoint and was frozen before candidate sequence
queries. Records labelled `Computational Method` are prohibited. This script performs
no candidate re-selection: it tests the 22 targeted variants and 61 unitigs exactly as
frozen, with Mash-PC score tests, one-sided Fisher evidence, and source/country
robustness summaries.

A statistical survivor is an independently replicated association, not a causal,
clinically validated, or novel resistance mechanism. Known mgrB disruption is retained
only as a positive control and is never counted as a novel discovery.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amr_discovery.run_structure_adjusted_score_gwas import (
    bh,
    contingency,
    odds_ratio_ci,
    pcoa,
    score_test_matrix,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True)
    p.add_argument("--candidate-panel", required=True)
    p.add_argument("--targeted-rtab", required=True)
    p.add_argument("--unitig-rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pc-dims", default="10,20,30")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-af", type=float, default=0.01)
    p.add_argument("--max-af", type=float, default=0.99)
    p.add_argument("--min-carriers", type=int, default=5)
    p.add_argument("--min-odds-ratio", type=float, default=2.0)
    p.add_argument("--min-groups", type=int, default=3)
    return p.parse_args()


def load_rtab(path: str | Path) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = x.index.astype(str)
    x.columns = x.columns.astype(str)
    return x.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)


def clean_group(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if text.lower() in {"", "nan", "none", "unknown", "unknown_source", "na", "n/a"}:
        return ""
    return text


def group_series(cohort: pd.DataFrame, candidates: list[str]) -> pd.Series:
    out = pd.Series([""] * len(cohort), index=cohort.index, dtype=str)
    for column in candidates:
        if column not in cohort:
            continue
        values = cohort[column].map(clean_group)
        out = out.where(out.str.len() > 0, values)
    return out


def random_effects_informative(
    x: np.ndarray,
    y: np.ndarray,
    groups: pd.Series,
) -> dict:
    effects: list[tuple[str, float, float, int, int, int, int]] = []
    gv = groups.map(clean_group).to_numpy(str)
    for group in sorted(set(gv) - {""}):
        idx = gv == group
        a, b, c, d = contingency(x[idx], y[idx])
        if min(a + c, b + d, a + b, c + d) == 0:
            continue
        aa, bb, cc, dd = map(float, (a, b, c, d))
        if min(aa, bb, cc, dd) == 0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5
        log_or = math.log((aa * dd) / (bb * cc))
        variance = 1 / aa + 1 / bb + 1 / cc + 1 / dd
        effects.append((group, log_or, variance, a, b, c, d))
    if not effects:
        return {"n_groups": 0, "details": []}
    yi = np.array([row[1] for row in effects], dtype=float)
    vi = np.array([row[2] for row in effects], dtype=float)
    wi = 1 / vi
    fixed = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * (yi - fixed) ** 2))
    df = len(effects) - 1
    cval = float(np.sum(wi) - np.sum(wi * wi) / np.sum(wi))
    tau2 = max(0.0, (q - df) / cval) if df > 0 and cval > 0 else 0.0
    wr = 1 / (vi + tau2)
    pooled = float(np.sum(wr * yi) / np.sum(wr))
    se = math.sqrt(1 / float(np.sum(wr)))
    z = pooled / se if se else 0.0
    return {
        "n_groups": len(effects),
        "odds_ratio": math.exp(pooled),
        "ci_low": math.exp(pooled - 1.96 * se),
        "ci_high": math.exp(pooled + 1.96 * se),
        "one_sided_p": float(norm.sf(z)),
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 and df > 0 else 0.0,
        "positive_groups": int(sum(row[1] > 0 for row in effects)),
        "details": [
            {
                "group": row[0], "log_or": row[1], "variance": row[2],
                "R_present": row[3], "S_present": row[4],
                "R_absent": row[5], "S_absent": row[6],
            }
            for row in effects
        ],
    }


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    dims = [int(value) for value in args.pc_dims.split(",") if value.strip()]

    cohort = pd.read_csv(args.cohort, dtype=str).fillna("")
    required = {"genome_id", "phenotype", "evidence"}
    if missing := required - set(cohort.columns):
        raise ValueError(f"Cohort missing {sorted(missing)}")
    cohort = cohort.drop_duplicates("genome_id").copy()
    if not cohort["evidence"].eq("Laboratory Method").all():
        bad = cohort.loc[~cohort["evidence"].eq("Laboratory Method"), ["genome_id", "evidence"]]
        raise ValueError(f"Non-laboratory phenotype records are prohibited: {bad.head().to_dict('records')}")
    if not cohort["phenotype"].isin(["R", "S"]).all():
        raise ValueError("Non-binary phenotype in external cohort")
    if "strict_bioproject_disjoint" in cohort and not cohort["strict_bioproject_disjoint"].astype(str).str.lower().eq("true").all():
        raise ValueError("External cohort is not uniformly BioProject-disjoint")

    panel = pd.read_csv(args.candidate_panel, dtype=str).fillna("")
    if panel["candidate_id"].duplicated().any() or len(panel) != 83:
        raise ValueError(f"Frozen candidate panel is not the expected unique 83-candidate set: n={len(panel)}")
    targeted = load_rtab(args.targeted_rtab)
    unitigs = load_rtab(args.unitig_rtab)
    distance = pd.read_csv(args.distance, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)

    ids = [sample for sample in cohort["genome_id"].astype(str) if sample in targeted.columns and sample in unitigs.columns and sample in distance.index]
    cohort = cohort.set_index("genome_id").loc[ids].reset_index()
    distance = distance.loc[ids, ids]
    if len(ids) < 100:
        raise ValueError(f"Too few exact-ID external genomes after alignment: {len(ids)}")
    if set(ids) != set(distance.index) or set(distance.index) != set(distance.columns):
        raise ValueError("External distance matrix IDs differ")

    vectors: list[np.ndarray] = []
    mapping_rows: list[dict] = []
    missing_candidates: list[dict] = []
    for _, row in panel.iterrows():
        candidate_id = str(row["candidate_id"])
        key = str(row["matrix_key"])
        kind = str(row["candidate_type"])
        if kind == "targeted_variant":
            if key not in targeted.index:
                missing_candidates.append(row.to_dict())
                continue
            vector = targeted.loc[key, ids].to_numpy(np.uint8)
        elif kind == "unitig":
            if candidate_id not in unitigs.index:
                missing_candidates.append(row.to_dict())
                continue
            vector = unitigs.loc[candidate_id, ids].to_numpy(np.uint8)
        else:
            raise ValueError(f"Unknown candidate type {kind}")
        vectors.append(vector)
        mapping_rows.append({
            "candidate_id": candidate_id,
            "matrix_key": key,
            "candidate_type": kind,
            "known_mgrB_axis": bool(kind == "targeted_variant" and key.startswith("mgrB:")),
        })
    mapping = pd.DataFrame(mapping_rows)
    (out / "MISSING_FROZEN_CANDIDATES.json").write_text(
        json.dumps(missing_candidates, indent=2, ensure_ascii=False) + "\n"
    )
    if not vectors:
        raise ValueError("No frozen candidate could be represented")
    x = np.column_stack(vectors).astype(np.uint8)
    y = cohort["phenotype"].eq("R").astype(np.uint8).to_numpy()
    af = x.mean(axis=0)
    mapping["external_carriers"] = x.sum(axis=0).astype(int)
    mapping["external_af"] = af
    mapping["externally_analysable"] = (af >= args.min_af) & (af <= args.max_af)
    mapping.to_csv(out / "EXTERNAL_CANDIDATE_MATRIX_MAPPING.csv", index=False)

    eligible_indices = np.flatnonzero(mapping["externally_analysable"].to_numpy())
    if len(eligible_indices) == 0:
        summary = {
            "status": "NO_FROZEN_CANDIDATE_VARIABLE_IN_EXTERNAL_COHORT",
            "n_external_samples": len(cohort),
            "external_counts": cohort["phenotype"].value_counts().to_dict(),
            "n_frozen_candidates": len(panel),
            "n_represented_candidates": len(mapping),
            "n_analysable_candidates": 0,
            "n_strict_external_replicates": 0,
            "strict_candidates": [],
        }
        (out / "BVBRC_EXTERNAL_VALIDATION_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return

    model_tables: list[pd.DataFrame] = []
    for dim in dims:
        pcs = pcoa(distance, dim)
        covariates = np.column_stack([np.ones(len(ids)), pcs])
        stats = score_test_matrix(x[:, eligible_indices], y, covariates)
        stats["candidate_id"] = mapping.iloc[eligible_indices[stats["feature_index"].to_numpy()]]["candidate_id"].to_numpy()
        stats["q_score"] = bh(stats["p_score"])
        stats["pc_dim"] = dim
        stats.to_csv(out / f"BVBRC_STRUCTURE_ADJUSTED_PC{dim}.csv", index=False)
        model_tables.append(stats[["candidate_id", "pc_dim", "beta_score", "se_score", "z_score", "p_score", "q_score"]])
    long = pd.concat(model_tables, ignore_index=True)
    wide = long.pivot(index="candidate_id", columns="pc_dim")
    wide.columns = [f"{metric}_pc{dimension}" for metric, dimension in wide.columns]
    wide = wide.reset_index()
    beta_columns = [f"beta_score_pc{dimension}" for dimension in dims]
    p_columns = [f"p_score_pc{dimension}" for dimension in dims]
    q_columns = [f"q_score_pc{dimension}" for dimension in dims]
    wide["beta_min"] = wide[beta_columns].min(axis=1)
    wide["p_max"] = wide[p_columns].max(axis=1)
    wide["q_max"] = wide[q_columns].max(axis=1)
    wide["structure_adjusted_replication"] = (wide["beta_min"] > 0) & (wide["q_max"] <= args.alpha)

    source_groups = group_series(cohort, ["bioproject_accession", "source_group", "pmid_text"])
    countries = group_series(cohort, ["country", "geographic_location"])
    evidence_rows: list[dict] = []
    meta_details: dict[str, dict] = {}
    id_to_col = {candidate: index for index, candidate in enumerate(mapping["candidate_id"])}
    for candidate_id in mapping.loc[mapping["externally_analysable"], "candidate_id"]:
        index = id_to_col[candidate_id]
        vector = x[:, index]
        a, b, c, d = contingency(vector, y)
        odds_ratio, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
        fisher_p = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        source = random_effects_informative(vector, y, source_groups)
        country = random_effects_informative(vector, y, countries)
        meta_details[candidate_id] = {"source": source, "country": country}
        evidence_rows.append({
            "candidate_id": candidate_id,
            "R_present": a,
            "S_present": b,
            "R_absent": c,
            "S_absent": d,
            "external_carriers": int(a + b),
            "external_odds_ratio": odds_ratio,
            "external_ci_low": ci_low,
            "external_ci_high": ci_high,
            "external_fisher_p": fisher_p,
            "source_n_groups": source.get("n_groups", 0),
            "source_odds_ratio": source.get("odds_ratio"),
            "source_ci_low": source.get("ci_low"),
            "source_ci_high": source.get("ci_high"),
            "source_p": source.get("one_sided_p"),
            "source_I2": source.get("I2_percent"),
            "country_n_groups": country.get("n_groups", 0),
            "country_odds_ratio": country.get("odds_ratio"),
            "country_ci_low": country.get("ci_low"),
            "country_ci_high": country.get("ci_high"),
            "country_p": country.get("one_sided_p"),
            "country_I2": country.get("I2_percent"),
        })
    evidence = pd.DataFrame(evidence_rows)
    evidence["external_fisher_q"] = bh(evidence["external_fisher_p"])
    evidence = evidence.merge(wide, on="candidate_id", how="left")
    evidence = evidence.merge(mapping, on="candidate_id", how="left")
    evidence["unadjusted_external_replication"] = (
        (evidence["external_carriers_x"] >= args.min_carriers)
        & (evidence["external_odds_ratio"] >= args.min_odds_ratio)
        & (evidence["external_ci_low"] > 1)
        & (evidence["external_fisher_q"] <= args.alpha)
    )
    evidence["source_replication"] = (
        (evidence["source_n_groups"] >= args.min_groups)
        & (pd.to_numeric(evidence["source_ci_low"], errors="coerce") > 1)
        & (pd.to_numeric(evidence["source_p"], errors="coerce") <= args.alpha)
    )
    evidence["country_replication"] = (
        (evidence["country_n_groups"] >= args.min_groups)
        & (pd.to_numeric(evidence["country_ci_low"], errors="coerce") > 1)
        & (pd.to_numeric(evidence["country_p"], errors="coerce") <= args.alpha)
    )
    evidence["strict_external_replication"] = (
        evidence["unadjusted_external_replication"]
        & evidence["structure_adjusted_replication"]
        & evidence["source_replication"]
        & evidence["country_replication"]
    )
    evidence["potentially_novel_external_replication"] = (
        evidence["strict_external_replication"] & ~evidence["known_mgrB_axis"]
    )
    evidence = evidence.sort_values(
        ["potentially_novel_external_replication", "strict_external_replication", "q_max", "external_fisher_q"],
        ascending=[False, False, True, True],
    )
    evidence.to_csv(out / "BVBRC_ALL_FROZEN_CANDIDATE_EVIDENCE.csv", index=False)
    strict = evidence[evidence["strict_external_replication"]].copy()
    potential = evidence[evidence["potentially_novel_external_replication"]].copy()
    strict.to_csv(out / "BVBRC_STRICT_EXTERNAL_REPLICATES.csv", index=False)
    potential.to_csv(out / "BVBRC_POTENTIALLY_NOVEL_EXTERNAL_REPLICATES.csv", index=False)
    (out / "BVBRC_SOURCE_COUNTRY_META_DETAILS.json").write_text(
        json.dumps(meta_details, indent=2, ensure_ascii=False, default=str) + "\n"
    )

    summary = {
        "status": (
            "NON_MGRB_EXTERNAL_REPLICATES_REQUIRE_CONTEXT_DATABASE_LITERATURE_AUDIT"
            if len(potential)
            else "NO_PREVIOUSLY_UNRECOGNIZED_MARKER_SURVIVED_COMPLETE_BVBRC_GATE"
        ),
        "n_external_samples": int(len(cohort)),
        "external_counts": cohort["phenotype"].value_counts().to_dict(),
        "n_external_bioprojects_informative": int(source_groups.replace("", np.nan).nunique(dropna=True)),
        "n_external_countries_informative": int(countries.replace("", np.nan).nunique(dropna=True)),
        "n_frozen_candidates": int(len(panel)),
        "n_represented_candidates": int(len(mapping)),
        "n_missing_candidates": int(len(missing_candidates)),
        "n_analysable_candidates": int(mapping["externally_analysable"].sum()),
        "n_unadjusted_external_replicates": int(evidence["unadjusted_external_replication"].sum()),
        "n_structure_adjusted_replicates": int(evidence["structure_adjusted_replication"].sum()),
        "n_source_replicates": int(evidence["source_replication"].sum()),
        "n_country_replicates": int(evidence["country_replication"].sum()),
        "n_strict_external_replicates": int(len(strict)),
        "strict_candidates": strict["candidate_id"].astype(str).tolist(),
        "n_potentially_novel_external_replicates": int(len(potential)),
        "potentially_novel_candidates": potential["candidate_id"].astype(str).tolist(),
        "known_mgrB_positive_controls": strict.loc[strict["known_mgrB_axis"], "candidate_id"].astype(str).tolist(),
        "boundary": (
            "Only laboratory-method BV-BRC records were tested. Statistical replication is not causality, "
            "clinical validity, or proof of novelty. Non-mgrB survivors require exact genomic context, "
            "known-mechanism, database, literature, and independent review; laboratory experiments are "
            "required to establish a resistance mechanism."
        ),
    }
    (out / "BVBRC_EXTERNAL_VALIDATION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    report = [
        "# Laboratory-only BV-BRC external validation of frozen candidates",
        "",
        f"- External genomes: **{len(cohort):,}** {summary['external_counts']}",
        f"- Frozen candidates: **{len(panel):,}**",
        f"- Represented / analysable: **{len(mapping):,} / {summary['n_analysable_candidates']:,}**",
        f"- Structure-adjusted replicates: **{summary['n_structure_adjusted_replicates']:,}**",
        f"- Complete external gate: **{len(strict):,}**",
        f"- Non-mgrB complete external gate: **{len(potential):,}**",
        "",
        "## Claim boundary",
        "",
        summary["boundary"],
    ]
    (out / "BVBRC_EXTERNAL_VALIDATION_REPORT.md").write_text("\n".join(report) + "\n")
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
