#!/usr/bin/env python3
"""Phenotype-blind lineage-stratified CMH discovery and validation (v2).

This version fixes a critical discovery/validation column-collision bug in the first
implementation. Discovery and validation statistics are given explicit, non-overlapping
column names before merging. Candidate selection uses discovery samples only; the
untouched validation split is queried only after the discovery candidate set is frozen.

A replicated association is not evidence of a novel resistance mechanism. Known-mechanism,
sequence-context, database/literature, clone/source, and biological-validation gates remain
required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.cluster import KMeans


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--feature-meta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--clusters", default="8,12,16")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-af", type=float, default=0.01)
    p.add_argument("--max-af", type=float, default=0.99)
    p.add_argument("--min-stratum-class", type=int, default=3)
    return p.parse_args()


def bh(values: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    if len(v):
        order = np.argsort(v)
        ranked = v[order]
        q = np.minimum.accumulate(
            (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
        )[::-1]
        out[np.flatnonzero(ok)[order]] = np.clip(q, 0.0, 1.0)
    return out


def pcoa(distance: pd.DataFrame, n_components: int = 40) -> np.ndarray:
    x = distance.to_numpy(dtype=float)
    x = (x + x.T) / 2.0
    np.fill_diagonal(x, 0.0)
    n = len(x)
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (x * x) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    keep = eigenvalues > max(eigenvalues[0] * 1e-12, 1e-15)
    eigenvalues = eigenvalues[keep]
    eigenvectors = eigenvectors[:, keep]
    k = min(n_components, len(eigenvalues))
    coordinates = eigenvectors[:, :k] * np.sqrt(eigenvalues[:k])
    scale = coordinates.std(axis=0)
    scale[scale == 0] = 1.0
    return (coordinates - coordinates.mean(axis=0)) / scale


def cmh_for_feature(
    feature: np.ndarray,
    phenotype: np.ndarray,
    cluster: np.ndarray,
    min_class: int,
) -> dict:
    numerator = 0.0
    variance_sum = 0.0
    or_numerator = 0.0
    or_denominator = 0.0
    informative = 0
    positive = 0
    details: list[dict] = []

    for group in sorted(set(cluster)):
        idx = cluster == group
        y = phenotype[idx]
        x = feature[idx]
        resistant = int(y.sum())
        susceptible = int(len(y) - resistant)
        if resistant < min_class or susceptible < min_class:
            continue

        a = int(((x == 1) & (y == 1)).sum())
        b = int(((x == 1) & (y == 0)).sum())
        c = int(((x == 0) & (y == 1)).sum())
        d = int(((x == 0) & (y == 0)).sum())
        n = a + b + c + d
        if n <= 1 or (a + b) == 0 or (c + d) == 0:
            continue

        expected = (a + b) * (a + c) / n
        variance = (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        if variance <= 0:
            continue

        numerator += a - expected
        variance_sum += variance
        or_numerator += a * d / n
        or_denominator += b * c / n
        informative += 1

        aa, bb, cc, dd = map(float, (a, b, c, d))
        if min(aa, bb, cc, dd) == 0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5
        log_or = math.log((aa * dd) / (bb * cc))
        positive += int(log_or > 0)
        details.append(
            {
                "cluster": int(group),
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "log_or": log_or,
            }
        )

    if informative == 0 or variance_sum <= 0:
        return {
            "p": np.nan,
            "z": np.nan,
            "or": np.nan,
            "n_strata": 0,
            "positive_strata": 0,
            "details": details,
        }

    z = numerator / math.sqrt(variance_sum)
    p = float(norm.sf(z))
    common_or = (or_numerator / or_denominator) if or_denominator > 0 else float("inf")
    return {
        "p": p,
        "z": z,
        "or": common_or,
        "n_strata": informative,
        "positive_strata": positive,
        "details": details,
    }


def run_matrix(
    matrix: np.ndarray,
    phenotype: np.ndarray,
    clusters: np.ndarray,
    features: np.ndarray,
    min_class: int,
) -> tuple[pd.DataFrame, dict]:
    rows: list[dict] = []
    details: dict = {}
    for index, feature_name in enumerate(features):
        result = cmh_for_feature(matrix[:, index], phenotype, clusters, min_class)
        details[str(feature_name)] = result.pop("details")
        result["feature"] = str(feature_name)
        rows.append(result)
    table = pd.DataFrame(rows)
    table["q"] = bh(table["p"])
    return table, details


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cluster_counts = [int(value) for value in args.clusters.split(",")]

    manifest = pd.read_csv(args.manifest, dtype={"assembly_ID": str}).drop_duplicates("assembly_ID")
    manifest["y"] = manifest["phenotype"].astype(str).eq("R").astype(int)
    sample_ids = manifest["assembly_ID"].astype(str).tolist()

    rtab = pd.read_csv(args.rtab, sep="\t", index_col=0)
    rtab.index = rtab.index.astype(str)
    rtab.columns = rtab.columns.astype(str)
    rtab = rtab.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)

    distance = pd.read_csv(args.distance, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)

    if (
        set(sample_ids) != set(rtab.columns)
        or set(sample_ids) != set(distance.index)
        or set(distance.index) != set(distance.columns)
    ):
        raise SystemExit("Sample ID mismatch")

    rtab = rtab.loc[:, sample_ids]
    distance = distance.loc[sample_ids, sample_ids]
    coordinates = pcoa(distance, 40)
    matrix = rtab.T.to_numpy(dtype=np.uint8)
    phenotype = manifest["y"].to_numpy(dtype=np.uint8)
    discovery = manifest["split"].eq("discovery").to_numpy()
    validation = manifest["split"].eq("validation").to_numpy()
    if not discovery.any() or not validation.any():
        raise SystemExit("Empty discovery or validation split")

    discovery_af = matrix[discovery].mean(axis=0)
    eligible = (discovery_af >= args.min_af) & (discovery_af <= args.max_af)
    features = rtab.index.to_numpy()[eligible]
    discovery_matrix = matrix[discovery][:, eligible]
    validation_matrix = matrix[validation][:, eligible]
    discovery_y = phenotype[discovery]
    validation_y = phenotype[validation]

    discovery_tables: list[pd.DataFrame] = []
    validation_tables: list[pd.DataFrame] = []
    all_details: dict = {}

    for k in cluster_counts:
        labels = KMeans(n_clusters=k, random_state=20260819, n_init=50).fit_predict(coordinates)
        manifest[f"genetic_cluster_k{k}"] = labels

        discovery_table, discovery_details = run_matrix(
            discovery_matrix,
            discovery_y,
            labels[discovery],
            features,
            args.min_stratum_class,
        )
        validation_table, validation_details = run_matrix(
            validation_matrix,
            validation_y,
            labels[validation],
            features,
            args.min_stratum_class,
        )

        discovery_table = discovery_table.rename(
            columns={
                "p": f"p_discovery_k{k}",
                "z": f"z_discovery_k{k}",
                "or": f"or_discovery_k{k}",
                "n_strata": f"n_strata_discovery_k{k}",
                "positive_strata": f"positive_strata_discovery_k{k}",
                "q": f"q_discovery_k{k}",
            }
        )
        validation_table = validation_table.rename(
            columns={
                "p": f"p_validation_k{k}",
                "z": f"z_validation_k{k}",
                "or": f"or_validation_k{k}",
                "n_strata": f"n_strata_validation_k{k}",
                "positive_strata": f"positive_strata_validation_k{k}",
                "q": f"q_validation_within_all_features_k{k}",
            }
        )
        discovery_table.to_csv(out / f"DISCOVERY_CMH_K{k}.csv", index=False)
        validation_table.to_csv(out / f"VALIDATION_CMH_K{k}.csv", index=False)
        discovery_tables.append(discovery_table)
        validation_tables.append(validation_table)
        all_details[f"k{k}"] = {
            "discovery": discovery_details,
            "validation": validation_details,
        }

    discovery_merged = discovery_tables[0]
    for table in discovery_tables[1:]:
        discovery_merged = discovery_merged.merge(table, on="feature", how="outer")

    discovery_p = [f"p_discovery_k{k}" for k in cluster_counts]
    discovery_q = [f"q_discovery_k{k}" for k in cluster_counts]
    discovery_or = [f"or_discovery_k{k}" for k in cluster_counts]
    discovery_n = [f"n_strata_discovery_k{k}" for k in cluster_counts]
    discovery_complete = discovery_merged[discovery_p + discovery_or + discovery_n].notna().all(axis=1)
    discovery_merged["p_max_discovery"] = discovery_merged[discovery_p].max(axis=1)
    discovery_merged["q_max_discovery"] = discovery_merged[discovery_q].max(axis=1)
    discovery_merged["or_min_discovery"] = discovery_merged[discovery_or].min(axis=1)
    discovery_merged["min_informative_strata_discovery"] = discovery_merged[discovery_n].min(axis=1)
    discovery_merged["discovery_model_complete"] = discovery_complete
    discovery_merged["discovery_stable"] = (
        discovery_complete
        & (discovery_merged["q_max_discovery"] <= args.alpha)
        & (discovery_merged["or_min_discovery"] > 1)
        & (discovery_merged["min_informative_strata_discovery"] >= 2)
    )
    frozen = discovery_merged[discovery_merged["discovery_stable"]].copy()
    frozen = frozen.sort_values(["q_max_discovery", "p_max_discovery"])
    frozen.to_csv(out / "FROZEN_LINEAGE_STRATIFIED_CANDIDATES.csv", index=False)

    validation_merged = validation_tables[0]
    for table in validation_tables[1:]:
        validation_merged = validation_merged.merge(table, on="feature", how="outer")
    evidence = frozen.merge(validation_merged, on="feature", how="left", validate="one_to_one")

    validation_p = [f"p_validation_k{k}" for k in cluster_counts]
    validation_or = [f"or_validation_k{k}" for k in cluster_counts]
    validation_n = [f"n_strata_validation_k{k}" for k in cluster_counts]
    validation_complete = evidence[validation_p + validation_or + validation_n].notna().all(axis=1)

    if len(evidence):
        evidence["validation_p_max"] = evidence[validation_p].max(axis=1)
        evidence["validation_q_across_frozen_candidates"] = bh(evidence["validation_p_max"])
        evidence["validation_or_min"] = evidence[validation_or].min(axis=1)
        evidence["validation_min_informative_strata"] = evidence[validation_n].min(axis=1)
        evidence["validation_model_complete"] = validation_complete
        evidence["heldout_lineage_replication"] = (
            validation_complete
            & (evidence["validation_q_across_frozen_candidates"] <= args.alpha)
            & (evidence["validation_or_min"] > 1)
            & (evidence["validation_min_informative_strata"] >= 2)
        )
    else:
        evidence["validation_p_max"] = pd.Series(dtype=float)
        evidence["validation_q_across_frozen_candidates"] = pd.Series(dtype=float)
        evidence["validation_or_min"] = pd.Series(dtype=float)
        evidence["validation_min_informative_strata"] = pd.Series(dtype=float)
        evidence["validation_model_complete"] = pd.Series(dtype=bool)
        evidence["heldout_lineage_replication"] = pd.Series(dtype=bool)

    evidence.to_csv(out / "LINEAGE_STRATIFIED_DISCOVERY_VALIDATION_EVIDENCE.csv", index=False)
    strict = evidence[evidence["heldout_lineage_replication"].eq(True)].copy()
    strict.to_csv(out / "STRICT_LINEAGE_STRATIFIED_REPLICATES.csv", index=False)
    manifest.to_csv(out / "PHENOTYPE_BLIND_GENETIC_CLUSTERS.csv", index=False)
    (out / "CMH_DETAILS.json").write_text(
        json.dumps(all_details, ensure_ascii=False, default=str) + "\n"
    )

    summary = {
        "analysis_version": "v2_explicit_discovery_validation_columns",
        "critical_fix": "Validation gates use only p_validation/or_validation/n_strata_validation columns. The v1 collision with discovery columns is removed.",
        "n_all": int(len(manifest)),
        "n_discovery": int(discovery.sum()),
        "n_validation": int(validation.sum()),
        "n_eligible_features": int(len(features)),
        "cluster_counts": cluster_counts,
        "n_discovery_stable": int(len(frozen)),
        "n_heldout_lineage_replicated": int(len(strict)),
        "strict_features": strict["feature"].astype(str).tolist() if len(strict) else [],
        "status": (
            "LINEAGE_STRATIFIED_CANDIDATES_REQUIRE_KNOWN_MECHANISM_AND_UNITIG_CONCORDANCE"
            if len(strict)
            else "NO_FEATURE_SURVIVED_LINEAGE_STRATIFIED_DISCOVERY_VALIDATION"
        ),
        "boundary": "Within-lineage statistical replication does not establish novelty or causality. Known mechanisms, sequence context, independent tools, literature and biological validation remain required.",
    }
    (out / "LINEAGE_STRATIFIED_CMH_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )

    report = [
        "# Lineage-stratified CMH sensitivity analysis (corrected v2)",
        "",
        "The v1 evaluator accidentally reused discovery columns when computing the validation gate. This version uses explicit `*_discovery_*` and `*_validation_*` columns and requires complete results across all requested cluster counts.",
        "",
        f"- Eligible features: **{len(features):,}**",
        f"- Discovery-stable features: **{len(frozen):,}**",
        f"- Untouched validation replicates: **{len(strict):,}**",
        "",
        summary["boundary"],
    ]
    if len(evidence):
        report += ["", "## Evidence-ranked candidates", "", evidence.head(50).to_markdown(index=False)]
    (out / "LINEAGE_STRATIFIED_CMH_REPORT.md").write_text("\n".join(report) + "\n")

    hashes = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}"
        for path in sorted(out.rglob("*"))
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
