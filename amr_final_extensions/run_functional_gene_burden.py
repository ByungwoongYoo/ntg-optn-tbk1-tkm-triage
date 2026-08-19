#!/usr/bin/env python3
"""Test phenotype-blind functional gene burdens with frozen BioProject validation.

Reference-relative variants are grouped before outcomes are inspected into prespecified
burden classes (missense, protein indel, severe disruption, promoter, and combined).
Discovery uses only the frozen discovery BioProjects. Candidate burdens must remain
positive across 10/20/30 Mash-PC models, replicate in untouched validation, and pass
BioProject- and country-stratified random-effects gates. A surviving burden is an
association, not a causal or novel resistance mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from amr_discovery.run_structure_adjusted_score_gwas import (
    bh,
    contingency,
    odds_ratio_ci,
    pcoa,
    random_effects,
    score_test_matrix,
)

BURDEN_CLASSES: dict[str, set[str]] = {
    "ANY_NONREFERENCE": {"gene_burden_any_nonreference"},
    "MISSENSE_ANY": {"aa_substitution"},
    "PROTEIN_INDEL": {"aa_insertion", "aa_deletion"},
    "SEVERE_DISRUPTION": {
        "potential_loss_of_function",
        "partial_or_divergent",
        "gene_absent_or_undetected",
        "premature_stop",
    },
    "PROMOTER_ANY": {
        "promoter_substitution",
        "promoter_deletion",
        "promoter_partial_or_structural",
        "promoter_undetected",
    },
    "PROMOTER_STRUCTURAL": {
        "promoter_deletion",
        "promoter_partial_or_structural",
        "promoter_undetected",
    },
    "ALL_FUNCTIONAL": {
        "aa_substitution",
        "aa_insertion",
        "aa_deletion",
        "potential_loss_of_function",
        "partial_or_divergent",
        "gene_absent_or_undetected",
        "premature_stop",
        "promoter_substitution",
        "promoter_deletion",
        "promoter_partial_or_structural",
        "promoter_undetected",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--feature-meta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pc-dims", default="10,20,30")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5)
    p.add_argument("--min-validation-odds-ratio", type=float, default=2.0)
    return p.parse_args()


def load_matrix(path: str | Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)


def construct_burdens(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame]:
    metadata = metadata.drop_duplicates("feature").set_index("feature")
    definitions: list[dict] = []
    vectors: list[np.ndarray] = []
    for gene in sorted(metadata["gene"].dropna().astype(str).unique()):
        gene_meta = metadata[metadata["gene"].astype(str).eq(gene)]
        for burden_class, feature_types in BURDEN_CLASSES.items():
            features = gene_meta[gene_meta["type"].isin(feature_types)].index.intersection(matrix.index)
            if len(features) == 0:
                continue
            vector = (matrix.loc[features].sum(axis=0) > 0).astype(np.uint8).to_numpy()
            if int(vector.sum()) in {0, len(vector)}:
                continue
            definitions.append({
                "burden_id": f"{gene}|{burden_class}",
                "gene": gene,
                "burden_class": burden_class,
                "n_component_features": int(len(features)),
                "component_features": ";".join(map(str, features)),
                "all_carriers": int(vector.sum()),
                "all_frequency": float(vector.mean()),
            })
            vectors.append(vector)
    if not vectors:
        return np.empty((matrix.shape[1], 0), dtype=np.uint8), pd.DataFrame()
    return np.column_stack(vectors).astype(np.uint8), pd.DataFrame(definitions)


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    dims = [int(value) for value in args.pc_dims.split(",") if value.strip()]

    manifest = pd.read_csv(args.manifest, dtype={"assembly_ID": str}).drop_duplicates("assembly_ID")
    required_manifest = {"assembly_ID", "phenotype", "split"}
    if missing := required_manifest - set(manifest.columns):
        raise ValueError(f"Manifest missing {sorted(missing)}")
    manifest["y"] = manifest["phenotype"].astype(str).eq("R").astype(np.uint8)
    ids = manifest["assembly_ID"].astype(str).tolist()

    matrix = load_matrix(args.rtab)
    distance = pd.read_csv(args.distance, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)
    if set(ids) != set(matrix.columns) or set(ids) != set(distance.index) or set(distance.index) != set(distance.columns):
        raise ValueError("Manifest, feature matrix, and Mash distance IDs differ")
    matrix = matrix.loc[:, ids]
    distance = distance.loc[ids, ids]
    metadata = pd.read_csv(args.feature_meta, dtype=str)
    for column in ["feature", "gene", "type"]:
        if column not in metadata.columns:
            raise ValueError(f"Feature metadata missing {column}")

    burden_matrix, definitions = construct_burdens(matrix, metadata)
    if definitions.empty:
        summary = {
            "status": "NO_CONSTRUCTABLE_FUNCTIONAL_BURDEN",
            "n_strict_statistically_replicated": 0,
            "boundary": "No burden met the phenotype-blind construction rules; no marker is claimed.",
        }
        (out / "FUNCTIONAL_GENE_BURDEN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return
    definitions.to_csv(out / "FUNCTIONAL_BURDEN_DEFINITIONS.csv", index=False)

    discovery_mask = manifest["split"].astype(str).eq("discovery").to_numpy()
    validation_mask = manifest["split"].astype(str).eq("validation").to_numpy()
    if not discovery_mask.any() or not validation_mask.any():
        raise ValueError("Empty discovery or validation split")
    if set(manifest.loc[discovery_mask, "assembly_ID"]) & set(manifest.loc[validation_mask, "assembly_ID"]):
        raise ValueError("Discovery and validation overlap")

    y = manifest["y"].to_numpy(np.uint8)
    discovery_ids = manifest.loc[discovery_mask, "assembly_ID"].astype(str).tolist()
    discovery_distance = distance.loc[discovery_ids, discovery_ids]
    model_tables: list[pd.DataFrame] = []
    for dimension in dims:
        coordinates = pcoa(discovery_distance, dimension)
        covariates = np.column_stack([np.ones(len(discovery_ids)), coordinates])
        stats = score_test_matrix(burden_matrix[discovery_mask], y[discovery_mask], covariates)
        stats["burden_id"] = definitions.iloc[stats["feature_index"].to_numpy()]["burden_id"].to_numpy()
        stats["q_score"] = bh(stats["p_score"])
        stats["pc_dim"] = dimension
        stats.to_csv(out / f"DISCOVERY_FUNCTIONAL_BURDEN_PC{dimension}.csv", index=False)
        model_tables.append(stats[["burden_id", "pc_dim", "beta_score", "p_score", "q_score"]])

    long = pd.concat(model_tables, ignore_index=True)
    wide = long.pivot(index="burden_id", columns="pc_dim")
    wide.columns = [f"{metric}_pc{dimension}" for metric, dimension in wide.columns]
    wide = wide.reset_index()
    beta_columns = [f"beta_score_pc{dimension}" for dimension in dims]
    p_columns = [f"p_score_pc{dimension}" for dimension in dims]
    q_columns = [f"q_score_pc{dimension}" for dimension in dims]
    wide["beta_min"] = wide[beta_columns].min(axis=1)
    wide["p_max"] = wide[p_columns].max(axis=1)
    wide["q_max"] = wide[q_columns].max(axis=1)
    wide["discovery_stable"] = (wide["beta_min"] > 0) & (wide["q_max"] <= args.alpha)
    frozen = wide[wide["discovery_stable"]].merge(definitions, on="burden_id", how="left")
    frozen = frozen.sort_values(["q_max", "p_max"])
    frozen.to_csv(out / "FROZEN_DISCOVERY_FUNCTIONAL_BURDENS.csv", index=False)

    validation_manifest = manifest.loc[validation_mask].reset_index(drop=True)
    validation_y = y[validation_mask]
    burden_index = {burden_id: index for index, burden_id in enumerate(definitions["burden_id"])}
    evidence_rows: list[dict] = []
    details: dict[str, dict] = {}
    for _, candidate in frozen.iterrows():
        burden_id = str(candidate["burden_id"])
        index = burden_index[burden_id]
        discovery_vector = burden_matrix[discovery_mask, index]
        validation_vector = burden_matrix[validation_mask, index]
        da, db, dc, dd = contingency(discovery_vector, y[discovery_mask])
        va, vb, vc, vd = contingency(validation_vector, validation_y)
        discovery_or, discovery_low, discovery_high = odds_ratio_ci(da, db, dc, dd)
        validation_or, validation_low, validation_high = odds_ratio_ci(va, vb, vc, vd)
        validation_p = float(fisher_exact([[va, vb], [vc, vd]], alternative="greater").pvalue)
        source_group = validation_manifest.get(
            "BioProject",
            validation_manifest.get("source_group", pd.Series(["UNKNOWN"] * len(validation_manifest))),
        )
        source = random_effects(validation_vector, validation_y, source_group)
        country = random_effects(
            validation_vector,
            validation_y,
            validation_manifest.get("ISO_country_code", pd.Series(["UNKNOWN"] * len(validation_manifest))),
        )
        details[burden_id] = {"source": source, "country": country}
        row = candidate.to_dict()
        row.update({
            "discovery_R_present": da,
            "discovery_S_present": db,
            "discovery_R_absent": dc,
            "discovery_S_absent": dd,
            "discovery_or": discovery_or,
            "discovery_ci_low": discovery_low,
            "discovery_ci_high": discovery_high,
            "validation_R_present": va,
            "validation_S_present": vb,
            "validation_R_absent": vc,
            "validation_S_absent": vd,
            "validation_or": validation_or,
            "validation_ci_low": validation_low,
            "validation_ci_high": validation_high,
            "validation_p": validation_p,
            "source_n_groups": source.get("n_groups", 0),
            "source_or": source.get("odds_ratio"),
            "source_ci_low": source.get("ci_low"),
            "source_p": source.get("one_sided_p"),
            "country_n_groups": country.get("n_groups", 0),
            "country_or": country.get("odds_ratio"),
            "country_ci_low": country.get("ci_low"),
            "country_p": country.get("one_sided_p"),
        })
        evidence_rows.append(row)

    evidence = pd.DataFrame(evidence_rows)
    if len(evidence):
        evidence["validation_q"] = bh(evidence["validation_p"])
        evidence["heldout_replication"] = (
            (evidence["validation_R_present"] + evidence["validation_S_present"] >= args.min_validation_present)
            & (evidence["validation_or"] >= args.min_validation_odds_ratio)
            & (evidence["validation_ci_low"] > 1)
            & (evidence["validation_q"] <= args.alpha)
        )
        evidence["source_replication"] = (
            (evidence["source_n_groups"] >= 3)
            & (pd.to_numeric(evidence["source_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["source_p"], errors="coerce") <= args.alpha)
        )
        evidence["country_replication"] = (
            (evidence["country_n_groups"] >= 3)
            & (pd.to_numeric(evidence["country_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["country_p"], errors="coerce") <= args.alpha)
        )
        evidence["strict_statistical_gate"] = (
            evidence["heldout_replication"]
            & evidence["source_replication"]
            & evidence["country_replication"]
        )
        evidence = evidence.sort_values(
            ["strict_statistical_gate", "validation_q", "q_max"],
            ascending=[False, True, True],
        )
    else:
        for column in [
            "validation_q", "heldout_replication", "source_replication",
            "country_replication", "strict_statistical_gate",
        ]:
            evidence[column] = pd.Series(dtype=float if column == "validation_q" else bool)
    evidence.to_csv(out / "FUNCTIONAL_BURDEN_DISCOVERY_VALIDATION_EVIDENCE.csv", index=False)
    strict = evidence[evidence["strict_statistical_gate"].eq(True)].copy() if len(evidence) else evidence.copy()
    strict.to_csv(out / "STRICT_REPLICATED_FUNCTIONAL_BURDENS.csv", index=False)
    (out / "FUNCTIONAL_BURDEN_META_ANALYSIS_DETAILS.json").write_text(
        json.dumps(details, indent=2, ensure_ascii=False, default=str) + "\n"
    )

    summary = {
        "n_all_samples": int(len(manifest)),
        "n_discovery": int(discovery_mask.sum()),
        "n_validation": int(validation_mask.sum()),
        "n_constructed_burdens": int(len(definitions)),
        "burden_class_counts": definitions["burden_class"].value_counts().to_dict(),
        "n_discovery_stable": int(len(frozen)),
        "n_heldout_replicated": int(evidence["heldout_replication"].sum()) if len(evidence) else 0,
        "n_source_replicated": int(evidence["source_replication"].sum()) if len(evidence) else 0,
        "n_country_replicated": int(evidence["country_replication"].sum()) if len(evidence) else 0,
        "n_strict_statistically_replicated": int(len(strict)),
        "strict_burdens": strict["burden_id"].astype(str).tolist() if len(strict) else [],
        "status": "FUNCTIONAL_BURDENS_REQUIRE_EXACT_SEQUENCE_AND_KNOWN_MECHANISM_AUDIT" if len(strict) else "NO_FUNCTIONAL_GENE_BURDEN_SURVIVED_COMPLETE_GATE",
        "boundary": "A replicated gene burden is not a novel mechanism. Component variants, assembly context, established mechanisms, external data, and causality must be audited independently.",
    }
    (out / "FUNCTIONAL_GENE_BURDEN_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    report = [
        "# Functional gene-burden validation",
        "",
        f"- Constructed phenotype-blind burdens: **{len(definitions):,}**",
        f"- Discovery-stable burdens: **{len(frozen):,}**",
        f"- Untouched validation replicates: **{int(evidence['heldout_replication'].sum()) if len(evidence) else 0:,}**",
        f"- Complete source/country gate: **{len(strict):,}**",
        "",
        summary["boundary"],
    ]
    if len(evidence):
        columns = [
            "burden_id", "n_component_features", "all_frequency", "beta_min", "q_max",
            "validation_or", "validation_ci_low", "validation_q", "source_n_groups",
            "source_or", "source_ci_low", "country_n_groups", "country_or",
            "country_ci_low", "strict_statistical_gate",
        ]
        report += ["", "## Evidence-ranked burdens", "", evidence.head(50)[columns].to_markdown(index=False)]
    (out / "FUNCTIONAL_GENE_BURDEN_REPORT.md").write_text("\n".join(report) + "\n")

    hashes = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}"
        for path in sorted(out.rglob("*"))
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
