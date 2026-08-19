#!/usr/bin/env python3
"""Rare-variant gene-burden GWAS with frozen discovery and untouched validation.

Variants are grouped by target gene using phenotype-free annotations. Discovery
carrier-frequency ceilings of 5% and 10% are evaluated as prespecified sensitivity
sets. Only discovery-defined nonsynonymous/disruptive variants contribute to each
burden in validation. This is an association screen, not a causal mechanism call.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm

from run_structure_adjusted_score_gwas import (
    bh,
    contingency,
    odds_ratio_ci,
    pcoa,
    random_effects,
    score_test_matrix,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--feature-meta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--frequency-ceilings", default="0.05,0.10")
    p.add_argument("--pc-dims", default="10,20,30")
    p.add_argument("--min-carriers", type=int, default=3)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5)
    return p.parse_args()


def infer_gene(feature: str, metadata: pd.Series | None) -> str:
    if metadata is not None:
        for column in ["gene", "gene_name", "target_gene", "locus", "locus_tag"]:
            if column in metadata.index:
                value = str(metadata[column]).strip()
                if value and value.lower() not in {"nan", "none", "na", "unknown"}:
                    return value
    token = str(feature).split(":", 1)[0].split("|", 1)[0].strip()
    return token or "UNKNOWN"


def is_nonsyn_or_disruptive(feature: str, metadata: pd.Series | None) -> bool:
    text = str(feature)
    if metadata is not None:
        text += " " + " ".join(str(v) for v in metadata.tolist())
    upper = text.upper()
    if any(term in upper for term in [
        "STOP", "NONSENSE", "FRAMESHIFT", "FRAME_SHIFT", "TRUNC", "DISRUPT",
        "INSERT", "DELETION", "INDEL", "MISSING", "PARTIAL", "PREMATURE",
        "LOSS_OF_FUNCTION", "LOF", "IS_ELEMENT", "INTERPROTEIN",
    ]):
        return True
    # Typical amino-acid notation, e.g. N354A. Same first/last residue is synonymous.
    matches = re.findall(r"(?<![A-Z])([A-Z*])([0-9]+)([A-Z*])(?![A-Z])", upper)
    if matches:
        return any(first != last for first, _, last in matches)
    if any(term in upper for term in ["SYNONYMOUS", "SILENT"]):
        return False
    # Gene absence/coverage/coding disruption features without an AA token are retained.
    if any(term in upper for term in ["ABSENT", "COVERAGE", "BROKEN", "INCOMPLETE"]):
        return True
    return False


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ceilings = [float(x) for x in a.frequency_ceilings.split(",") if x.strip()]
    dims = [int(x) for x in a.pc_dims.split(",") if x.strip()]

    manifest = pd.read_csv(a.manifest, dtype={"assembly_ID": str}).drop_duplicates("assembly_ID")
    manifest["y"] = manifest["phenotype"].astype(str).eq("R").astype(int)
    ids = manifest["assembly_ID"].astype(str).tolist()

    rtab = pd.read_csv(a.rtab, sep="\t", index_col=0)
    rtab.index = rtab.index.astype(str)
    rtab.columns = rtab.columns.astype(str)
    rtab = rtab.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)

    distance = pd.read_csv(a.distance, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)
    if set(ids) != set(rtab.columns) or set(ids) != set(distance.index) or set(distance.index) != set(distance.columns):
        raise SystemExit("Manifest, feature matrix and distance matrix sample IDs differ")
    rtab = rtab.loc[:, ids]
    distance = distance.loc[ids, ids]

    meta = pd.read_csv(a.feature_meta, dtype=str)
    key = next((c for c in ["variant", "feature", "feature_id", "name"] if c in meta.columns), None)
    if key is not None:
        meta = meta.drop_duplicates(key).set_index(key)
    else:
        meta = pd.DataFrame(index=rtab.index)

    discovery_mask = manifest["split"].eq("discovery").to_numpy()
    validation_mask = manifest["split"].eq("validation").to_numpy()
    if not discovery_mask.any() or not validation_mask.any():
        raise SystemExit("Empty discovery or validation split")
    if set(manifest.loc[discovery_mask, "assembly_ID"]) & set(manifest.loc[validation_mask, "assembly_ID"]):
        raise SystemExit("Discovery and validation overlap")

    x_all = rtab.T.to_numpy(dtype=np.uint8)
    y_all = manifest["y"].to_numpy(dtype=np.uint8)
    xd = x_all[discovery_mask]
    yd = y_all[discovery_mask]
    xv = x_all[validation_mask]
    yv = y_all[validation_mask]

    annotations = []
    for feature in rtab.index:
        row = meta.loc[feature] if feature in meta.index else None
        annotations.append({
            "feature": feature,
            "gene": infer_gene(feature, row),
            "eligible_functional_class": is_nonsyn_or_disruptive(feature, row),
        })
    annotations = pd.DataFrame(annotations)
    annotations["discovery_carriers"] = xd.sum(axis=0)
    annotations["discovery_af"] = annotations["discovery_carriers"] / len(xd)
    annotations.to_csv(out / "FEATURE_BURDEN_ANNOTATION_AUDIT.csv", index=False)

    burdens = []
    burden_meta = []
    for ceiling in ceilings:
        eligible = annotations[
            annotations["eligible_functional_class"]
            & (annotations["discovery_carriers"] >= a.min_carriers)
            & (annotations["discovery_af"] <= ceiling)
        ]
        for gene, sub in eligible.groupby("gene"):
            feature_indices = [rtab.index.get_loc(f) for f in sub["feature"]]
            if not feature_indices:
                continue
            vector = (x_all[:, feature_indices].sum(axis=1) > 0).astype(np.uint8)
            burden_id = f"{gene}|rare_le_{ceiling:.2f}"
            burdens.append(vector)
            burden_meta.append({
                "burden_id": burden_id,
                "gene": gene,
                "frequency_ceiling": ceiling,
                "n_component_variants": len(feature_indices),
                "component_features": ";".join(sub["feature"].astype(str)),
                "discovery_burden_carriers": int(vector[discovery_mask].sum()),
                "validation_burden_carriers": int(vector[validation_mask].sum()),
            })
    if not burdens:
        summary = {
            "status": "NO_ELIGIBLE_RARE_VARIANT_BURDENS",
            "n_strict_statistically_replicated": 0,
            "claim_boundary": "No burden met the prespecified construction rules; no marker is claimed.",
        }
        (out / "RARE_VARIANT_BURDEN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return

    burden_matrix = np.column_stack(burdens).astype(np.uint8)
    burden_meta = pd.DataFrame(burden_meta)
    burden_meta.to_csv(out / "BURDEN_DEFINITIONS.csv", index=False)

    discovery_ids = manifest.loc[discovery_mask, "assembly_ID"].astype(str).tolist()
    discovery_distance = distance.loc[discovery_ids, discovery_ids]
    model_tables = []
    for dim in dims:
        coords = pcoa(discovery_distance, dim)
        covariates = np.column_stack([np.ones(len(discovery_ids)), coords])
        stats = score_test_matrix(burden_matrix[discovery_mask], yd, covariates)
        stats["burden_id"] = burden_meta.iloc[stats["feature_index"].to_numpy()]["burden_id"].to_numpy()
        stats["q_score"] = bh(stats["p_score"])
        stats["pc_dim"] = dim
        stats.to_csv(out / f"DISCOVERY_BURDEN_SCORE_PC{dim}.csv", index=False)
        model_tables.append(stats[["burden_id", "pc_dim", "beta_score", "p_score", "q_score"]])

    long = pd.concat(model_tables, ignore_index=True)
    wide = long.pivot(index="burden_id", columns="pc_dim")
    wide.columns = [f"{metric}_pc{dim}" for metric, dim in wide.columns]
    wide = wide.reset_index()
    beta_cols = [f"beta_score_pc{d}" for d in dims]
    p_cols = [f"p_score_pc{d}" for d in dims]
    q_cols = [f"q_score_pc{d}" for d in dims]
    wide["beta_min"] = wide[beta_cols].min(axis=1)
    wide["p_max"] = wide[p_cols].max(axis=1)
    wide["q_max"] = wide[q_cols].max(axis=1)
    wide["discovery_stable"] = (wide["beta_min"] > 0) & (wide["q_max"] <= a.alpha)
    frozen = wide[wide["discovery_stable"]].merge(burden_meta, on="burden_id", how="left")
    frozen = frozen.sort_values(["q_max", "p_max"])
    frozen.to_csv(out / "FROZEN_DISCOVERY_BURDENS.csv", index=False)

    val_manifest = manifest.loc[validation_mask].reset_index(drop=True)
    evidence_rows = []
    detail = {}
    burden_index = {bid: i for i, bid in enumerate(burden_meta["burden_id"])}
    for _, candidate in frozen.iterrows():
        bid = candidate["burden_id"]
        j = burden_index[bid]
        vd = burden_matrix[discovery_mask, j]
        vv = burden_matrix[validation_mask, j]
        va, vb, vc, vd0 = contingency(vv, yv)
        da, db, dc, dd = contingency(vd, yd)
        vor, vlo, vhi = odds_ratio_ci(va, vb, vc, vd0)
        dor, dlo, dhi = odds_ratio_ci(da, db, dc, dd)
        vp = float(fisher_exact([[va, vb], [vc, vd0]], alternative="greater").pvalue)
        source = random_effects(vv, yv, val_manifest.get("source_group", pd.Series(["UNKNOWN"] * len(val_manifest))))
        country = random_effects(vv, yv, val_manifest.get("ISO_country_code", pd.Series(["UNKNOWN"] * len(val_manifest))))
        detail[bid] = {"source": source, "country": country}
        row = candidate.to_dict()
        row.update({
            "discovery_R_present": da, "discovery_S_present": db, "discovery_R_absent": dc, "discovery_S_absent": dd,
            "discovery_or": dor, "discovery_ci_low": dlo, "discovery_ci_high": dhi,
            "validation_R_present": va, "validation_S_present": vb, "validation_R_absent": vc, "validation_S_absent": vd0,
            "validation_or": vor, "validation_ci_low": vlo, "validation_ci_high": vhi, "validation_p": vp,
            "source_n_groups": source.get("n_groups", 0), "source_or": source.get("odds_ratio"), "source_ci_low": source.get("ci_low"), "source_p": source.get("one_sided_p"),
            "country_n_groups": country.get("n_groups", 0), "country_or": country.get("odds_ratio"), "country_ci_low": country.get("ci_low"), "country_p": country.get("one_sided_p"),
        })
        evidence_rows.append(row)

    evidence = pd.DataFrame(evidence_rows)
    if len(evidence):
        evidence["validation_q"] = bh(evidence["validation_p"])
        evidence["heldout_replication"] = (
            (evidence["validation_R_present"] + evidence["validation_S_present"] >= a.min_validation_present)
            & (evidence["validation_or"] > 1)
            & (evidence["validation_ci_low"] > 1)
            & (evidence["validation_q"] <= a.alpha)
        )
        evidence["source_replication"] = (
            (evidence["source_n_groups"] >= 3)
            & (pd.to_numeric(evidence["source_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["source_p"], errors="coerce") <= a.alpha)
        )
        evidence["country_replication"] = (
            (evidence["country_n_groups"] >= 3)
            & (pd.to_numeric(evidence["country_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["country_p"], errors="coerce") <= a.alpha)
        )
        evidence["strict_statistical_gate"] = evidence["heldout_replication"] & evidence["source_replication"] & evidence["country_replication"]
        evidence = evidence.sort_values(["strict_statistical_gate", "validation_q", "q_max"], ascending=[False, True, True])
    else:
        for column in ["validation_q", "heldout_replication", "source_replication", "country_replication", "strict_statistical_gate"]:
            evidence[column] = pd.Series(dtype=float if column == "validation_q" else bool)
    evidence.to_csv(out / "RARE_BURDEN_DISCOVERY_VALIDATION_EVIDENCE.csv", index=False)
    strict = evidence[evidence["strict_statistical_gate"].eq(True)].copy() if len(evidence) else evidence.copy()
    strict.to_csv(out / "STRICT_REPLICATED_RARE_BURDENS.csv", index=False)
    (out / "BURDEN_META_ANALYSIS_DETAILS.json").write_text(json.dumps(detail, indent=2, ensure_ascii=False, default=str) + "\n")

    summary = {
        "n_all_samples": int(len(manifest)),
        "n_discovery": int(discovery_mask.sum()),
        "n_validation": int(validation_mask.sum()),
        "frequency_ceilings": ceilings,
        "n_constructed_burdens": int(len(burden_meta)),
        "n_discovery_stable": int(len(frozen)),
        "n_heldout_replicated": int(evidence["heldout_replication"].sum()) if len(evidence) else 0,
        "n_source_replicated": int(evidence["source_replication"].sum()) if len(evidence) else 0,
        "n_country_replicated": int(evidence["country_replication"].sum()) if len(evidence) else 0,
        "n_strict_statistically_replicated": int(len(strict)),
        "strict_burdens": strict["burden_id"].astype(str).tolist() if len(strict) else [],
        "status": "RARE_BURDEN_CANDIDATES_REQUIRE_KNOWN_MECHANISM_AND_CONTEXT_AUDIT" if len(strict) else "NO_RARE_VARIANT_BURDEN_SURVIVED_COMPLETE_GATE",
        "claim_boundary": "A replicated burden is not a novel resistance mechanism. Component variants and sequence context must be audited against known mechanisms and independent data, and causality still requires laboratory validation.",
    }
    (out / "RARE_VARIANT_BURDEN_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    report = [
        "# Rare-variant gene-burden colistin analysis", "",
        f"- Constructed burdens: **{summary['n_constructed_burdens']:,}**",
        f"- Discovery-stable burdens: **{summary['n_discovery_stable']:,}**",
        f"- Untouched validation replicates: **{summary['n_heldout_replicated']:,}**",
        f"- Complete statistical gate: **{summary['n_strict_statistically_replicated']:,}**",
        "", "## Claim boundary", "", summary["claim_boundary"],
    ]
    if len(evidence):
        columns = [c for c in ["burden_id", "n_component_variants", "q_max", "validation_or", "validation_ci_low", "validation_ci_high", "validation_q", "source_n_groups", "source_ci_low", "country_n_groups", "country_ci_low", "strict_statistical_gate"] if c in evidence.columns]
        report += ["", "## Evidence-ranked burdens", "", evidence.head(50)[columns].to_markdown(index=False)]
    (out / "RARE_VARIANT_BURDEN_REPORT.md").write_text("\n".join(report) + "\n")

    hashes = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}"
        for path in sorted(out.rglob("*")) if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
