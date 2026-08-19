#!/usr/bin/env python3
"""Prepare phenotype-blind, exact-ID panels for phylogenetic convergence analysis.

Targeted and selected-unitig features are filtered by discovery carrier frequency and
collapsed only when their full-cohort occurrence vectors are identical. Discovery and
validation remain BioProject-disjoint. No phenotype-dependent threshold is optimized in
this preparation step.
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
    p.add_argument("--extended-root", required=True)
    p.add_argument("--unitig-root", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-discovery-af", type=float, default=0.01)
    p.add_argument("--max-discovery-af", type=float, default=0.99)
    return p.parse_args()


def find_one(root: str | Path, name: str, contains: str | None = None) -> Path:
    hits = []
    for p in Path(root).rglob(name):
        normalized = str(p).replace("\\", "/")
        if contains is None or contains.replace("\\", "/") in normalized:
            hits.append(p)
    if not hits:
        raise FileNotFoundError(f"{name} not found under {root} contains={contains}")
    return sorted(hits, key=lambda p: (len(p.parts), str(p)))[0]


def load_rtab(path: Path) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = x.index.astype(str)
    x.columns = x.columns.astype(str)
    return x.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)


def write_panel(
    *,
    name: str,
    manifest: pd.DataFrame,
    rtab: pd.DataFrame,
    distance: pd.DataFrame,
    feature_meta: pd.DataFrame,
    out: Path,
    min_af: float,
    max_af: float,
) -> dict:
    panel = out / name
    panel.mkdir(parents=True, exist_ok=True)
    manifest = manifest.drop_duplicates("assembly_ID").copy()
    manifest["assembly_ID"] = manifest["assembly_ID"].astype(str)
    ids = [x for x in manifest["assembly_ID"] if x in rtab.columns and x in distance.index]
    manifest = manifest.set_index("assembly_ID").loc[ids].reset_index()
    rtab = rtab.loc[:, ids]
    distance = distance.loc[ids, ids]
    if set(ids) != set(rtab.columns) or set(ids) != set(distance.index) or set(distance.index) != set(distance.columns):
        raise ValueError(f"{name}: exact sample-ID mismatch")
    discovery = manifest["split"].astype(str).eq("discovery").to_numpy()
    validation = manifest["split"].astype(str).eq("validation").to_numpy()
    if not discovery.any() or not validation.any():
        raise ValueError(f"{name}: empty discovery or validation")
    if set(manifest.loc[discovery, "assembly_ID"]) & set(manifest.loc[validation, "assembly_ID"]):
        raise ValueError(f"{name}: sample overlap")
    if "BioProject" in manifest:
        dproj = set(manifest.loc[discovery, "BioProject"].dropna().astype(str)) - {""}
        vproj = set(manifest.loc[validation, "BioProject"].dropna().astype(str)) - {""}
        if dproj & vproj:
            raise ValueError(f"{name}: BioProject overlap: {sorted(dproj & vproj)[:20]}")

    af = rtab.loc[:, manifest.loc[discovery, "assembly_ID"]].mean(axis=1)
    eligible = rtab.index[(af >= min_af) & (af <= max_af)]
    rtab = rtab.loc[eligible]

    # Collapse only exact full-cohort occurrence patterns. This prevents duplicate
    # sequence fragments from multiplying a single biological signal.
    pattern_to_features: dict[str, list[str]] = {}
    pattern_to_vector: dict[str, np.ndarray] = {}
    for feature, row in rtab.iterrows():
        vector = row.to_numpy(np.uint8)
        key = hashlib.sha256(vector.tobytes()).hexdigest()
        pattern_to_features.setdefault(key, []).append(str(feature))
        pattern_to_vector[key] = vector

    rows = []
    vectors = []
    meta_key = None
    for candidate in ["feature", "variant", "canonical_sequence", "Unitig_sequence"]:
        if candidate in feature_meta.columns:
            meta_key = candidate
            break
    meta_lookup = (
        feature_meta.drop_duplicates(meta_key).set_index(meta_key)
        if meta_key is not None
        else pd.DataFrame()
    )
    for i, key in enumerate(sorted(pattern_to_features), start=1):
        members = sorted(pattern_to_features[key])
        representative = members[0]
        pattern_id = f"{name.upper()}_PATTERN_{i:05d}"
        vector = pattern_to_vector[key]
        rec = {
            "pattern_id": pattern_id,
            "representative_feature": representative,
            "n_member_features": len(members),
            "member_features": ";".join(members),
            "occurrence_sha256": key,
            "all_carriers": int(vector.sum()),
            "discovery_carriers": int(vector[discovery].sum()),
            "validation_carriers": int(vector[validation].sum()),
            "discovery_af": float(vector[discovery].mean()),
            "validation_af": float(vector[validation].mean()),
        }
        if not meta_lookup.empty and representative in meta_lookup.index:
            for col, value in meta_lookup.loc[representative].items():
                if col not in rec:
                    rec[f"meta_{col}"] = value
        rows.append(rec)
        vectors.append(vector)

    pattern_meta = pd.DataFrame(rows)
    if not vectors:
        raise ValueError(f"{name}: no eligible feature patterns")
    pattern_matrix = pd.DataFrame(
        np.column_stack(vectors),
        index=ids,
        columns=pattern_meta["pattern_id"].astype(str),
        dtype=np.uint8,
    )
    manifest.to_csv(panel / "manifest.csv", index=False)
    pattern_meta.to_csv(panel / "pattern_metadata.csv", index=False)
    distance.to_csv(panel / "distance.tsv", sep="\t")
    for split, mask in [("discovery", discovery), ("validation", validation), ("all", np.ones(len(ids), dtype=bool))]:
        sub_ids = manifest.loc[mask, "assembly_ID"].astype(str).tolist()
        x = pattern_matrix.loc[sub_ids]
        x.insert(0, "sample_id", x.index)
        x.to_csv(panel / f"{split}_genotypes.tsv", sep="\t", index=False)
        ph = manifest.loc[mask, ["assembly_ID", "phenotype"]].copy()
        ph.columns = ["sample_id", "phenotype"]
        ph["phenotype_binary"] = ph["phenotype"].astype(str).eq("R").astype(int)
        ph.to_csv(panel / f"{split}_phenotypes.tsv", sep="\t", index=False)
        distance.loc[sub_ids, sub_ids].to_csv(panel / f"{split}_distance.tsv", sep="\t")

    summary = {
        "panel": name,
        "n_samples": int(len(manifest)),
        "n_discovery": int(discovery.sum()),
        "n_validation": int(validation.sum()),
        "discovery_counts": manifest.loc[discovery, "phenotype"].value_counts().to_dict(),
        "validation_counts": manifest.loc[validation, "phenotype"].value_counts().to_dict(),
        "n_input_features": int(len(af)),
        "n_eligible_features": int(len(eligible)),
        "n_unique_occurrence_patterns": int(len(pattern_meta)),
        "n_collapsed_duplicates": int(len(eligible) - len(pattern_meta)),
        "min_discovery_af": min_af,
        "max_discovery_af": max_af,
        "bioproject_disjoint": True,
        "boundary": "Features were filtered and collapsed without inspecting phenotype associations. Pattern significance and replication remain separate tests.",
    }
    (panel / "PANEL_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # GitHub Actions artifacts are extracted at their artifact root, so the
    # top-level workflow directory name is not retained. Resolve unique filenames
    # inside each separately downloaded artifact rather than assuming a prefix.
    extended_manifest = pd.read_csv(find_one(args.extended_root, "gwas_sample_manifest.csv"), dtype=str)
    extended_rtab = load_rtab(find_one(args.extended_root, "all_variants.Rtab"))
    extended_distance_path = find_one(args.extended_root, "all_mash.tsv")
    extended_distance = pd.read_csv(extended_distance_path, sep="\t", index_col=0)
    extended_distance.index = extended_distance.index.astype(str)
    extended_distance.columns = extended_distance.columns.astype(str)
    extended_meta = pd.read_csv(find_one(args.extended_root, "targeted_variant_metadata.csv"), dtype=str)

    unitig_manifest = pd.read_csv(find_one(args.unitig_root, "gwas_sample_manifest.csv"), dtype=str)
    unitig_rtab = load_rtab(find_one(args.unitig_root, "all.rtab"))
    unitig_distance = pd.read_csv(find_one(args.unitig_root, "all_mash.tsv"), sep="\t", index_col=0)
    unitig_distance.index = unitig_distance.index.astype(str)
    unitig_distance.columns = unitig_distance.columns.astype(str)
    unitig_selection = pd.read_csv(find_one(args.unitig_root, "SELECTED_DISCOVERY_UNITIGS.csv"), dtype=str)
    unitig_meta = unitig_selection.rename(columns={"canonical_sequence": "feature"}).copy()

    summaries = [
        write_panel(
            name="targeted",
            manifest=extended_manifest,
            rtab=extended_rtab,
            distance=extended_distance,
            feature_meta=extended_meta,
            out=out,
            min_af=args.min_discovery_af,
            max_af=args.max_discovery_af,
        ),
        write_panel(
            name="unitig",
            manifest=unitig_manifest,
            rtab=unitig_rtab,
            distance=unitig_distance,
            feature_meta=unitig_meta,
            out=out,
            min_af=args.min_discovery_af,
            max_af=args.max_discovery_af,
        ),
    ]
    (out / "TREEWAS_PANEL_PREPARATION_SUMMARY.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False) + "\n"
    )
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
