#!/usr/bin/env python3
"""Prepare source-held-out discovery/validation inputs for targeted and unitig GWAS.

The split is deterministic and group-disjoint at the primary BioProject level whenever
metadata permit. No candidate feature is examined while choosing the split.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")
PRJ_RE = re.compile(r"PRJ[A-Z]+\d+")
BIOSAMPLE_RE = re.compile(r"SAM[NED][A-Z]?\d+")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--assembly-report", required=True)
    p.add_argument("--assemblies-dir", required=True)
    p.add_argument("--out", default="gwas_inputs")
    p.add_argument("--validation-fraction", type=float, default=0.30)
    p.add_argument("--seed", type=int, default=20260819)
    p.add_argument("--min-class", type=int, default=25)
    return p.parse_args()


def recursive_values(obj: Any, key_hint: str) -> list[str]:
    vals: list[str] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if key_hint.lower() in str(k).lower():
                if isinstance(v, (str, int, float)):
                    vals.append(str(v))
                elif isinstance(v, list):
                    for x in v:
                        if isinstance(x, (str, int, float)):
                            vals.append(str(x))
                        elif isinstance(x, dict):
                            vals.extend(str(y) for y in x.values() if isinstance(y, (str, int, float)))
                elif isinstance(v, dict):
                    vals.extend(str(x) for x in v.values() if isinstance(x, (str, int, float)))
            vals.extend(recursive_values(v, key_hint))
    elif isinstance(obj, list):
        for x in obj:
            vals.extend(recursive_values(x, key_hint))
    return vals


def _primary_bioproject(obj: dict[str, Any], all_projects: list[str]) -> str | None:
    """Prefer NCBI's explicit assembly-level primary BioProject, then BioSample projects."""
    info = obj.get("assemblyInfo") or {}
    primary = info.get("bioprojectAccession")
    if isinstance(primary, str) and PRJ_RE.fullmatch(primary):
        return primary
    biosample = info.get("biosample") or {}
    for item in biosample.get("bioprojects") or []:
        if isinstance(item, dict):
            value = item.get("accession")
            if isinstance(value, str) and PRJ_RE.fullmatch(value):
                return value
    return all_projects[0] if all_projects else None


def parse_assembly_report(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            obj = json.loads(line)
            # The top-level accession is the assembly actually returned for the requested
            # accession. Do not replace a requested GCA with the paired GCF accession.
            acc = obj.get("accession") or obj.get("currentAccession")
            if not isinstance(acc, str) or not ACC_RE.fullmatch(acc):
                text = json.dumps(obj)
                accs = ACC_RE.findall(text)
                if not accs:
                    continue
                acc = accs[0]
            text = json.dumps(obj)
            prjs = sorted(set(PRJ_RE.findall(text)))
            primary_project = _primary_bioproject(obj, prjs)
            info = obj.get("assemblyInfo") or {}
            biosample_obj = info.get("biosample") or {}
            biosample = biosample_obj.get("accession")
            if not isinstance(biosample, str) or not BIOSAMPLE_RE.fullmatch(biosample):
                bios = BIOSAMPLE_RE.findall(text)
                biosample = bios[0] if bios else None
            submitter = info.get("submitter")
            if not submitter:
                values = recursive_values(obj, "submitter") + recursive_values(obj, "organization")
                submitter = values[0] if values else None
            rows.append({
                "assembly_ID": acc,
                "BioProject": primary_project,
                "BioProjects_all": ";".join(prjs),
                "BioSample_report": biosample,
                "Submitter": str(submitter)[:300] if submitter else None,
            })
    if not rows:
        return pd.DataFrame(columns=["assembly_ID", "BioProject", "BioProjects_all", "BioSample_report", "Submitter"])
    result = pd.DataFrame(rows).drop_duplicates("assembly_ID")
    if result.assembly_ID.duplicated().any():
        raise RuntimeError("Duplicate assembly accessions after NCBI report parsing")
    return result


def find_assembly_paths(meta: pd.DataFrame, directory: Path) -> dict[str, str]:
    paths = list(directory.rglob("*.fna")); by_base: dict[str, list[Path]] = {}
    for p in paths:
        m = ACC_RE.search(str(p))
        if m: by_base.setdefault(m.group(0).split(".")[0], []).append(p)
    out = {}
    for acc in meta.assembly_ID.astype(str):
        exact = [p for p in paths if acc in str(p)]
        if exact: out[acc] = str(sorted(exact, key=lambda z: len(str(z)))[0].resolve())
        elif acc.split(".")[0] in by_base: out[acc] = str(sorted(by_base[acc.split(".")[0]], key=lambda z: len(str(z)))[0].resolve())
    return out


def choose_split(meta: pd.DataFrame, frac: float, seed: int, min_class: int) -> tuple[np.ndarray, np.ndarray, str]:
    y = (meta.phenotype.astype(str) == "R").astype(int).to_numpy()
    fallback = meta.get("AMR_associated_publications", pd.Series(index=meta.index, dtype=object)).fillna("").astype(str) + "|" + meta.get("ISO_country_code", pd.Series(index=meta.index, dtype=object)).fillna("UNK").astype(str) + "|" + meta.get("collection_year", pd.Series(index=meta.index, dtype=object)).fillna("UNK").astype(str)
    groups = meta.BioProject.fillna("").astype(str)
    n_project_rows = int(groups.str.len().gt(0).sum())
    groups = groups.where(groups.str.len() > 0, fallback)
    if n_project_rows < max(100, int(0.5 * len(meta))) or groups.nunique() < 4:
        groups, basis = fallback, "country-year-publication fallback"
    else:
        basis = "primary BioProject with fallback only for missing projects"
    splitter = GroupShuffleSplit(n_splits=5000, test_size=frac, random_state=seed)
    best = None; target_n = len(meta) * frac; target_r = y.sum() * frac; target_s = (len(y) - y.sum()) * frac
    for tr, va in splitter.split(meta, y, groups):
        yr = y[va]; r = int(yr.sum()); s = int(len(yr) - r); dr = int(y[tr].sum()); ds = int(len(tr) - dr)
        if min(r, s, dr, ds) < min_class: continue
        if set(groups.iloc[tr]).intersection(set(groups.iloc[va])): continue
        score = abs(len(va) - target_n) / max(target_n, 1) + abs(r - target_r) / max(target_r, 1) + abs(s - target_s) / max(target_s, 1) + 0.05 / max(groups.iloc[va].nunique(), 1)
        if best is None or score < best[0]: best = (score, tr, va)
    if best is None: raise RuntimeError("Could not construct a group-disjoint split with adequate R/S counts")
    return best[1], best[2], basis


def write_rtab_subset(df: pd.DataFrame, samples: list[str], path: Path) -> None:
    df.loc[:, [df.columns[0], *samples]].to_csv(path, sep="\t", index=False)


def write_pheno(meta: pd.DataFrame, path: Path) -> None:
    x = meta[["assembly_ID", "phenotype"]].copy(); x["phenotype"] = (x.phenotype.astype(str) == "R").astype(int); x.columns = ["samples", "phenotype"]; x.to_csv(path, sep="\t", index=False)


def main() -> None:
    a = parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(a.metadata, dtype={"assembly_ID": str})
    if "qc_pass" in meta.columns: meta = meta[meta.qc_pass.astype(bool)].copy()
    meta = meta.drop_duplicates("assembly_ID")
    # The NCBI assembly report is the sole authoritative source for these fields.
    # Input manifests may already contain earlier parsed copies; dropping them prevents
    # BioProject_x/BioProject_y suffixes and makes the split provenance unambiguous.
    meta = meta.drop(columns=["BioProject", "BioProjects_all", "BioSample_report", "Submitter"], errors="ignore")
    report = parse_assembly_report(Path(a.assembly_report)); report["base"] = report.assembly_ID.astype(str).str.split(".").str[0]
    meta["base"] = meta.assembly_ID.astype(str).str.split(".").str[0]; meta = meta.merge(report.drop(columns=["assembly_ID"]), on="base", how="left").drop(columns=["base"])
    rtab = pd.read_csv(a.rtab, sep="\t", dtype=str); sample_cols = set(rtab.columns[1:]); meta = meta[meta.assembly_ID.isin(sample_cols)].reset_index(drop=True)
    paths = find_assembly_paths(meta, Path(a.assemblies_dir)); meta["assembly_path"] = meta.assembly_ID.map(paths); meta = meta[meta.assembly_path.notna()].reset_index(drop=True)
    if len(meta) < 100: raise RuntimeError(f"Only {len(meta)} metadata rows have matching assembly paths")
    _, va, basis = choose_split(meta, a.validation_fraction, a.seed, a.min_class)
    meta["split"] = "discovery"; meta.loc[va, "split"] = "validation"
    fb = meta.get("AMR_associated_publications", pd.Series(index=meta.index, dtype=object)).fillna("").astype(str) + "|" + meta.get("ISO_country_code", pd.Series(index=meta.index, dtype=object)).fillna("UNK").astype(str) + "|" + meta.get("collection_year", pd.Series(index=meta.index, dtype=object)).fillna("UNK").astype(str)
    meta["source_group"] = meta.BioProject.fillna("").astype(str); meta["source_group"] = meta.source_group.where(meta.source_group.str.len() > 0, fb); meta.to_csv(out / "gwas_sample_manifest.csv", index=False)
    all_samples = meta.assembly_ID.astype(str).tolist(); discovery = meta.loc[meta.split == "discovery", "assembly_ID"].astype(str).tolist(); validation = meta.loc[meta.split == "validation", "assembly_ID"].astype(str).tolist()
    write_rtab_subset(rtab, all_samples, out / "all_variants.Rtab"); write_rtab_subset(rtab, discovery, out / "discovery_variants.Rtab"); write_rtab_subset(rtab, validation, out / "validation_variants.Rtab")
    write_pheno(meta, out / "all_phenotypes.tsv"); write_pheno(meta[meta.split == "discovery"], out / "discovery_phenotypes.tsv"); write_pheno(meta[meta.split == "validation"], out / "validation_phenotypes.tsv")
    for name, frame in [("all", meta), ("discovery", meta[meta.split == "discovery"]), ("validation", meta[meta.split == "validation"])]:
        (out / f"{name}_refs.txt").write_text("\n".join(frame.assembly_path.astype(str)) + "\n"); (out / f"{name}_samples.txt").write_text("\n".join(frame.assembly_ID.astype(str)) + "\n")
    summary = {
        "n_all": int(len(meta)), "n_discovery": int((meta.split == "discovery").sum()), "n_validation": int((meta.split == "validation").sum()),
        "class_counts_all": meta.phenotype.value_counts().to_dict(), "class_counts_discovery": meta[meta.split == "discovery"].phenotype.value_counts().to_dict(), "class_counts_validation": meta[meta.split == "validation"].phenotype.value_counts().to_dict(),
        "n_primary_bioproject_rows": int(meta.BioProject.notna().sum()), "primary_bioproject_coverage": float(meta.BioProject.notna().mean()), "n_primary_bioprojects": int(meta.BioProject.nunique(dropna=True)),
        "n_source_groups_all": int(meta.source_group.nunique()), "n_source_groups_discovery": int(meta[meta.split == "discovery"].source_group.nunique()), "n_source_groups_validation": int(meta[meta.split == "validation"].source_group.nunique()),
        "group_disjoint": not bool(set(meta[meta.split == "discovery"].source_group) & set(meta[meta.split == "validation"].source_group)), "split_basis": basis, "seed": a.seed, "validation_fraction_requested": a.validation_fraction,
        "boundary": "The split was chosen from source-group and phenotype counts only; no genomic feature was examined."
    }
    (out / "GWAS_INPUT_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    hashes = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n"); print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
