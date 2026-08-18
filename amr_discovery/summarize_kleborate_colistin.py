#!/usr/bin/env python3
"""Merge Kleborate v3 results with strict colistin phenotype labels.

The script treats acquired mcr-family genes and Kleborate-reported MgrB/PmrB
inactivating mutations as *known-mechanism evidence*. It does not assume that their
absence excludes PhoPQ/PmrAB/CrrAB substitutions or other known mechanisms; the
remaining resistant assemblies are only a residual sequence-audit cohort.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--labels", required=True)
    p.add_argument("--kleborate-dir", required=True)
    p.add_argument("--out", default="kp_colistin_kleborate_summary")
    return p.parse_args()


def is_nonempty_marker(value: Any) -> bool:
    if pd.isna(value):
        return False
    x = str(value).strip()
    return x not in {"", "-", ".", "None", "nan", "0"}


def extract_accession(row: pd.Series) -> str | None:
    # Kleborate normally reports the input filename, but this is intentionally
    # schema-agnostic across v3 point releases.
    preferred = [c for c in row.index if c.lower() in {"strain", "assembly", "input", "file", "filename"}]
    for c in preferred + list(row.index):
        m = ACC_RE.search(str(row[c]))
        if m:
            return m.group(0)
    return None


def find_column(columns: list[str], exact: list[str], contains: list[str]) -> str | None:
    low = {c.lower(): c for c in columns}
    for e in exact:
        if e.lower() in low:
            return low[e.lower()]
    for c in columns:
        lc = c.lower()
        if all(x.lower() in lc for x in contains):
            return c
    return None


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str) + "\n")


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(a.labels)

    files = sorted(Path(a.kleborate_dir).glob("*.txt"))
    if not files:
        raise SystemExit(f"No Kleborate TSV output found under {a.kleborate_dir}")
    frames = []
    for p in files:
        try:
            x = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            continue
        if len(x):
            x["_kleborate_output_file"] = p.name
            frames.append(x)
    if not frames:
        raise SystemExit("Kleborate output files contained no rows")
    k = pd.concat(frames, ignore_index=True, sort=False)
    k["assembly_ID"] = k.apply(extract_accession, axis=1)
    missing = k.assembly_ID.isna().sum()
    if missing:
        k[k.assembly_ID.isna()].to_csv(out / "unresolved_kleborate_rows.tsv", sep="\t", index=False)
    k = k.dropna(subset=["assembly_ID"]).copy()
    k.to_csv(out / "kleborate_all_outputs.tsv", sep="\t", index=False)

    cols = list(k.columns)
    st_col = find_column(cols, ["ST"], ["st"])
    species_col = find_column(cols, ["species"], ["species"])
    colmut_col = find_column(cols, ["Col_mutations"], ["col", "mutation"])

    object_text = k.fillna("").astype(str).agg(" | ".join, axis=1)
    k["has_mcr_family_call"] = object_text.str.contains(
        r"(?i)(^|[^A-Za-z0-9])mcr[-_]?\d", regex=True, na=False
    )
    if colmut_col:
        k["has_kleborate_col_mutation"] = k[colmut_col].map(is_nonempty_marker)
    else:
        k["has_kleborate_col_mutation"] = False
    k["has_kleborate_known_colistin_evidence"] = (
        k.has_mcr_family_call | k.has_kleborate_col_mutation
    )

    keep = ["assembly_ID", "has_mcr_family_call", "has_kleborate_col_mutation", "has_kleborate_known_colistin_evidence"]
    if st_col:
        k = k.rename(columns={st_col: "Kleborate_ST"})
        keep.append("Kleborate_ST")
    if species_col:
        k = k.rename(columns={species_col: "Kleborate_species"})
        keep.append("Kleborate_species")
    if colmut_col and colmut_col != "Col_mutations":
        k = k.rename(columns={colmut_col: "Col_mutations"})
    if "Col_mutations" in k.columns:
        keep.append("Col_mutations")

    # Multiple species-output rows for an accession are unexpected; collapse
    # deterministically while preserving all marker strings.
    def join_unique(s: pd.Series) -> str:
        vals = sorted({str(x) for x in s if is_nonempty_marker(x)})
        return ";".join(vals) if vals else "-"

    agg: dict[str, Any] = {
        "has_mcr_family_call": "max",
        "has_kleborate_col_mutation": "max",
        "has_kleborate_known_colistin_evidence": "max",
    }
    for c in ["Kleborate_ST", "Kleborate_species", "Col_mutations"]:
        if c in k.columns:
            agg[c] = join_unique
    kk = k[keep].groupby("assembly_ID", as_index=False).agg(agg)

    merged = labels.merge(kk, on="assembly_ID", how="left", validate="one_to_one")
    for c in ["has_mcr_family_call", "has_kleborate_col_mutation", "has_kleborate_known_colistin_evidence"]:
        merged[c] = merged[c].fillna(False).astype(bool)
    merged["kleborate_result_present"] = merged.assembly_ID.isin(set(kk.assembly_ID))
    merged.to_csv(out / "labels_with_kleborate.csv", index=False)

    portal_r0 = (merged.phenotype == "R") & merged.portal_unexplained.astype(bool)
    portal_s0 = (merged.phenotype == "S") & merged.portal_unexplained.astype(bool)
    residual_r = portal_r0 & ~merged.has_kleborate_known_colistin_evidence
    residual_s = portal_s0 & ~merged.has_kleborate_known_colistin_evidence
    merged.loc[residual_r].to_csv(out / "residual_resistant_after_kleborate.csv", index=False)
    merged.loc[residual_s].to_csv(out / "residual_susceptible_after_kleborate.csv", index=False)

    # Match susceptible controls by ST first, then country/year. This is a cohort
    # construction step only; the downstream association model still uses genomic PCs.
    r = merged.loc[residual_r].copy()
    s = merged.loc[residual_s].copy()
    st = "Kleborate_ST" if "Kleborate_ST" in merged.columns else None
    chosen = []
    used: set[str] = set()
    for _, rr in r.sort_values([c for c in [st, "ISO_country_code", "collection_year", "assembly_ID"] if c]).iterrows():
        pool = s[~s.assembly_ID.isin(used)]
        if st and is_nonempty_marker(rr.get(st)):
            same = pool[pool[st] == rr[st]]
            if len(same):
                pool = same
        same_country = pool[pool.ISO_country_code.fillna("UNK") == (rr.get("ISO_country_code") if pd.notna(rr.get("ISO_country_code")) else "UNK")]
        if len(same_country):
            pool = same_country
        if len(pool) == 0:
            continue
        yr = rr.get("collection_year")
        temp = pool.copy()
        if pd.notna(yr):
            temp["_year_distance"] = (pd.to_numeric(temp.collection_year, errors="coerce") - float(yr)).abs().fillna(9999)
        else:
            temp["_year_distance"] = 9999
        pick = temp.sort_values(["_year_distance", "assembly_ID"]).iloc[0]
        used.add(pick.assembly_ID)
        chosen.append(pick.drop(labels=["_year_distance"]).to_dict())
    matched_s = pd.DataFrame(chosen)
    residual_cohort = pd.concat([r, matched_s], ignore_index=True, sort=False)
    residual_cohort.to_csv(out / "residual_matched_sequence_audit_cohort.csv", index=False)
    (out / "residual_sequence_audit_accessions.txt").write_text(
        "\n".join(residual_cohort.assembly_ID.astype(str)) + "\n"
    )

    summary = {
        "kleborate_output_files": [p.name for p in files],
        "kleborate_rows": int(len(k)),
        "resolved_assemblies": int(kk.assembly_ID.nunique()),
        "unresolved_rows": int(missing),
        "portal_unexplained_resistant": int(portal_r0.sum()),
        "portal_unexplained_susceptible": int(portal_s0.sum()),
        "r0_with_mcr_family_call": int((portal_r0 & merged.has_mcr_family_call).sum()),
        "r0_with_mgrb_or_pmrb_inactivation": int((portal_r0 & merged.has_kleborate_col_mutation).sum()),
        "r0_with_any_kleborate_known_colistin_evidence": int((portal_r0 & merged.has_kleborate_known_colistin_evidence).sum()),
        "residual_r_after_kleborate": int(residual_r.sum()),
        "residual_s_after_kleborate": int(residual_s.sum()),
        "matched_residual_s": int(len(matched_s)),
        "residual_sequence_audit_cohort": int(len(residual_cohort)),
        "boundary": (
            "Residual means no mcr-family call and no Kleborate-detected MgrB/PmrB inactivation. "
            "It does not exclude known missense substitutions, promoter effects, heteroresistance, "
            "phenotype error, or other established chromosomal mechanisms."
        ),
    }
    write_json(out / "KLEBORATE_COLISTIN_SUMMARY.json", summary)

    if st and st in merged.columns:
        sttab = merged.groupby([st, "phenotype", "has_kleborate_known_colistin_evidence"], dropna=False).size().reset_index(name="n")
        sttab.to_csv(out / "ST_phenotype_mechanism_counts.csv", index=False)

    report = [
        "# K. pneumoniae colistin sequence-level known-mechanism audit",
        "",
        f"- Assemblies with Kleborate output: **{summary['resolved_assemblies']:,}**",
        f"- Portal-unexplained resistant isolates: **{summary['portal_unexplained_resistant']:,}**",
        f"- Resistant isolates with an mcr-family call: **{summary['r0_with_mcr_family_call']:,}**",
        f"- Resistant isolates with a Kleborate MgrB/PmrB inactivation call: **{summary['r0_with_mgrb_or_pmrb_inactivation']:,}**",
        f"- Resistant isolates remaining for the expanded chromosomal audit: **{summary['residual_r_after_kleborate']:,}**",
        "",
        "## Interpretation boundary",
        "",
        summary["boundary"],
        "",
    ]
    (out / "KLEBORATE_COLISTIN_REPORT.md").write_text("\n".join(report))

    hashes = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print((out / "KLEBORATE_COLISTIN_REPORT.md").read_text())


if __name__ == "__main__":
    main()
