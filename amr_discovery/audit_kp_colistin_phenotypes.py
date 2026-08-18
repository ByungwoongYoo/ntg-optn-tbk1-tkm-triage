#!/usr/bin/env python3
"""Audit K. pneumoniae colistin phenotype quality in the AMR Portal.

This script does not reinterpret susceptibility breakpoints. It preserves all source
records, classifies the deposited test-method strings, checks whether raw and portal-
updated S/R labels disagree, identifies assembly-level conflicts, and emits a strict
broth-microdilution-supported cohort for sensitivity analyses.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd
import requests

ROOT = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases"
ORGANISM = "Klebsiella pneumoniae"
ANTIBIOTIC = "colistin"
ONTOLOGY = "ARO_0000067"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="kp_colistin_phenotype_qc")
    p.add_argument("--release", default=None)
    return p.parse_args()


def latest_release() -> str:
    r = requests.get(f"{ROOT}/", timeout=90)
    r.raise_for_status()
    vals = sorted(set(re.findall(r'href=["\'](20\d{2}-\d{2})/?["\']', r.text)))
    if not vals:
        raise RuntimeError("No AMR Portal release directory found")
    return vals[-1]


def classify_method(value: Any) -> str:
    s = "" if pd.isna(value) else str(value).strip().lower()
    if not s:
        return "unknown"
    if re.search(r"broth.{0,20}micro.?dilution|micro.?broth|micro.?dilution|sensititre|micronaut", s):
        return "reference_or_commercial_bmd"
    if re.search(r"agar.{0,20}dilution", s):
        return "agar_dilution"
    if re.search(r"vitek|phoenix|microscan|automated", s):
        return "automated_system"
    if re.search(r"e.?test|etest|epsilometer|gradient", s):
        return "gradient_diffusion"
    if re.search(r"disc|disk|kirby", s):
        return "disk_diffusion"
    if re.search(r"pcr|sequenc|genotyp", s):
        return "non_ast_or_genotypic"
    return "other_or_unresolved"


def norm_sir(value: Any) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip().lower()
    if s in {"resistant", "r", "non-susceptible", "nonsusceptible", "non susceptible"}:
        return "R"
    if s in {"susceptible", "s", "sensitive"}:
        return "S"
    if s in {"intermediate", "i"}:
        return "I"
    return None


def numeric_measurement(value: Any) -> float | None:
    if pd.isna(value):
        return None
    s = str(value).strip().replace(",", ".")
    m = re.search(r"[-+]?\d+(?:\.\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def first_nonempty(values: pd.Series) -> str | None:
    for x in values:
        if pd.notna(x) and str(x).strip():
            return str(x)
    return None


def join_unique(values: pd.Series) -> str:
    vals = sorted({str(x).strip() for x in values if pd.notna(x) and str(x).strip()})
    return ";".join(vals)


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str) + "\n")


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    release = args.release or latest_release()
    url = f"{ROOT}/{release}/phenotype.parquet"

    con = duckdb.connect()
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("SET enable_progress_bar=false")
    rows = con.execute(
        f"""
        SELECT *
        FROM read_parquet('{url}')
        WHERE organism=?
          AND lower(trim(antibiotic_name))=?
          AND antibiotic_ontology=?
          AND assembly_ID IS NOT NULL
        ORDER BY assembly_ID, BioSample_ID
        """,
        [ORGANISM, ANTIBIOTIC, ONTOLOGY],
    ).fetchdf()
    con.close()
    if rows.empty:
        raise RuntimeError("No target phenotype rows found")

    rows["raw_sir"] = rows["resistance_phenotype"].map(norm_sir)
    rows["clsi_sir"] = rows["Updated_phenotype_CLSI"].map(norm_sir)
    rows["eucast_sir"] = rows["Updated_phenotype_EUCAST"].map(norm_sir)
    rows["method_class"] = rows["laboratory_typing_method"].map(classify_method)
    rows["measurement_numeric"] = rows["measurement"].map(numeric_measurement)
    rows["has_numeric_mic_like_measurement"] = (
        rows["measurement_numeric"].notna()
        & rows["measurement_units"].fillna("").astype(str).str.lower().str.contains(r"mg|µg|ug", regex=True)
    )
    rows["raw_vs_clsi_conflict"] = (
        rows.raw_sir.isin(["R", "S"])
        & rows.clsi_sir.isin(["R", "S"])
        & (rows.raw_sir != rows.clsi_sir)
    )
    rows["raw_vs_eucast_conflict"] = (
        rows.raw_sir.isin(["R", "S"])
        & rows.eucast_sir.isin(["R", "S"])
        & (rows.raw_sir != rows.eucast_sir)
    )
    rows["portal_updated_conflict"] = rows.raw_vs_clsi_conflict | rows.raw_vs_eucast_conflict
    rows["bmd_supported_row"] = (
        rows.method_class.eq("reference_or_commercial_bmd")
        & rows.raw_sir.isin(["R", "S"])
    )
    rows.to_csv(out / "all_target_phenotype_rows.csv", index=False)

    method_summary = (
        rows.groupby(
            ["method_class", "laboratory_typing_method", "raw_sir"], dropna=False
        )
        .agg(
            rows=("assembly_ID", "size"),
            assemblies=("assembly_ID", "nunique"),
            biosamples=("BioSample_ID", "nunique"),
            numeric_measurements=("measurement_numeric", lambda x: int(x.notna().sum())),
            portal_updated_conflicts=("portal_updated_conflict", "sum"),
        )
        .reset_index()
        .sort_values(["method_class", "assemblies"], ascending=[True, False])
    )
    method_summary.to_csv(out / "method_summary.csv", index=False)

    publication_summary = (
        rows.assign(publication=rows.AMR_associated_publications.fillna("UNSPECIFIED"))
        .groupby(["publication", "method_class", "raw_sir"], dropna=False)
        .agg(
            rows=("assembly_ID", "size"),
            assemblies=("assembly_ID", "nunique"),
            countries=("ISO_country_code", "nunique"),
            years=("collection_year", "nunique"),
            portal_updated_conflicts=("portal_updated_conflict", "sum"),
        )
        .reset_index()
        .sort_values(["assemblies", "publication"], ascending=[False, True])
    )
    publication_summary.to_csv(out / "publication_method_summary.csv", index=False)

    records: list[dict[str, Any]] = []
    for (assembly, biosample), g in rows.groupby(["assembly_ID", "BioSample_ID"], dropna=False):
        raw_set = sorted(set(g.raw_sir.dropna()) & {"R", "S"})
        clsi_set = sorted(set(g.clsi_sir.dropna()) & {"R", "S"})
        eucast_set = sorted(set(g.eucast_sir.dropna()) & {"R", "S"})
        strict_raw = raw_set[0] if len(raw_set) == 1 else None
        any_updated_conflict = bool(g.portal_updated_conflict.any())
        bmd = g[g.bmd_supported_row]
        bmd_set = sorted(set(bmd.raw_sir.dropna()) & {"R", "S"})
        bmd_label = bmd_set[0] if len(bmd_set) == 1 else None
        records.append(
            {
                "assembly_ID": assembly,
                "BioSample_ID": biosample,
                "strict_raw_label": strict_raw,
                "raw_label_set": ";".join(raw_set),
                "clsi_label_set": ";".join(clsi_set),
                "eucast_label_set": ";".join(eucast_set),
                "source_row_count": int(len(g)),
                "assembly_label_conflict": len(raw_set) > 1,
                "portal_updated_conflict": any_updated_conflict,
                "bmd_supported_label": bmd_label,
                "bmd_supported_rows": int(len(bmd)),
                "bmd_label_conflict": len(bmd_set) > 1,
                "has_numeric_measurement": bool(g.measurement_numeric.notna().any()),
                "methods": join_unique(g.laboratory_typing_method),
                "method_classes": join_unique(g.method_class),
                "ast_standards": join_unique(g.ast_standard),
                "publications": join_unique(g.AMR_associated_publications),
                "collection_year": first_nonempty(g.collection_year),
                "ISO_country_code": first_nonempty(g.ISO_country_code),
                "country": first_nonempty(g.country),
                "isolation_source": first_nonempty(g.isolation_source),
            }
        )
    assembly = pd.DataFrame(records)
    assembly["strict_portal_label"] = assembly.strict_raw_label.where(
        assembly.strict_raw_label.isin(["R", "S"])
        & ~assembly.assembly_label_conflict
        & ~assembly.portal_updated_conflict
    )
    assembly["strict_bmd_label"] = assembly.bmd_supported_label.where(
        assembly.bmd_supported_label.isin(["R", "S"])
        & ~assembly.bmd_label_conflict
        & ~assembly.assembly_label_conflict
        & ~assembly.portal_updated_conflict
    )
    assembly.to_csv(out / "assembly_level_phenotype_qc.csv", index=False)
    assembly[assembly.strict_portal_label.isin(["R", "S"])].to_csv(
        out / "strict_portal_labels.csv", index=False
    )
    assembly[assembly.strict_bmd_label.isin(["R", "S"])].to_csv(
        out / "strict_bmd_supported_labels.csv", index=False
    )

    counts = {
        "release": release,
        "target": {"organism": ORGANISM, "antibiotic": ANTIBIOTIC, "ontology": ONTOLOGY},
        "source_rows": int(len(rows)),
        "assemblies": int(rows.assembly_ID.nunique()),
        "biosamples": int(rows.BioSample_ID.nunique()),
        "raw_R_rows": int((rows.raw_sir == "R").sum()),
        "raw_S_rows": int((rows.raw_sir == "S").sum()),
        "rows_with_reference_or_commercial_bmd_method": int(rows.bmd_supported_row.sum()),
        "assemblies_with_bmd_supported_row": int(rows.loc[rows.bmd_supported_row, "assembly_ID"].nunique()),
        "assembly_level_raw_label_conflicts": int(assembly.assembly_label_conflict.sum()),
        "assembly_level_portal_updated_conflicts": int(assembly.portal_updated_conflict.sum()),
        "strict_portal_R": int((assembly.strict_portal_label == "R").sum()),
        "strict_portal_S": int((assembly.strict_portal_label == "S").sum()),
        "strict_bmd_R": int((assembly.strict_bmd_label == "R").sum()),
        "strict_bmd_S": int((assembly.strict_bmd_label == "S").sum()),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "boundary": (
            "Method classes are based on deposited free-text labels. No MIC breakpoint was reinterpreted. "
            "BMD-supported means that at least one deposited method string matched a broth-microdilution or commercial-BMD term; it does not independently verify laboratory execution."
        ),
    }
    write_json(out / "PHENOTYPE_QC_SUMMARY.json", counts)

    report = [
        "# K. pneumoniae–colistin phenotype quality audit",
        "",
        f"- AMR Portal release: **{release}**",
        f"- Source rows / distinct assemblies: **{counts['source_rows']:,} / {counts['assemblies']:,}**",
        f"- Strict portal labels R / S: **{counts['strict_portal_R']:,} / {counts['strict_portal_S']:,}**",
        f"- Strict BMD-supported labels R / S: **{counts['strict_bmd_R']:,} / {counts['strict_bmd_S']:,}**",
        f"- Assembly-level raw-label conflicts: **{counts['assembly_level_raw_label_conflicts']:,}**",
        f"- Assembly-level raw-versus-updated conflicts: **{counts['assembly_level_portal_updated_conflicts']:,}**",
        "",
        "## Claim boundary",
        "",
        counts["boundary"],
        "",
        "## Deposited method summary",
        "",
        method_summary.head(50).to_markdown(index=False),
        "",
    ]
    (out / "PHENOTYPE_QC_REPORT.md").write_text("\n".join(report))

    hashes = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print((out / "PHENOTYPE_QC_REPORT.md").read_text())


if __name__ == "__main__":
    main()
