#!/usr/bin/env python3
"""Focused audit for portal-unexplained colistin resistance in K. pneumoniae.

This stage does not claim biological novelty. It reconstructs a strict phenotype set,
extracts every AMR Portal genotype row for the same assemblies, reports exact
same-antibiotic joins, searches all genotype fields for known polymyxin/colistin
mechanism terms, and creates a balanced assembly manifest for sequence-level audit.
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
TARGET_ORGANISM = "Klebsiella pneumoniae"
TARGET_ANTIBIOTIC = "colistin"
TARGET_ONTOLOGY = "ARO_0000067"
KNOWN_TERMS = [
    "colistin", "polymyxin", "mcr-", "mcr_", "mgrb", "phoP", "phoQ",
    "pmrA", "pmrB", "pmrC", "pmrD", "crrA", "crrB", "arnA", "arnB",
    "arnC", "arnD", "arnT", "pbg", "lpxA", "lpxC", "lpxD", "lpxM",
    "eptA", "ugd", "pagP", "ramA", "acrAB", "oqxAB",
]


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="kp_colistin_portal_audit")
    p.add_argument("--release", default=None)
    p.add_argument("--per-class", type=int, default=358)
    return p.parse_args()


def latest_release() -> str:
    txt = requests.get(f"{ROOT}/", timeout=90).text
    vals = sorted(set(re.findall(r'href=["\'](20\d{2}-\d{2})/?["\']', txt)))
    if not vals:
        raise RuntimeError("Could not discover AMR Portal release")
    return vals[-1]


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str) + "\n")


def main() -> None:
    a = args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    release = a.release or latest_release()
    base = f"{ROOT}/{release}"
    phen = f"{base}/phenotype.parquet"
    geno = f"{base}/genotype.parquet"
    merged = f"{base}/phenotype_genotype_merged.parquet"

    con = duckdb.connect(str(out / "target_audit.duckdb"))
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("SET enable_progress_bar=false")
    con.execute("SET preserve_insertion_order=false")
    con.execute("SET threads=4")

    for label, url in [("phenotype", phen), ("genotype", geno), ("merged", merged)]:
        schema = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{url}')").fetchdf()
        schema.to_csv(out / f"{label}_schema.tsv", sep="\t", index=False)

    con.execute(
        f"""
        CREATE OR REPLACE TABLE target_labels AS
        WITH normalized AS (
          SELECT
            BioSample_ID, assembly_ID, organism, species, genus,
            antibiotic_name, antibiotic_ontology,
            collection_year, ISO_country_code, country,
            geographical_region, geographical_subregion,
            host, isolation_source, isolation_source_category,
            AMR_associated_publications,
            CASE
              WHEN lower(trim(resistance_phenotype))='resistant' THEN 'R'
              WHEN lower(trim(resistance_phenotype))='susceptible' THEN 'S'
              ELSE NULL
            END AS y
          FROM read_parquet('{phen}')
          WHERE organism=? AND lower(trim(antibiotic_name))=?
            AND antibiotic_ontology=? AND assembly_ID IS NOT NULL
        ), agg AS (
          SELECT
            BioSample_ID, assembly_ID, organism, species, genus,
            antibiotic_name, antibiotic_ontology,
            count(DISTINCT y) FILTER (WHERE y IS NOT NULL) AS n_labels,
            max(y) FILTER (WHERE y IS NOT NULL) AS phenotype,
            any_value(collection_year) AS collection_year,
            any_value(ISO_country_code) AS ISO_country_code,
            any_value(country) AS country,
            any_value(geographical_region) AS geographical_region,
            any_value(geographical_subregion) AS geographical_subregion,
            any_value(host) AS host,
            any_value(isolation_source) AS isolation_source,
            any_value(isolation_source_category) AS isolation_source_category,
            any_value(AMR_associated_publications) AS AMR_associated_publications,
            count(*) AS phenotype_source_rows
          FROM normalized GROUP BY ALL
        ), exact AS (
          SELECT DISTINCT BioSample_ID, assembly_ID, antibiotic_ontology
          FROM read_parquet('{merged}')
          WHERE antibiotic_ontology=?
        )
        SELECT a.* EXCLUDE(n_labels),
               (e.assembly_ID IS NOT NULL) AS has_exact_same_antibiotic_determinant,
               NOT (e.assembly_ID IS NOT NULL) AS portal_unexplained
        FROM agg a
        LEFT JOIN exact e USING (BioSample_ID, assembly_ID, antibiotic_ontology)
        WHERE n_labels=1 AND phenotype IN ('R','S')
        """,
        [TARGET_ORGANISM, TARGET_ANTIBIOTIC, TARGET_ONTOLOGY, TARGET_ONTOLOGY],
    )
    labels = con.execute("SELECT * FROM target_labels ORDER BY phenotype, assembly_ID").fetchdf()
    labels.to_csv(out / "target_labels.csv", index=False)

    con.execute(
        """
        CREATE OR REPLACE TABLE balanced_manifest AS
        WITH x AS (
          SELECT *,
            row_number() OVER (
              PARTITION BY phenotype, coalesce(ISO_country_code,'UNK'), coalesce(collection_year,-1)
              ORDER BY hash(assembly_ID, BioSample_ID)
            ) AS within_stratum
          FROM target_labels
          WHERE portal_unexplained AND phenotype IN ('R','S')
        ), y AS (
          SELECT *, row_number() OVER (
            PARTITION BY phenotype ORDER BY within_stratum, hash(assembly_ID, BioSample_ID)
          ) AS class_rank
          FROM x
        )
        SELECT * EXCLUDE(within_stratum, class_rank)
        FROM y WHERE class_rank <= ?
        ORDER BY phenotype, assembly_ID
        """,
        [a.per_class],
    )
    manifest = con.execute("SELECT * FROM balanced_manifest").fetchdf()
    manifest.to_csv(out / "balanced_manifest.csv", index=False)
    (out / "assembly_accessions.txt").write_text("\n".join(manifest.assembly_ID.astype(str)) + "\n")

    # Extract all genotype records for every target assembly, not merely the balanced subset.
    con.execute(
        f"""
        COPY (
          SELECT g.*
          FROM read_parquet('{geno}') g
          SEMI JOIN (SELECT DISTINCT assembly_ID FROM target_labels) t USING (assembly_ID)
        ) TO '{(out / 'target_all_genotypes.parquet').as_posix()}'
        (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )
    genotypes = con.execute(
        f"SELECT * FROM read_parquet('{(out / 'target_all_genotypes.parquet').as_posix()}')"
    ).fetchdf()
    genotypes.head(200).to_csv(out / "target_genotype_preview.tsv", sep="\t", index=False)

    object_cols = [c for c in genotypes.columns if genotypes[c].dtype == "object"]
    text = genotypes[object_cols].fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
    mask = pd.Series(False, index=genotypes.index)
    term_counts: dict[str, int] = {}
    for term in KNOWN_TERMS:
        hit = text.str.contains(re.escape(term.lower()), regex=True, na=False)
        term_counts[term] = int(hit.sum())
        mask |= hit
    known_rows = genotypes.loc[mask].copy()
    known_rows.to_csv(out / "known_colistin_keyword_genotypes.tsv", sep="\t", index=False)

    # Assembly-level summary of all portal genotype rows and known-keyword rows.
    all_counts = genotypes.groupby("assembly_ID").size().rename("n_all_portal_genotype_rows")
    known_counts = known_rows.groupby("assembly_ID").size().rename("n_known_colistin_keyword_rows")
    audit = labels.set_index("assembly_ID").join(all_counts).join(known_counts).fillna(
        {"n_all_portal_genotype_rows": 0, "n_known_colistin_keyword_rows": 0}
    ).reset_index()
    audit.to_csv(out / "assembly_level_portal_audit.csv", index=False)

    merged_col = con.execute(
        f"""
        SELECT * FROM read_parquet('{merged}')
        WHERE organism=? AND lower(trim(antibiotic_name))=? AND antibiotic_ontology=?
        ORDER BY resistance_phenotype, assembly_ID
        """,
        [TARGET_ORGANISM, TARGET_ANTIBIOTIC, TARGET_ONTOLOGY],
    ).fetchdf()
    merged_col.to_csv(out / "exact_colistin_merged_rows.tsv", sep="\t", index=False)

    summary = {
        "release": release,
        "target": {
            "organism": TARGET_ORGANISM,
            "antibiotic": TARGET_ANTIBIOTIC,
            "antibiotic_ontology": TARGET_ONTOLOGY,
        },
        "strict_label_counts": labels.groupby(
            ["phenotype", "has_exact_same_antibiotic_determinant", "portal_unexplained"],
            dropna=False,
        ).size().reset_index(name="n").to_dict(orient="records"),
        "balanced_manifest_counts": manifest.groupby("phenotype").size().to_dict(),
        "target_assemblies": int(labels.assembly_ID.nunique()),
        "portal_genotype_rows_for_target_assemblies": int(len(genotypes)),
        "assemblies_with_any_portal_genotype_row": int(genotypes.assembly_ID.nunique()),
        "known_keyword_genotype_rows": int(len(known_rows)),
        "assemblies_with_known_keyword_row": int(known_rows.assembly_ID.nunique()),
        "keyword_row_counts": term_counts,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "boundary": (
            "Keyword and exact-ontology audits only. Sequence-level reannotation of mgrB, "
            "PhoPQ/PmrAB/CrrAB and other chromosomal mechanisms remains required."
        ),
    }
    write_json(out / "TARGET_PORTAL_AUDIT_SUMMARY.json", summary)

    # Compact column-wise value audit to understand the genotype schema without assumptions.
    column_audit = []
    for c in genotypes.columns:
        s = genotypes[c]
        column_audit.append({
            "column": c,
            "dtype": str(s.dtype),
            "nonnull": int(s.notna().sum()),
            "nunique": int(s.nunique(dropna=True)),
            "examples": [str(x)[:300] for x in s.dropna().drop_duplicates().head(10).tolist()],
        })
    write_json(out / "GENOTYPE_COLUMN_AUDIT.json", column_audit)

    report = [
        "# Focused K. pneumoniae–colistin portal audit",
        "",
        f"- AMR Portal release: **{release}**",
        f"- Strict labeled assembly–colistin pairs: **{len(labels):,}**",
        f"- Portal-unexplained resistant: **{int(((labels.phenotype=='R') & labels.portal_unexplained).sum()):,}**",
        f"- Portal-unexplained susceptible: **{int(((labels.phenotype=='S') & labels.portal_unexplained).sum()):,}**",
        f"- All genotype rows for target assemblies: **{len(genotypes):,}**",
        f"- Assemblies with known-colistin keyword rows: **{known_rows.assembly_ID.nunique():,}**",
        "",
        "## Interpretation boundary",
        "",
        "This stage only checks the portal's exact antibiotic-ontology join and text-searches all genotype fields. It cannot detect disruptive mgrB insertions/deletions, promoter changes, unlisted chromosomal substitutions, expression effects, or low-quality assembly artefacts. The balanced manifest therefore advances to raw sequence reannotation and population-structure-aware analysis; no marker is called novel here.",
        "",
        "## Keyword row counts",
        "",
        pd.DataFrame(sorted(term_counts.items(), key=lambda kv: (-kv[1], kv[0])), columns=["term", "rows"]).to_markdown(index=False),
        "",
    ]
    (out / "TARGET_PORTAL_AUDIT_REPORT.md").write_text("\n".join(report))

    hashes = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    con.close()
    print((out / "TARGET_PORTAL_AUDIT_REPORT.md").read_text())


if __name__ == "__main__":
    main()
