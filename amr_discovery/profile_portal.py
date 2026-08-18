#!/usr/bin/env python3
"""Profile public AMR phenotype/genotype data and select a defensible discovery target.

The script operates directly on the latest EMBL-EBI AMR Portal Parquet release with
DuckDB. It creates a strict, deduplicated assembly–antibiotic phenotype table, marks
whether the portal contains an AMRFinderPlus determinant joined to the exact same
antibiotic ontology, ranks clinically important organism–antibiotic pairs, and emits a
balanced discovery manifest for downstream genome analysis.

`portal_unexplained = True` is an operational data state, not proof that all known
biological mechanisms have been excluded.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd
import requests

RELEASE_ROOT = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases"

ANTIBIOTIC_WEIGHTS = {
    "cefiderocol": 12.0,
    "ceftazidime-avibactam": 11.0,
    "ceftazidime avibactam": 11.0,
    "colistin": 11.0,
    "polymyxin b": 10.5,
    "daptomycin": 10.0,
    "linezolid": 10.0,
    "vancomycin": 9.5,
    "meropenem": 9.5,
    "imipenem": 9.0,
    "ertapenem": 9.0,
    "ceftriaxone": 9.0,
    "cefotaxime": 8.5,
    "ceftazidime": 8.5,
    "cefepime": 8.5,
    "azithromycin": 8.5,
    "fosfomycin": 8.5,
    "tigecycline": 8.0,
    "ciprofloxacin": 8.0,
    "levofloxacin": 7.5,
    "amikacin": 7.0,
    "gentamicin": 6.5,
    "trimethoprim-sulfamethoxazole": 6.5,
    "trimethoprim sulfamethoxazole": 6.5,
}

ORGANISM_WEIGHTS = {
    "klebsiella pneumoniae": 6.0,
    "acinetobacter baumannii": 6.0,
    "pseudomonas aeruginosa": 6.0,
    "staphylococcus aureus": 6.0,
    "enterococcus faecium": 6.0,
    "neisseria gonorrhoeae": 6.0,
    "escherichia coli": 5.0,
    "salmonella enterica": 4.5,
    "campylobacter jejuni": 4.0,
    "campylobacter coli": 3.5,
    "streptococcus pneumoniae": 4.0,
    "enterococcus faecalis": 4.0,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="amr_profile_results")
    p.add_argument("--release", default=None, help="YYYY-MM; default: discover latest")
    p.add_argument("--min-resistant", type=int, default=40)
    p.add_argument("--min-susceptible", type=int, default=80)
    p.add_argument("--min-unexplained-resistant", type=int, default=20)
    p.add_argument("--manifest-per-class", type=int, default=600)
    return p.parse_args()


def discover_latest_release() -> tuple[str, str]:
    r = requests.get(f"{RELEASE_ROOT}/", timeout=90)
    r.raise_for_status()
    releases = sorted(set(re.findall(r'href=["\'](20\d{2}-\d{2})/?["\']', r.text)))
    if not releases:
        raise RuntimeError("No YYYY-MM AMR Portal release directory found")
    return releases[-1], r.text


def weight_by_substring(value: str, mapping: dict[str, float]) -> tuple[float, str]:
    low = (value or "").lower().strip()
    hits = [(w, key) for key, w in mapping.items() if key in low]
    return max(hits, default=(0.0, ""))


def clinical_score(row: pd.Series) -> tuple[float, dict[str, Any]]:
    abx_weight, abx_key = weight_by_substring(str(row.antibiotic_name), ANTIBIOTIC_WEIGHTS)
    org_weight, org_key = weight_by_substring(str(row.organism), ORGANISM_WEIGHTS)
    r0 = int(row.n_r_without_known)
    s0 = int(row.n_s_without_known)
    nr = int(row.n_r)
    exact_known_fraction = float(row.n_r_with_known) / nr if nr else 0.0
    unexplained_fraction = float(r0) / nr if nr else 0.0

    score = (
        3.2 * math.log10(r0 + 1)
        + 2.2 * math.log10(s0 + 1)
        + 0.30 * min(int(row.n_countries_r_without_known), 12)
        + 0.18 * min(int(row.n_years_r_without_known), 15)
        + abx_weight
        + org_weight
    )
    # If virtually every resistant isolate is "unexplained" and exact ontology
    # matches are absent, an ontology mismatch is more plausible than a discovery.
    ontology_penalty = 0.0
    if unexplained_fraction >= 0.95 and exact_known_fraction <= 0.02:
        ontology_penalty = 8.0
    elif unexplained_fraction >= 0.85 and exact_known_fraction <= 0.05:
        ontology_penalty = 4.0
    score -= ontology_penalty

    eligible_high_impact = abx_weight > 0 and org_weight > 0
    return score, {
        "antibiotic_priority_key": abx_key,
        "organism_priority_key": org_key,
        "antibiotic_priority_weight": abx_weight,
        "organism_priority_weight": org_weight,
        "ontology_mismatch_penalty": ontology_penalty,
        "high_impact_pair": eligible_high_impact,
        "r_exact_known_fraction": exact_known_fraction,
        "r_portal_unexplained_fraction": unexplained_fraction,
    }


def q(con: duckdb.DuckDBPyConnection, sql: str, params: list[Any] | None = None) -> pd.DataFrame:
    return con.execute(sql, params or []).fetchdf()


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.release:
        release = args.release
        listing = "release supplied explicitly"
    else:
        release, listing = discover_latest_release()
    base = f"{RELEASE_ROOT}/{release}"
    phen = f"{base}/phenotype.parquet"
    geno = f"{base}/genotype.parquet"
    merged = f"{base}/phenotype_genotype_merged.parquet"

    con = duckdb.connect(str(out / "profile.duckdb"))
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("SET enable_progress_bar=false")
    con.execute("SET preserve_insertion_order=false")
    con.execute("SET threads=4")

    # Strict phenotype consensus: one assembly/antibiotic key must have only R or only S.
    con.execute(
        f"""
        CREATE OR REPLACE TABLE phen_strict AS
        WITH normalized AS (
          SELECT
            BioSample_ID,
            assembly_ID,
            organism,
            species,
            genus,
            antibiotic_name,
            antibiotic_ontology,
            collection_year,
            ISO_country_code,
            country,
            geographical_region,
            geographical_subregion,
            host,
            isolation_source,
            isolation_source_category,
            AMR_associated_publications,
            CASE
              WHEN lower(trim(resistance_phenotype)) = 'resistant' THEN 'R'
              WHEN lower(trim(resistance_phenotype)) = 'susceptible' THEN 'S'
              ELSE NULL
            END AS y
          FROM read_parquet('{phen}')
          WHERE assembly_ID IS NOT NULL
            AND antibiotic_ontology IS NOT NULL
            AND antibiotic_name IS NOT NULL
        ), agg AS (
          SELECT
            BioSample_ID,
            assembly_ID,
            organism,
            species,
            genus,
            antibiotic_name,
            antibiotic_ontology,
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
            count(*) AS source_rows
          FROM normalized
          GROUP BY ALL
        )
        SELECT * EXCLUDE(n_labels)
        FROM agg
        WHERE n_labels = 1 AND phenotype IN ('R','S')
        """
    )

    # The portal merged view embodies its exact BioSample + assembly + antibiotic ontology join.
    con.execute(
        f"""
        CREATE OR REPLACE TABLE exact_known_keys AS
        SELECT DISTINCT BioSample_ID, assembly_ID, antibiotic_ontology
        FROM read_parquet('{merged}')
        WHERE BioSample_ID IS NOT NULL
          AND assembly_ID IS NOT NULL
          AND antibiotic_ontology IS NOT NULL
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE strict_labeled AS
        SELECT p.*,
               (k.assembly_ID IS NOT NULL) AS has_exact_known_determinant,
               NOT (k.assembly_ID IS NOT NULL) AS portal_unexplained
        FROM phen_strict p
        LEFT JOIN exact_known_keys k
          USING (BioSample_ID, assembly_ID, antibiotic_ontology)
        """
    )

    pair_stats = q(
        con,
        """
        SELECT
          organism,
          species,
          genus,
          antibiotic_name,
          antibiotic_ontology,
          count(*) FILTER (WHERE phenotype='R') AS n_r,
          count(*) FILTER (WHERE phenotype='S') AS n_s,
          count(*) FILTER (WHERE phenotype='R' AND has_exact_known_determinant) AS n_r_with_known,
          count(*) FILTER (WHERE phenotype='R' AND portal_unexplained) AS n_r_without_known,
          count(*) FILTER (WHERE phenotype='S' AND has_exact_known_determinant) AS n_s_with_known,
          count(*) FILTER (WHERE phenotype='S' AND portal_unexplained) AS n_s_without_known,
          count(DISTINCT ISO_country_code) FILTER (WHERE phenotype='R' AND portal_unexplained AND ISO_country_code IS NOT NULL) AS n_countries_r_without_known,
          count(DISTINCT collection_year) FILTER (WHERE phenotype='R' AND portal_unexplained AND collection_year IS NOT NULL) AS n_years_r_without_known,
          min(collection_year) FILTER (WHERE phenotype='R' AND portal_unexplained) AS first_year_r_without_known,
          max(collection_year) FILTER (WHERE phenotype='R' AND portal_unexplained) AS last_year_r_without_known,
          count(DISTINCT AMR_associated_publications) FILTER (WHERE phenotype='R' AND portal_unexplained AND AMR_associated_publications IS NOT NULL) AS n_publication_strings_r_without_known
        FROM strict_labeled
        GROUP BY ALL
        HAVING n_r >= ? AND n_s >= ? AND n_r_without_known >= ?
        ORDER BY n_r_without_known DESC, n_s_without_known DESC
        """,
        [args.min_resistant, args.min_susceptible, args.min_unexplained_resistant],
    )
    if pair_stats.empty:
        raise RuntimeError("No organism–antibiotic pair passed the pre-specified minimum counts")

    scored_rows: list[dict[str, Any]] = []
    for _, row in pair_stats.iterrows():
        score, extra = clinical_score(row)
        d = row.to_dict()
        d.update(extra)
        d["selection_score"] = score
        scored_rows.append(d)
    ranking = pd.DataFrame(scored_rows).sort_values(
        ["high_impact_pair", "selection_score", "n_r_without_known"],
        ascending=[False, False, False],
    )
    ranking.to_csv(out / "pair_ranking.csv", index=False)

    eligible = ranking[ranking["high_impact_pair"]]
    if eligible.empty:
        eligible = ranking
    target = eligible.iloc[0].to_dict()
    target_key = [target["organism"], target["antibiotic_name"], target["antibiotic_ontology"]]

    target_rows = q(
        con,
        """
        SELECT *
        FROM strict_labeled
        WHERE organism=? AND antibiotic_name=? AND antibiotic_ontology=?
        ORDER BY phenotype, portal_unexplained DESC, collection_year, ISO_country_code, assembly_ID
        """,
        target_key,
    )
    target_rows.to_csv(out / "selected_pair_all_records.csv", index=False)

    # Round-robin country/year strata before the deterministic hash order, so a large
    # single study cannot automatically dominate the genome manifest.
    manifest = q(
        con,
        """
        WITH eligible AS (
          SELECT *,
                 row_number() OVER (
                   PARTITION BY phenotype,
                     coalesce(ISO_country_code, 'UNK'),
                     coalesce(collection_year, -1)
                   ORDER BY hash(assembly_ID, BioSample_ID)
                 ) AS stratum_rank
          FROM strict_labeled
          WHERE organism=? AND antibiotic_name=? AND antibiotic_ontology=?
            AND portal_unexplained
            AND phenotype IN ('R','S')
        ), ranked AS (
          SELECT *,
                 row_number() OVER (
                   PARTITION BY phenotype
                   ORDER BY stratum_rank, hash(assembly_ID, BioSample_ID)
                 ) AS class_rank
          FROM eligible
        )
        SELECT * EXCLUDE(stratum_rank, class_rank)
        FROM ranked
        WHERE class_rank <= ?
        ORDER BY phenotype, assembly_ID
        """,
        target_key + [args.manifest_per_class],
    )
    manifest.to_csv(out / "genome_manifest_stage2.csv", index=False)

    global_counts = q(
        con,
        f"""
        SELECT
          (SELECT count(*) FROM read_parquet('{phen}')) AS phenotype_rows,
          (SELECT count(DISTINCT BioSample_ID) FROM read_parquet('{phen}')) AS phenotype_biosamples,
          (SELECT count(DISTINCT assembly_ID) FROM read_parquet('{phen}') WHERE assembly_ID IS NOT NULL) AS phenotype_assemblies,
          (SELECT count(*) FROM phen_strict) AS strict_assembly_antibiotic_pairs,
          (SELECT count(DISTINCT assembly_ID) FROM phen_strict) AS strict_assemblies,
          (SELECT count(*) FROM read_parquet('{geno}')) AS genotype_rows,
          (SELECT count(DISTINCT assembly_ID) FROM read_parquet('{geno}')) AS genotype_assemblies,
          (SELECT count(*) FROM read_parquet('{merged}')) AS merged_rows
        """,
    ).iloc[0].to_dict()

    phenotype_qc = q(
        con,
        f"""
        WITH n AS (
          SELECT BioSample_ID, assembly_ID, antibiotic_ontology, antibiotic_name,
                 count(DISTINCT CASE
                   WHEN lower(trim(resistance_phenotype))='resistant' THEN 'R'
                   WHEN lower(trim(resistance_phenotype))='susceptible' THEN 'S'
                 END) AS n_labels
          FROM read_parquet('{phen}')
          WHERE assembly_ID IS NOT NULL AND antibiotic_ontology IS NOT NULL
          GROUP BY ALL
        )
        SELECT
          count(*) AS raw_keys,
          count(*) FILTER (WHERE n_labels=1) AS unambiguous_keys,
          count(*) FILTER (WHERE n_labels>1) AS conflicting_keys,
          count(*) FILTER (WHERE n_labels=0) AS keys_without_strict_r_or_s
        FROM n
        """,
    ).iloc[0].to_dict()

    release_info = {
        "release": release,
        "release_root": base,
        "phenotype_url": phen,
        "genotype_url": geno,
        "merged_url": merged,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "listing_sha256_not_recorded": True,
        "listing_excerpt": listing[:1000],
    }
    (out / "release.json").write_text(json.dumps(release_info, indent=2, default=str) + "\n")
    (out / "global_counts.json").write_text(json.dumps(global_counts, indent=2, default=str) + "\n")
    (out / "phenotype_qc.json").write_text(json.dumps(phenotype_qc, indent=2, default=str) + "\n")
    (out / "selected_target.json").write_text(json.dumps(target, indent=2, default=str) + "\n")

    top = ranking.head(20).copy()
    cols = [
        "organism", "antibiotic_name", "antibiotic_ontology", "n_r", "n_s",
        "n_r_with_known", "n_r_without_known", "n_s_without_known",
        "n_countries_r_without_known", "n_years_r_without_known",
        "r_portal_unexplained_fraction", "ontology_mismatch_penalty", "selection_score",
    ]
    report = [
        "# Global public-AMR discordance census",
        "",
        f"- EMBL-EBI AMR Portal release: **{release}**",
        f"- Strict assembly–antibiotic phenotype pairs: **{int(global_counts['strict_assembly_antibiotic_pairs']):,}**",
        f"- Assemblies represented by strict R/S phenotypes: **{int(global_counts['strict_assemblies']):,}**",
        f"- Candidate genome manifest: **{len(manifest):,}** assemblies",
        "",
        "## Automatically selected high-impact pair",
        "",
        f"- Organism: **{target['organism']}**",
        f"- Antibiotic: **{target['antibiotic_name']}** (`{target['antibiotic_ontology']}`)",
        f"- Resistant / susceptible: **{int(target['n_r']):,} / {int(target['n_s']):,}**",
        f"- Resistant without an exact same-antibiotic portal genotype: **{int(target['n_r_without_known']):,}**",
        f"- Susceptible without an exact same-antibiotic portal genotype: **{int(target['n_s_without_known']):,}**",
        f"- Resistant portal-unexplained fraction: **{float(target['r_portal_unexplained_fraction']):.3f}**",
        f"- Countries / collection years among portal-unexplained resistant isolates: **{int(target['n_countries_r_without_known'])} / {int(target['n_years_r_without_known'])}**",
        "",
        "## Claim boundary",
        "",
        "The anti-join identifies resistance not explained by an **exact antibiotic-ontology match** in the portal's AMRFinderPlus genotype view. It does not prove absence of all known biological mechanisms. The selected pair advances only to a population-structure-aware genome analysis and known-mechanism re-audit.",
        "",
        "## Top ranked pairs",
        "",
        top[cols].to_markdown(index=False),
        "",
    ]
    (out / "PROFILE_REPORT.md").write_text("\n".join(report))

    provenance = {
        "python": __import__("sys").version,
        "duckdb": duckdb.__version__,
        "pandas": pd.__version__,
        "requests": requests.__version__,
        "arguments": vars(args),
        "operational_definition": "strict R phenotype and no exact same-antibiotic ontology match in AMR Portal merged genotype view",
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")
    con.close()
    print((out / "PROFILE_REPORT.md").read_text())


if __name__ == "__main__":
    main()
