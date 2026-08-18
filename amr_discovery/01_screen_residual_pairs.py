#!/usr/bin/env python3
"""Screen a pinned public AMR release for phenotype–genotype residuals.

This stage identifies organism–antibiotic pairs suitable for whole-genome marker
discovery. It does not claim a novel marker or a biologically unexplained phenotype.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Any

import duckdb

R_STRICT = ("resistant", "r")
R_BROAD = R_STRICT + ("non-susceptible", "nonsusceptible", "non susceptible", "ns")
S_VALUES = ("susceptible", "sensitive", "s")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--release", default="2026-07")
    p.add_argument("--data-dir", default="amr_data")
    p.add_argument("--out-dir", default="amr_results/screen")
    p.add_argument("--min-total", type=int, default=200)
    p.add_argument("--min-resistant", type=int, default=50)
    p.add_argument("--min-susceptible", type=int, default=50)
    p.add_argument("--min-residual-resistant", type=int, default=25)
    p.add_argument("--min-known-resistant", type=int, default=10)
    p.add_argument("--min-known-or", type=float, default=2.0)
    p.add_argument("--top-n", type=int, default=30)
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sql_list(values: tuple[str, ...]) -> str:
    return ",".join("'" + v.replace("'", "''") + "'" for v in values)


def qstr(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def relation_to_dicts(rel: duckdb.DuckDBPyRelation) -> list[dict[str, Any]]:
    cols = [d[0] for d in rel.description]
    return [dict(zip(cols, row)) for row in rel.fetchall()]


def copy_query(con: duckdb.DuckDBPyConnection, query: str, path: Path) -> None:
    con.execute(f"COPY ({query}) TO {qstr(path)} (HEADER, DELIMITER '\\t')")


def main() -> int:
    args = parse_args()
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    phenotype = data_dir / "phenotype.parquet"
    genotype = data_dir / "genotype.parquet"
    merged = data_dir / "phenotype_genotype_merged.parquet"
    for path in (phenotype, genotype, merged):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required file: {path}")

    con = duckdb.connect(str(out_dir / "screen.duckdb"))
    con.execute("PRAGMA threads=4")
    con.execute("PRAGMA memory_limit='10GB'")
    con.execute("PRAGMA preserve_insertion_order=false")
    con.execute(f"CREATE OR REPLACE VIEW phenotype_raw AS SELECT * FROM read_parquet({qstr(phenotype)})")
    con.execute(f"CREATE OR REPLACE VIEW genotype_raw AS SELECT * FROM read_parquet({qstr(genotype)})")
    con.execute(f"CREATE OR REPLACE VIEW merged_raw AS SELECT * FROM read_parquet({qstr(merged)})")

    census = relation_to_dicts(con.sql("""
        SELECT
          (SELECT count(*) FROM phenotype_raw) AS phenotype_rows,
          (SELECT count(DISTINCT BioSample_ID) FROM phenotype_raw) AS phenotype_biosamples,
          (SELECT count(DISTINCT assembly_ID) FROM phenotype_raw WHERE assembly_ID IS NOT NULL) AS phenotype_assemblies,
          (SELECT count(*) FROM genotype_raw) AS genotype_rows,
          (SELECT count(DISTINCT BioSample_ID) FROM genotype_raw) AS genotype_biosamples,
          (SELECT count(DISTINCT assembly_ID) FROM genotype_raw) AS genotype_assemblies,
          (SELECT count(*) FROM merged_raw) AS merged_rows,
          (SELECT count(DISTINCT BioSample_ID) FROM merged_raw) AS merged_biosamples,
          (SELECT count(DISTINCT assembly_ID) FROM merged_raw) AS merged_assemblies
    """))[0]
    census.update({
        "release": args.release,
        "input_files": {
            x.name: {"bytes": x.stat().st_size, "sha256": sha256_file(x)}
            for x in (phenotype, genotype, merged)
        },
        "runtime": {
            "python": sys.version,
            "duckdb": duckdb.__version__,
            "platform": platform.platform(),
        },
    })
    write_json(out_dir / "release_census.json", census)

    copy_query(con, """
        SELECT lower(trim(coalesce(resistance_phenotype, '<NULL>'))) AS raw_value,
               count(*) AS n
        FROM phenotype_raw
        GROUP BY 1 ORDER BY n DESC
    """, out_dir / "phenotype_value_counts.tsv")

    # Presence in genotype.parquet proves at least one AMRFinderPlus result exists
    # for the assembly. This is a deliberately conservative calibration subset.
    con.execute("""
        CREATE OR REPLACE TABLE genotyped_assembly AS
        SELECT DISTINCT BioSample_ID, assembly_ID
        FROM genotype_raw
        WHERE assembly_ID IS NOT NULL AND BioSample_ID IS NOT NULL
    """)
    con.execute("""
        CREATE OR REPLACE TABLE known_determinant AS
        SELECT DISTINCT BioSample_ID, assembly_ID, antibiotic_ontology
        FROM genotype_raw
        WHERE assembly_ID IS NOT NULL
          AND BioSample_ID IS NOT NULL
          AND antibiotic_ontology IS NOT NULL
          AND trim(antibiotic_ontology) <> ''
    """)

    def build_screen(label: str, resistant_values: tuple[str, ...]) -> None:
        calls_table = f"calls_{label}"
        pair_table = f"pairs_{label}"
        rvals = sql_list(resistant_values)
        svals = sql_list(S_VALUES)

        con.execute(f"""
            CREATE OR REPLACE TABLE {calls_table} AS
            WITH mapped AS (
              SELECT
                p.BioSample_ID,
                p.assembly_ID,
                p.species,
                p.organism,
                p.antibiotic_ontology,
                p.antibiotic_name,
                p.collection_year,
                p.ISO_country_code,
                p.country,
                p.geographical_region,
                p.host,
                p.isolation_source_category,
                p.AMR_associated_publications,
                CASE
                  WHEN lower(trim(p.resistance_phenotype)) IN ({rvals}) THEN 'R'
                  WHEN lower(trim(p.resistance_phenotype)) IN ({svals}) THEN 'S'
                  ELSE NULL
                END AS call
              FROM phenotype_raw p
              INNER JOIN genotyped_assembly ga
                ON p.BioSample_ID = ga.BioSample_ID
               AND p.assembly_ID = ga.assembly_ID
              WHERE p.assembly_ID IS NOT NULL
                AND p.BioSample_ID IS NOT NULL
                AND p.species IS NOT NULL
                AND p.antibiotic_ontology IS NOT NULL
                AND trim(p.antibiotic_ontology) <> ''
            ), grouped AS (
              SELECT
                BioSample_ID, assembly_ID, species,
                any_value(organism) AS organism,
                antibiotic_ontology,
                any_value(antibiotic_name) AS antibiotic_name,
                any_value(collection_year) AS collection_year,
                any_value(ISO_country_code) AS ISO_country_code,
                any_value(country) AS country,
                any_value(geographical_region) AS geographical_region,
                any_value(host) AS host,
                any_value(isolation_source_category) AS isolation_source_category,
                any_value(AMR_associated_publications) AS AMR_associated_publications,
                count(*) AS source_rows,
                count(DISTINCT call) FILTER (WHERE call IS NOT NULL) AS n_distinct_calls,
                min(call) FILTER (WHERE call IS NOT NULL) AS consensus_call
              FROM mapped
              GROUP BY BioSample_ID, assembly_ID, species, antibiotic_ontology
            )
            SELECT
              g.*,
              CASE WHEN kd.BioSample_ID IS NULL THEN false ELSE true END AS has_known_determinant
            FROM grouped g
            LEFT JOIN known_determinant kd
              ON g.BioSample_ID = kd.BioSample_ID
             AND g.assembly_ID = kd.assembly_ID
             AND g.antibiotic_ontology = kd.antibiotic_ontology
            WHERE g.n_distinct_calls = 1
              AND g.consensus_call IS NOT NULL
        """)

        con.execute(f"""
            CREATE OR REPLACE TABLE {pair_table} AS
            WITH agg AS (
              SELECT
                species,
                any_value(organism) AS organism,
                antibiotic_ontology,
                any_value(antibiotic_name) AS antibiotic_name,
                count(*) AS n_total,
                count(*) FILTER (WHERE consensus_call='R') AS n_R,
                count(*) FILTER (WHERE consensus_call='S') AS n_S,
                count(*) FILTER (WHERE consensus_call='R' AND has_known_determinant) AS a_R_known,
                count(*) FILTER (WHERE consensus_call='S' AND has_known_determinant) AS b_S_known,
                count(*) FILTER (WHERE consensus_call='R' AND NOT has_known_determinant) AS c_R_residual,
                count(*) FILTER (WHERE consensus_call='S' AND NOT has_known_determinant) AS d_S_no_known,
                count(DISTINCT ISO_country_code) FILTER (WHERE ISO_country_code IS NOT NULL) AS countries,
                count(DISTINCT collection_year) FILTER (WHERE collection_year IS NOT NULL) AS years,
                min(collection_year) FILTER (WHERE collection_year IS NOT NULL) AS first_year,
                max(collection_year) FILTER (WHERE collection_year IS NOT NULL) AS last_year,
                count(DISTINCT geographical_region) FILTER (WHERE geographical_region IS NOT NULL) AS regions,
                count(DISTINCT AMR_associated_publications) FILTER (
                  WHERE AMR_associated_publications IS NOT NULL
                    AND trim(AMR_associated_publications)<>''
                ) AS publication_strings,
                count(*) FILTER (
                  WHERE lower(coalesce(host,'')) LIKE '%homo sapiens%'
                     OR lower(coalesce(host,''))='human'
                ) AS human_host_n
              FROM {calls_table}
              GROUP BY species, antibiotic_ontology
            )
            SELECT
              *,
              c_R_residual::DOUBLE / nullif(n_R,0) AS residual_fraction_R,
              a_R_known::DOUBLE / nullif(n_R,0) AS known_explanation_fraction_R,
              ((a_R_known + 0.5) * (d_S_no_known + 0.5)) /
                ((b_S_known + 0.5) * (c_R_residual + 0.5)) AS known_determinant_OR,
              a_R_known::DOUBLE / nullif(a_R_known + c_R_residual,0) AS known_determinant_sensitivity,
              d_S_no_known::DOUBLE / nullif(b_S_known + d_S_no_known,0) AS known_determinant_specificity,
              least(n_R,n_S)::DOUBLE / greatest(n_R,n_S) AS class_balance,
              sqrt(c_R_residual::DOUBLE * n_S::DOUBLE)
                * log2(2 + countries)
                * log2(2 + years)
                * sqrt(least(n_R,n_S)::DOUBLE / greatest(n_R,n_S)) AS discovery_score
            FROM agg
            WHERE n_total >= {args.min_total}
              AND n_R >= {args.min_resistant}
              AND n_S >= {args.min_susceptible}
              AND c_R_residual >= {args.min_residual_resistant}
              AND a_R_known >= {args.min_known_resistant}
              AND ((a_R_known + 0.5) * (d_S_no_known + 0.5)) /
                  ((b_S_known + 0.5) * (c_R_residual + 0.5)) >= {args.min_known_or}
            ORDER BY discovery_score DESC
        """)

        copy_query(con, f"SELECT * FROM {pair_table}", out_dir / f"screened_pairs_{label}.tsv")
        copy_query(con, f"""
            SELECT c.*
            FROM {calls_table} c
            INNER JOIN (
              SELECT species, antibiotic_ontology
              FROM {pair_table}
              ORDER BY discovery_score DESC
              LIMIT {args.top_n}
            ) t USING (species, antibiotic_ontology)
            ORDER BY c.species, c.antibiotic_ontology, c.consensus_call, c.assembly_ID
        """, out_dir / f"top_pair_samples_{label}.tsv")

    build_screen("strict", R_STRICT)
    build_screen("broad", R_BROAD)

    top = relation_to_dicts(con.sql(
        "SELECT * FROM pairs_strict ORDER BY discovery_score DESC LIMIT 1"
    ))
    selection_mode = "strict"
    if not top:
        top = relation_to_dicts(con.sql(
            "SELECT * FROM pairs_broad ORDER BY discovery_score DESC LIMIT 1"
        ))
        selection_mode = "broad_fallback"

    top_obj = {
        "selection_mode": selection_mode,
        "selection_gates": {
            "min_total": args.min_total,
            "min_resistant": args.min_resistant,
            "min_susceptible": args.min_susceptible,
            "min_residual_resistant": args.min_residual_resistant,
            "min_known_resistant": args.min_known_resistant,
            "min_known_determinant_OR": args.min_known_or,
        },
        "selected_pair": top[0] if top else None,
        "claim_boundary": (
            "This screen identifies a candidate phenotype-genotype residual set. "
            "It does not establish that the phenotype is biologically unexplained, "
            "nor that an uncatalogued causal marker exists."
        ),
    }
    write_json(out_dir / "TOP_PAIR.json", top_obj)

    copy_query(con, """
        SELECT species, count(*) AS n_calls,
               count(DISTINCT antibiotic_ontology) AS antibiotics,
               count(DISTINCT assembly_ID) AS assemblies
        FROM calls_strict
        GROUP BY species
        ORDER BY n_calls DESC
    """, out_dir / "species_coverage_strict.tsv")

    methods = f"""# Stage 1: residual-pair screen

Pinned release: EMBL-EBI AMR Portal `{args.release}`.

Duplicate phenotype rows are collapsed to one unambiguous R/S call per BioSample–assembly–antibiotic ontology. The conservative screen requires the assembly to occur in the genotype release and treats a determinant as catalogued only when the genotype table contains the same antibiotic ontology.

Candidate pairs require at least {args.min_total} labeled assemblies, {args.min_resistant} resistant, {args.min_susceptible} susceptible, {args.min_residual_resistant} resistant without an exact ontology-matched determinant, and {args.min_known_resistant} resistant with a matched determinant. The known determinant must itself associate with resistance at OR >= {args.min_known_or}.

This is a triage output, not a causal or novelty claim. Marker discovery still requires whole-genome features, lineage correction, temporal/geographic validation, and database/literature exclusion.
"""
    (out_dir / "METHODS.md").write_text(methods, encoding="utf-8")

    manifest: list[dict[str, Any]] = []
    for path in sorted(out_dir.iterdir()):
        if path.is_file() and path.name not in {
            "screen.duckdb", "SHA256SUMS.txt", "manifest.json"
        }:
            manifest.append({
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    write_json(out_dir / "manifest.json", {"release": args.release, "files": manifest})
    (out_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{x['sha256']}  {x['path']}\n" for x in manifest),
        encoding="utf-8",
    )

    print(json.dumps(top_obj, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
