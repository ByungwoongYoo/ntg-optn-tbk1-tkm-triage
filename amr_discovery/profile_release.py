#!/usr/bin/env python3
"""Profile the EMBL-EBI AMR Portal release and rank organism-drug pairs
with phenotypic resistance not matched to a known AMRFinderPlus determinant.

This is a discovery triage only. "Unexplained" means no genotype row joined on
BioSample_ID, assembly_ID and antibiotic_ontology in the same data release; it
does not mean that no biological mechanism is known.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import urllib.request

import duckdb
import pandas as pd

RELEASE = "2026-07"
BASE = f"https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/{RELEASE}"
PHENO = f"{BASE}/phenotype.parquet"
GENO = f"{BASE}/genotype.parquet"
MERGED = f"{BASE}/phenotype_genotype_merged.parquet"
OUT = pathlib.Path("amr_profile")
OUT.mkdir(exist_ok=True)


def fetch(url: str, path: pathlib.Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "AMR-unexplained-resistance-audit/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(path, "wb") as f:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


# Download the compact Parquet files so the exact audited bytes are preserved.
for name in ["phenotype.parquet", "genotype.parquet", "phenotype_genotype_merged.parquet", "manifest.yaml", "md5sums"]:
    fetch(f"{BASE}/{name}", OUT / name)

con = duckdb.connect()
con.execute("PRAGMA threads=4")

p = str((OUT / "phenotype.parquet").resolve())
g = str((OUT / "genotype.parquet").resolve())
m = str((OUT / "phenotype_genotype_merged.parquet").resolve())

# Raw release dimensions.
summary = {
    "release": RELEASE,
    "source": BASE,
    "files": {x.name: {"bytes": x.stat().st_size, "sha256": sha256(x)} for x in OUT.iterdir() if x.is_file()},
}
summary["phenotype_rows"] = con.execute(f"SELECT count(*) FROM read_parquet('{p}')").fetchone()[0]
summary["phenotype_biosamples"] = con.execute(f"SELECT count(DISTINCT BioSample_ID) FROM read_parquet('{p}')").fetchone()[0]
summary["phenotype_assemblies"] = con.execute(f"SELECT count(DISTINCT assembly_ID) FROM read_parquet('{p}') WHERE assembly_ID IS NOT NULL").fetchone()[0]
summary["genotype_rows"] = con.execute(f"SELECT count(*) FROM read_parquet('{g}')").fetchone()[0]
summary["genotype_biosamples"] = con.execute(f"SELECT count(DISTINCT BioSample_ID) FROM read_parquet('{g}')").fetchone()[0]
summary["merged_rows"] = con.execute(f"SELECT count(*) FROM read_parquet('{m}')").fetchone()[0]

# Preserve all phenotype spellings before normalization.
values = con.execute(f"""
    SELECT lower(trim(resistance_phenotype)) AS phenotype, count(*) AS n
    FROM read_parquet('{p}')
    GROUP BY 1 ORDER BY n DESC
""").df()
values.to_csv(OUT / "phenotype_value_counts.csv", index=False)

# One non-conflicting phenotype per sample/assembly/organism/antibiotic.
# R-class: resistant, non-susceptible/nonsusceptible; S-class: susceptible/sensitive.
# Intermediate and SDD are excluded from the primary binary screen.
con.execute(f"""
CREATE OR REPLACE TEMP VIEW pheno_binary AS
WITH x AS (
  SELECT
    BioSample_ID,
    assembly_ID,
    organism,
    species,
    antibiotic_name,
    antibiotic_ontology,
    collection_year,
    country,
    geographical_region,
    CASE
      WHEN lower(trim(resistance_phenotype)) IN ('resistant','r','non-susceptible','nonsusceptible','ns') THEN 'R'
      WHEN lower(trim(resistance_phenotype)) IN ('susceptible','sensitive','s') THEN 'S'
      ELSE NULL
    END AS y
  FROM read_parquet('{p}')
  WHERE assembly_ID IS NOT NULL
    AND antibiotic_ontology IS NOT NULL
), d AS (
  SELECT *,
         count(DISTINCT y) OVER (PARTITION BY BioSample_ID, assembly_ID, organism, antibiotic_ontology) AS n_labels,
         row_number() OVER (PARTITION BY BioSample_ID, assembly_ID, organism, antibiotic_ontology ORDER BY collection_year NULLS LAST) AS rn
  FROM x WHERE y IS NOT NULL
)
SELECT BioSample_ID, assembly_ID, organism, species, antibiotic_name, antibiotic_ontology,
       collection_year, country, geographical_region, y
FROM d WHERE n_labels=1 AND rn=1
""")

# Known AMRFinderPlus match for the same ontology in the same release.
con.execute(f"""
CREATE OR REPLACE TEMP VIEW known_same_drug AS
SELECT DISTINCT BioSample_ID, assembly_ID, antibiotic_ontology
FROM read_parquet('{g}')
WHERE antibiotic_ontology IS NOT NULL
""")

pair = con.execute("""
WITH j AS (
  SELECT p.*,
         CASE WHEN k.BioSample_ID IS NULL THEN 0 ELSE 1 END AS known_same_drug
  FROM pheno_binary p
  LEFT JOIN known_same_drug k USING (BioSample_ID, assembly_ID, antibiotic_ontology)
)
SELECT organism, antibiotic_name, antibiotic_ontology,
       count(*) AS n_total,
       count(*) FILTER (WHERE y='R') AS n_R,
       count(*) FILTER (WHERE y='S') AS n_S,
       count(*) FILTER (WHERE y='R' AND known_same_drug=0) AS n_R_unexplained,
       count(*) FILTER (WHERE y='R' AND known_same_drug=1) AS n_R_with_known,
       round(100.0 * count(*) FILTER (WHERE y='R' AND known_same_drug=0) /
             NULLIF(count(*) FILTER (WHERE y='R'),0), 2) AS pct_R_unexplained,
       count(DISTINCT country) AS n_countries,
       min(collection_year) AS first_year,
       max(collection_year) AS last_year
FROM j
GROUP BY ALL
HAVING n_R >= 20 AND n_S >= 20
ORDER BY n_R_unexplained DESC, n_total DESC
""").df()
pair.to_csv(OUT / "organism_antibiotic_screen.csv", index=False)

# A stricter shortlist for downstream whole-genome association.
short = pair[(pair.n_R >= 80) & (pair.n_S >= 80) & (pair.n_R_unexplained >= 30)].copy()
short["balance"] = short[["n_R", "n_S"]].min(axis=1) / short[["n_R", "n_S"]].max(axis=1)
short["triage_score"] = (
    short["n_R_unexplained"].clip(upper=1000).pow(0.5)
    * short["balance"].clip(lower=0.05).pow(0.5)
    * (short["n_countries"].clip(lower=1).pow(0.25))
)
short = short.sort_values(["triage_score", "n_R_unexplained"], ascending=False)
short.to_csv(OUT / "downstream_shortlist.csv", index=False)

# Export isolate-level labels for the top 20 pairs.
top = short.head(20)[["organism", "antibiotic_ontology"]]
con.register("top_pairs", top)
labels = con.execute("""
WITH j AS (
  SELECT p.*, CASE WHEN k.BioSample_ID IS NULL THEN 0 ELSE 1 END AS known_same_drug
  FROM pheno_binary p
  LEFT JOIN known_same_drug k USING (BioSample_ID, assembly_ID, antibiotic_ontology)
)
SELECT j.* FROM j
JOIN top_pairs t USING (organism, antibiotic_ontology)
ORDER BY organism, antibiotic_ontology, y, BioSample_ID
""").df()
labels.to_csv(OUT / "top_pair_isolate_labels.csv", index=False)

summary["binary_nonconflicting_rows"] = int(con.execute("SELECT count(*) FROM pheno_binary").fetchone()[0])
summary["screened_pairs_n20_each"] = int(len(pair))
summary["strict_shortlist_count"] = int(len(short))
summary["top_shortlist"] = short.head(20).to_dict(orient="records")
(OUT / "PROFILE_SUMMARY.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

md = [
    f"# EMBL-EBI AMR Portal {RELEASE} profile",
    "",
    "## Release scale",
    f"- Phenotype rows: **{summary['phenotype_rows']:,}**",
    f"- Phenotyped BioSamples: **{summary['phenotype_biosamples']:,}**",
    f"- Phenotyped assemblies: **{summary['phenotype_assemblies']:,}**",
    f"- Genotype rows: **{summary['genotype_rows']:,}**",
    f"- Genotyped BioSamples: **{summary['genotype_biosamples']:,}**",
    f"- Non-conflicting binary sample-drug labels: **{summary['binary_nonconflicting_rows']:,}**",
    "",
    "## Strict downstream shortlist",
    "Criteria: R>=80, S>=80, phenotypically resistant with no same-drug AMRFinderPlus match>=30.",
    "",
]
if len(short):
    md.append(short.head(20).to_markdown(index=False))
else:
    md.append("No pair met the strict prespecified screen.")
md += [
    "",
    "## Interpretation boundary",
    "`n_R_unexplained` means no AMRFinderPlus genotype row joined to the same antibiotic ontology in this release. It is a triage label, not proof that the resistance mechanism is biologically unknown.",
]
(OUT / "PROFILE_REPORT.md").write_text("\n".join(md) + "\n")

# Hash all outputs last, excluding the hash file itself.
with open(OUT / "SHA256SUMS.txt", "w") as f:
    for path in sorted(OUT.iterdir()):
        if path.name != "SHA256SUMS.txt" and path.is_file():
            f.write(f"{sha256(path)}  {path.name}\n")
print((OUT / "PROFILE_REPORT.md").read_text())
