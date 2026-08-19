#!/usr/bin/env python3
"""Rank organism–drug pairs after broad AMRFinderPlus mechanism exclusion.

The initial exact-antibiotic ontology join is intentionally not used as a
novelty claim. This stage maps each drug to the broad AMRFinderPlus class and,
where appropriate, drug-specific subclass/symbol patterns. The primary
`R_no_relevant_known` stratum additionally requires at least one genotype row,
so missing portal annotation is not mistaken for unexplained biology.
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
OUT = pathlib.Path("amr_target_ranking")
OUT.mkdir(exist_ok=True)

for name in ["phenotype.parquet", "genotype.parquet", "manifest.yaml", "md5sums"]:
    path = OUT / name
    if not path.exists():
        req = urllib.request.Request(f"{BASE}/{name}", headers={"User-Agent": "AMR-target-ranking/1.0"})
        with urllib.request.urlopen(req, timeout=180) as r, open(path, "wb") as f:
            while chunk := r.read(1024 * 1024):
                f.write(chunk)

p = str((OUT / "phenotype.parquet").resolve())
g = str((OUT / "genotype.parquet").resolve())
con = duckdb.connect()
con.execute("PRAGMA threads=4")

# Prespecified candidate pairs chosen for sample size, public-health relevance,
# and the existence of both resistant and susceptible assemblies.
pairs = pd.DataFrame([
    ("Klebsiella pneumoniae", "tobramycin", "AMINOGLYCOSIDE"),
    ("Klebsiella pneumoniae", "gentamicin", "AMINOGLYCOSIDE"),
    ("Klebsiella pneumoniae", "amikacin", "AMINOGLYCOSIDE"),
    ("Acinetobacter baumannii", "tobramycin", "AMINOGLYCOSIDE"),
    ("Acinetobacter baumannii", "gentamicin", "AMINOGLYCOSIDE"),
    ("Acinetobacter baumannii", "amikacin", "AMINOGLYCOSIDE"),
    ("Pseudomonas aeruginosa", "tobramycin", "AMINOGLYCOSIDE"),
    ("Pseudomonas aeruginosa", "gentamicin", "AMINOGLYCOSIDE"),
    ("Klebsiella pneumoniae", "colistin", "COLISTIN"),
    ("Acinetobacter baumannii", "colistin", "COLISTIN"),
    ("Pseudomonas aeruginosa", "colistin", "COLISTIN"),
    ("Neisseria gonorrhoeae", "azithromycin", "MACROLIDE_AZITHROMYCIN"),
    ("Neisseria gonorrhoeae", "ciprofloxacin", "QUINOLONE"),
    ("Escherichia coli", "trimethoprim", "TRIMETHOPRIM"),
    ("Escherichia coli", "ampicillin", "BETA_LACTAM"),
    ("Klebsiella pneumoniae", "meropenem", "BETA_LACTAM"),
    ("Acinetobacter baumannii", "imipenem", "BETA_LACTAM"),
], columns=["organism", "antibiotic_name", "target_class"])
con.register("candidate_pairs", pairs)

con.execute(f"""
CREATE OR REPLACE TEMP VIEW phenotype_binary AS
WITH raw AS (
  SELECT
    ph.*,
    lower(trim(ph.antibiotic_name)) AS drug,
    CASE
      WHEN lower(trim(resistance_phenotype)) IN ('resistant','r','non-susceptible','nonsusceptible','ns') THEN 'R'
      WHEN lower(trim(resistance_phenotype)) IN ('susceptible','sensitive','s') THEN 'S'
      ELSE NULL
    END AS y
  FROM read_parquet('{p}') ph
  JOIN candidate_pairs cp
    ON ph.organism=cp.organism AND lower(trim(ph.antibiotic_name))=cp.antibiotic_name
  WHERE ph.assembly_ID IS NOT NULL
), labels AS (
  SELECT BioSample_ID, assembly_ID, organism, drug,
         count(DISTINCT y) FILTER (WHERE y IS NOT NULL) AS n_binary_labels
  FROM raw GROUP BY ALL
), ranked AS (
  SELECT r.*,
         row_number() OVER (
           PARTITION BY r.BioSample_ID,r.assembly_ID,r.organism,r.drug
           ORDER BY
             CASE WHEN r.measurement IS NOT NULL THEN 0 ELSE 1 END,
             CASE WHEN r.ast_standard IS NOT NULL THEN 0 ELSE 1 END,
             r.collection_year DESC NULLS LAST
         ) rn
  FROM raw r JOIN labels l USING(BioSample_ID,assembly_ID,organism,drug)
  WHERE r.y IS NOT NULL AND l.n_binary_labels=1
)
SELECT BioSample_ID, assembly_ID, organism, drug AS antibiotic_name, y,
       collection_year, country, ISO_country_code, geographical_region,
       host, isolation_source, isolation_source_category,
       AMR_associated_publications, ast_standard, laboratory_typing_method,
       measurement, measurement_sign, measurement_units, SRA_accession
FROM ranked WHERE rn=1
""")

# Flag relevant known determinants. Mixed class strings are matched by
# substring because AMRFinderPlus may emit combinations separated by '/'.
con.execute(f"""
CREATE OR REPLACE TEMP VIEW genotype_pair_flags AS
SELECT p.organism, p.antibiotic_name, p.target_class,
       g.BioSample_ID, g.assembly_ID,
       count(*) AS genotype_rows,
       max(CASE
         WHEN p.target_class='AMINOGLYCOSIDE' AND (
           upper(coalesce(g.class,'')) LIKE '%AMINOGLYCOSIDE%'
           OR upper(coalesce(g.subclass,'')) LIKE '%' || upper(p.antibiotic_name) || '%'
           OR upper(coalesce(g.split_subclass,''))=upper(p.antibiotic_name)
         ) THEN 1
         WHEN p.target_class='COLISTIN' AND (
           upper(coalesce(g.class,'')) LIKE '%COLISTIN%'
           OR upper(coalesce(g.subclass,'')) LIKE '%COLISTIN%'
           OR upper(coalesce(g.split_subclass,''))='COLISTIN'
         ) THEN 1
         WHEN p.target_class='MACROLIDE_AZITHROMYCIN' AND (
           upper(coalesce(g.class,'')) LIKE '%MACROLIDE%'
           OR upper(coalesce(g.subclass,'')) LIKE '%AZITHROMYCIN%'
           OR upper(coalesce(g.split_subclass,''))='AZITHROMYCIN'
           OR lower(coalesce(g.amr_element_symbol,'')) SIMILAR TO '(mtrr|mtrr_%|23s%|erm%|mef%|msr%|mph%|ere%)'
         ) THEN 1
         WHEN p.target_class='QUINOLONE' AND (
           upper(coalesce(g.class,'')) LIKE '%QUINOLONE%'
           OR upper(coalesce(g.class,'')) LIKE '%FLUOROQUINOLONE%'
         ) THEN 1
         WHEN p.target_class='TRIMETHOPRIM' AND upper(coalesce(g.class,'')) LIKE '%TRIMETHOPRIM%' THEN 1
         WHEN p.target_class='BETA_LACTAM' AND upper(coalesce(g.class,'')) LIKE '%BETA-LACTAM%' THEN 1
         ELSE 0 END) AS has_relevant_known,
       string_agg(DISTINCT coalesce(g.amr_element_symbol,''), ';') AS all_elements,
       string_agg(DISTINCT coalesce(g.class,''), ';') AS all_classes,
       string_agg(DISTINCT coalesce(g.subclass,''), ';') AS all_subclasses
FROM candidate_pairs p
JOIN read_parquet('{g}') g ON g.organism=p.organism
GROUP BY ALL
""")

labels = con.execute("""
SELECT ph.*, cp.target_class,
       coalesce(gf.genotype_rows,0) AS genotype_rows,
       coalesce(gf.has_relevant_known,0) AS has_relevant_known,
       coalesce(gf.all_elements,'') AS all_elements,
       coalesce(gf.all_classes,'') AS all_classes,
       coalesce(gf.all_subclasses,'') AS all_subclasses,
       CASE WHEN coalesce(gf.genotype_rows,0)>0 AND coalesce(gf.has_relevant_known,0)=0 THEN 1 ELSE 0 END AS no_relevant_known_annotated
FROM phenotype_binary ph
JOIN candidate_pairs cp USING(organism,antibiotic_name)
LEFT JOIN genotype_pair_flags gf USING(organism,antibiotic_name,target_class,BioSample_ID,assembly_ID)
ORDER BY organism,antibiotic_name,y,BioSample_ID
""").df()
labels.to_csv(OUT / "candidate_pair_isolate_labels.csv", index=False)

rank = con.execute("""
WITH x AS (
  SELECT ph.*, cp.target_class,
         coalesce(gf.genotype_rows,0) AS genotype_rows,
         coalesce(gf.has_relevant_known,0) AS has_relevant_known,
         CASE WHEN coalesce(gf.genotype_rows,0)>0 AND coalesce(gf.has_relevant_known,0)=0 THEN 1 ELSE 0 END AS no_relevant_known_annotated
  FROM phenotype_binary ph
  JOIN candidate_pairs cp USING(organism,antibiotic_name)
  LEFT JOIN genotype_pair_flags gf USING(organism,antibiotic_name,target_class,BioSample_ID,assembly_ID)
)
SELECT organism, antibiotic_name, target_class,
       count(*) AS n_total,
       count(*) FILTER (WHERE y='R') AS n_R,
       count(*) FILTER (WHERE y='S') AS n_S,
       count(*) FILTER (WHERE y='R' AND genotype_rows=0) AS n_R_missing_all_genotype,
       count(*) FILTER (WHERE y='R' AND has_relevant_known=1) AS n_R_relevant_known,
       count(*) FILTER (WHERE y='R' AND no_relevant_known_annotated=1) AS n_R_no_relevant_known,
       count(*) FILTER (WHERE y='S' AND no_relevant_known_annotated=1) AS n_S_no_relevant_known,
       count(DISTINCT country) FILTER (WHERE y='R' AND no_relevant_known_annotated=1) AS unexplained_R_countries,
       min(collection_year) FILTER (WHERE y='R' AND no_relevant_known_annotated=1) AS unexplained_R_first_year,
       max(collection_year) FILTER (WHERE y='R' AND no_relevant_known_annotated=1) AS unexplained_R_last_year,
       count(*) FILTER (WHERE y='R' AND no_relevant_known_annotated=1 AND measurement IS NOT NULL) AS unexplained_R_with_measurement,
       count(*) FILTER (WHERE y='R' AND no_relevant_known_annotated=1 AND ast_standard IS NOT NULL) AS unexplained_R_with_standard
FROM x GROUP BY ALL
ORDER BY n_R_no_relevant_known DESC, n_S_no_relevant_known DESC
""").df()

rank["balance_no_known"] = rank[["n_R_no_relevant_known","n_S_no_relevant_known"]].min(axis=1) / rank[["n_R_no_relevant_known","n_S_no_relevant_known"]].max(axis=1).replace(0, pd.NA)
rank["evidence_score"] = (
    rank["n_R_no_relevant_known"].clip(upper=1000).pow(0.5)
    * rank["n_S_no_relevant_known"].clip(upper=2000).pow(0.25)
    * rank["balance_no_known"].fillna(0).clip(lower=0.02).pow(0.25)
    * rank["unexplained_R_countries"].clip(lower=1).pow(0.15)
)
rank = rank.sort_values(["evidence_score","n_R_no_relevant_known"], ascending=False)
rank.to_csv(OUT / "target_pair_ranking.csv", index=False)

summary = {
    "release": RELEASE,
    "operational_definition": "R with >=1 genotype row and no prespecified relevant AMRFinderPlus class/subclass/symbol match",
    "candidate_pairs": len(pairs),
    "binary_nonconflicting_labels": int(len(labels)),
    "ranking": rank.to_dict(orient="records"),
}
(OUT / "TARGET_RANKING_SUMMARY.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

md = [
    f"# Mechanistically filtered AMR target ranking ({RELEASE})",
    "",
    "Primary operational stratum: phenotypically resistant assembly with at least one AMRFinderPlus genotype row but no prespecified relevant class/subclass/symbol match.",
    "",
    rank.to_markdown(index=False),
    "",
    "## Boundary",
    "Absence of a relevant AMRFinderPlus match is not proof of a novel biological mechanism. Known mechanisms outside the AMRFinderPlus catalogue, structural variation, regulatory changes, AST error, and population structure remain possible and must be tested downstream.",
]
(OUT / "TARGET_RANKING_REPORT.md").write_text("\n".join(md) + "\n")

with open(OUT / "SHA256SUMS.txt", "w") as f:
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            h = hashlib.sha256(path.read_bytes()).hexdigest()
            f.write(f"{h}  {path.name}\n")
print((OUT / "TARGET_RANKING_REPORT.md").read_text())
