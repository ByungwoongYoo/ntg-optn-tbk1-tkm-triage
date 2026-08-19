#!/usr/bin/env python3
"""Audit how AMRFinderPlus genotype class/subclass fields map to selected AST drugs."""
from __future__ import annotations
import pathlib, urllib.request, hashlib
import duckdb, pandas as pd

RELEASE='2026-07'
BASE=f'https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/{RELEASE}'
OUT=pathlib.Path('amr_mapping_audit'); OUT.mkdir(exist_ok=True)
for name in ['phenotype.parquet','genotype.parquet']:
    path=OUT/name
    if not path.exists():
        req=urllib.request.Request(f'{BASE}/{name}',headers={'User-Agent':'AMR-mapping-audit/1.0'})
        with urllib.request.urlopen(req,timeout=120) as r, open(path,'wb') as f:
            while b:=r.read(1024*1024): f.write(b)

p=str((OUT/'phenotype.parquet').resolve()); g=str((OUT/'genotype.parquet').resolve())
con=duckdb.connect(); con.execute('PRAGMA threads=4')

pairs=[
('Escherichia coli','ampicillin'),
('Neisseria gonorrhoeae','ciprofloxacin'),
('Klebsiella pneumoniae','cefoxitin'),
('Klebsiella pneumoniae','ertapenem'),
('Mycobacterium tuberculosis','rifampin'),
('Klebsiella pneumoniae','imipenem'),
('Klebsiella pneumoniae','piperacillin-tazobactam'),
('Streptococcus pneumoniae','trimethoprim-sulfamethoxazole'),
('Acinetobacter baumannii','imipenem'),
('Klebsiella pneumoniae','meropenem'),
('Klebsiella pneumoniae','tobramycin'),
]
rows=[]
for org,drug in pairs:
    df=con.execute(f"""
    WITH ids AS (
      SELECT DISTINCT BioSample_ID, assembly_ID
      FROM read_parquet('{p}')
      WHERE organism=? AND lower(antibiotic_name)=?
        AND lower(trim(resistance_phenotype)) IN ('resistant','r','non-susceptible','nonsusceptible','ns')
        AND assembly_ID IS NOT NULL
    )
    SELECT ? AS organism, ? AS phenotype_drug,
           coalesce(class,'') AS class, coalesce(subclass,'') AS subclass,
           coalesce(split_subclass,'') AS split_subclass,
           coalesce(antibiotic_name,'') AS genotype_antibiotic_name,
           coalesce(antibiotic_ontology,'') AS genotype_antibiotic_ontology,
           coalesce(amr_element_symbol,'') AS amr_element_symbol,
           coalesce(amrfinderplus_method,'') AS amrfinderplus_method,
           count(DISTINCT g.BioSample_ID) AS resistant_biosamples
    FROM read_parquet('{g}') g JOIN ids USING(BioSample_ID,assembly_ID)
    GROUP BY ALL ORDER BY resistant_biosamples DESC
    LIMIT 100
    """,[org,drug,org,drug]).df()
    rows.append(df)
out=pd.concat(rows,ignore_index=True)
out.to_csv(OUT/'top_genotype_rows_for_resistant_pairs.csv',index=False)

# Global vocabularies, with counts, to support a transparent map rather than exact-ontology matching.
for col in ['class','subclass','split_subclass','antibiotic_name','amrfinderplus_method']:
    safe=col.replace('/','_')
    con.execute(f"""
      COPY (SELECT coalesce({col},'') AS value, count(*) n, count(DISTINCT BioSample_ID) biosamples
            FROM read_parquet('{g}') GROUP BY 1 ORDER BY biosamples DESC)
      TO '{str((OUT/(safe+'_vocabulary.csv')).resolve())}' (HEADER, DELIMITER ',')
    """)

# Per-isolate count of any genotype rows for the selected pairs; useful to separate no-annotation from no-relevant-annotation.
labels=[]
for org,drug in pairs:
    df=con.execute(f"""
    WITH ph AS (
      SELECT DISTINCT BioSample_ID, assembly_ID,
             CASE WHEN lower(trim(resistance_phenotype)) IN ('resistant','r','non-susceptible','nonsusceptible','ns') THEN 'R'
                  WHEN lower(trim(resistance_phenotype)) IN ('susceptible','sensitive','s') THEN 'S' END y
      FROM read_parquet('{p}') WHERE organism=? AND lower(antibiotic_name)=? AND assembly_ID IS NOT NULL
    ), gc AS (
      SELECT BioSample_ID, assembly_ID, count(*) genotype_rows,
             string_agg(DISTINCT class, ';') AS classes,
             string_agg(DISTINCT subclass, ';') AS subclasses,
             string_agg(DISTINCT amr_element_symbol, ';') AS elements
      FROM read_parquet('{g}') GROUP BY 1,2
    )
    SELECT ? AS organism, ? AS antibiotic_name, ph.*, coalesce(genotype_rows,0) genotype_rows,
           coalesce(classes,'') classes, coalesce(subclasses,'') subclasses, coalesce(elements,'') elements
    FROM ph LEFT JOIN gc USING(BioSample_ID,assembly_ID) WHERE y IS NOT NULL
    """,[org,drug,org,drug]).df()
    labels.append(df)
pd.concat(labels,ignore_index=True).to_csv(OUT/'selected_pair_isolate_genotype_context.csv',index=False)

with open(OUT/'SHA256SUMS.txt','w') as f:
    for x in sorted(OUT.iterdir()):
        if x.is_file() and x.name!='SHA256SUMS.txt':
            h=hashlib.sha256(x.read_bytes()).hexdigest(); f.write(f'{h}  {x.name}\n')
print(out.groupby(['organism','phenotype_drug']).head(20).to_string(index=False))
