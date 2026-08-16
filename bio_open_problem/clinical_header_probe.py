#!/usr/bin/env python3
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd

OUT = Path('artifact/clinical_header_probe')
OUT.mkdir(parents=True, exist_ok=True)

ref_path = Path('pg_clinical_probe/clinical_substitutions.csv')
data_zip_path = Path('pg_clinical_probe/clinical_ProteinGym_substitutions.zip')
score_zip_path = Path('pg_clinical_probe/zero_shot_clinical_substitutions_scores.zip')

ref = pd.read_csv(ref_path, low_memory=False)
report = {
    'reference_shape': list(ref.shape),
    'reference_columns': [str(c) for c in ref.columns],
    'reference_head': ref.head(8).where(pd.notna(ref), None).to_dict(orient='records'),
    'reference_dtypes': {str(c): str(t) for c, t in ref.dtypes.items()},
    'reference_unique_small': {},
}
for c in ref.columns:
    n = ref[c].nunique(dropna=True)
    if n <= 20:
        report['reference_unique_small'][str(c)] = [None if pd.isna(x) else str(x) for x in ref[c].drop_duplicates().head(30)]

for label, zp in [('data_zip', data_zip_path), ('score_zip', score_zip_path)]:
    with zipfile.ZipFile(zp) as z:
        csvs = [n for n in z.namelist() if n.lower().endswith('.csv')]
        report[label] = {
            'n_csv_files': len(csvs),
            'first_files': csvs[:10],
        }
        previews = []
        for n in csvs[:3]:
            with z.open(n) as f:
                d = pd.read_csv(f, nrows=8, low_memory=False)
            previews.append({
                'file': n,
                'columns': [str(c) for c in d.columns],
                'head': d.where(pd.notna(d), None).to_dict(orient='records'),
            })
        report[label]['previews'] = previews

# Try likely linkage keys against first score file.
with zipfile.ZipFile(score_zip_path) as z:
    n = next(x for x in z.namelist() if x.lower().endswith('.csv'))
    stem = Path(n).stem
report['first_score_stem'] = stem
for c in ref.columns:
    vals = ref[c].astype(str)
    matches = int((vals == stem).sum())
    if matches:
        report.setdefault('stem_matches', {})[str(c)] = matches

(OUT / 'clinical_header_probe.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2)[:50000])
