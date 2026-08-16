#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, zipfile
from collections import Counter
from pathlib import Path
import pandas as pd

ALLOWED=["PoET","TranceptEVE_L","GEMME","EVE","ESM1b"]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--clinical-score-zip',required=True);ap.add_argument('--out-dir',required=True)
    args=ap.parse_args();out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    labels=Counter();total=0;complete=0;files=0;examples=[]
    with zipfile.ZipFile(args.clinical_score_zip) as z:
        members=[n for n in z.namelist() if n.lower().endswith('.csv')]
        for i,n in enumerate(members,1):
            with z.open(n) as f:d=pd.read_csv(f,low_memory=False)
            files+=1;total+=len(d)
            if 'DMS_bin_score' in d.columns:
                c=d['DMS_bin_score'].astype(str).str.strip().value_counts(dropna=False)
                labels.update(c.to_dict())
            if set(ALLOWED).issubset(d.columns):complete+=int(d[ALLOWED].apply(pd.to_numeric,errors='coerce').notna().all(axis=1).sum())
            if i<=3: examples.append({'file':Path(n).name,'columns':list(d.columns),'label_examples':d.get('DMS_bin_score',pd.Series(dtype=str)).astype(str).head(20).tolist()})
    result={'files':files,'rows':total,'label_counts':dict(labels),'complete_case_rows_five_models':complete,'examples':examples,
            'vus_rows_present':any('uncertain' in str(k).lower() or 'vus' in str(k).lower() for k in labels)}
    (out/'vus_probe.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
