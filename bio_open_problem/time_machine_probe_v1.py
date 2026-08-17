#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, zipfile
from pathlib import Path
import pandas as pd


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--raw-clinical-zip',required=True);ap.add_argument('--score-zip',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    result={}
    with zipfile.ZipFile(a.raw_clinical_zip) as z:
        members=[n for n in z.namelist() if n.lower().endswith(('.csv','.tsv','.txt'))]
        result['raw_member_count']=len(members);result['raw_members']=members[:100]
        previews=[]
        for n in members[:20]:
            try:
                with z.open(n) as f:
                    df=pd.read_csv(f,sep=None,engine='python',nrows=5,low_memory=False)
                previews.append({'member':n,'columns':list(df.columns),'rows':df.astype(str).head(3).to_dict(orient='records')})
            except Exception as e: previews.append({'member':n,'error':repr(e)})
        result['raw_previews']=previews
    with zipfile.ZipFile(a.score_zip) as z:
        members=[n for n in z.namelist() if n.lower().endswith('.csv')]
        result['score_member_count']=len(members)
        n=members[0]
        with z.open(n) as f:df=pd.read_csv(f,nrows=5,low_memory=False)
        result['score_example']={'member':n,'columns':list(df.columns),'rows':df.astype(str).head(3).to_dict(orient='records')}
    (out/'probe.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'raw_member_count':result['raw_member_count'],'raw_members':result['raw_members'][:20],'first_columns':result['raw_previews'][0].get('columns') if result['raw_previews'] else None},indent=2))
if __name__=='__main__':main()
