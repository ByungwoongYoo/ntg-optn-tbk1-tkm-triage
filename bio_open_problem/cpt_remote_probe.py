#!/usr/bin/env python3
from __future__ import annotations
import json, os, re
from pathlib import Path
import requests
from remotezip import RemoteZip

OUT=Path(os.environ.get('OUT_DIR','artifact/cpt_probe'));OUT.mkdir(parents=True,exist_ok=True)
record=requests.get('https://zenodo.org/api/records/7954657',timeout=120).json()
targets=['ACADM','BRCA1','TP53','PTEN','PSEN1','APOE']
results=[]
for f in record.get('files',[]):
    name=f['key']
    if not name.endswith('.zip'):continue
    url=f.get('links',{}).get('content') or f.get('links',{}).get('self')
    row={'archive':name,'bytes':f.get('size'),'url':url}
    try:
        rz=RemoteZip(url)
        names=rz.namelist()
        row['n_members']=len(names);row['first_members']=names[:30]
        row['target_matches']={t:[n for n in names if re.search(rf'(^|[/_.-]){re.escape(t)}([/_.-]|$)',n,re.I)][:10] for t in targets}
        rz.close()
    except Exception as e:row['error']=repr(e)
    results.append(row)
(OUT/'probe.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
