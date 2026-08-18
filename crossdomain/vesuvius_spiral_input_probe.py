#!/usr/bin/env python3
"""Probe the current PHercParis4 spiral-input annotation bundle without bulk volume download."""
from __future__ import annotations
import hashlib, json, os, re
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

BASE='https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/'
TARGETS=['umbilicus.json','same_windings.json','relative_windings.json','abs_winding.json']

def sha(b:bytes)->str:return hashlib.sha256(b).hexdigest()

def shape_summary(x,depth=0):
    if depth>=4:return {'type':type(x).__name__}
    if isinstance(x,dict):
        keys=list(x.keys())
        return {'type':'dict','length':len(x),'keys_first30':[str(k) for k in keys[:30]],'values_first5':{str(k):shape_summary(x[k],depth+1) for k in keys[:5]}}
    if isinstance(x,list):return {'type':'list','length':len(x),'items_first5':[shape_summary(v,depth+1) for v in x[:5]]}
    return {'type':type(x).__name__,'value':x if isinstance(x,(str,int,float,bool)) or x is None else repr(x)}

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/vesuvius_spiral_probe'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers['User-Agent']='Vesuvius-constraint-audit/0.1'
    root=s.get(BASE,timeout=120);root.raise_for_status();(out/'ROOT_INDEX.html').write_bytes(root.content)
    soup=BeautifulSoup(root.text,'html.parser')
    hrefs=sorted({a.get('href','') for a in soup.find_all('a') if a.get('href')})
    result={'base_url':BASE,'root_status':root.status_code,'root_sha256':sha(root.content),'hrefs':hrefs,'targets':{},'status':'RUNNING'}
    # Try root and a few common subdirectories while preserving every response status.
    prefixes=['','annotations/','constraints/','graphs/']
    for name in TARGETS:
        attempts=[];found=False
        for pref in prefixes:
            url=urljoin(BASE,pref+name)
            r=s.get(url,timeout=120)
            attempts.append({'url':url,'status':r.status_code,'bytes':len(r.content),'content_type':r.headers.get('content-type')})
            if r.status_code==200 and r.content.strip().startswith((b'{',b'[')):
                p=out/name;p.write_bytes(r.content)
                data=json.loads(r.content)
                result['targets'][name]={'found':True,'url':url,'bytes':len(r.content),'sha256':sha(r.content),'structure':shape_summary(data),'attempts':attempts}
                found=True;break
        if not found:result['targets'][name]={'found':False,'attempts':attempts}
    result['status']='CURRENT_CONSTRAINT_FILES_ACQUIRED' if all(v.get('found') for v in result['targets'].values()) else 'PARTIAL_OR_MISSING_CONSTRAINT_FILES'
    result['claim_boundary']='This is a metadata and constraint-file acquisition step. It does not fit a spiral, unwrap a scroll, detect ink, or recover text.'
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({'status':result['status'],'href_count':len(hrefs),'targets':{k:v.get('found') for k,v in result['targets'].items()}},indent=2))
if __name__=='__main__':main()
