#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();p.add_argument('--root',action='append',required=True);p.add_argument('--out',required=True);p.add_argument('--boundary',required=True);p.add_argument('--include-large',action='store_true');return p.parse_args()

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    a=parse_args();excluded={'.bam','.bai','.fastq','.fq'};rows=[]
    for root_text in a.root:
        root=Path(root_text);paths=[root] if root.is_file() else sorted(root.rglob('*'))
        for p in paths:
            if not p.is_file():continue
            if not a.include_large and (p.suffix in excluded or any(str(p).endswith(x+'.gz') for x in ('.fastq','.fq'))):continue
            rows.append({'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'sha256':digest(p)})
    try:git=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    except Exception:git=None
    out={'freeze_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'git_head':git,'truth_accessed':False,'n_files':len(rows),'total_hashed_bytes':sum(r['bytes'] for r in rows),'files':rows,'boundary':a.boundary,'excluded_large_sequence_and_alignment_bytes':not a.include_large}
    Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files'},indent=2))
if __name__=='__main__':main()
