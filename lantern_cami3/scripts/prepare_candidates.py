#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from fasta_utils import read_fasta,sanitize_seq,safe_id,seq_sha,write_record

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument('--assembly',action='append',required=True,help='SOURCE_ID=FASTA')
    p.add_argument('--source-manifest',required=True); p.add_argument('--min-length',type=int,default=1000); p.add_argument('--out',required=True); return p.parse_args()

def main():
    a=parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    sources={}
    with open(a.source_manifest,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'): sources[r['source_id']]=r
    rows=[]; seen=set(); exact={}
    with open(out/'combined_candidates.fasta','w') as fo:
        for spec in a.assembly:
            sid,path=spec.split('=',1)
            if sid not in sources: raise ValueError(f'{sid} absent from source manifest')
            n=0
            for orig,seq in read_fasta(path):
                seq=sanitize_seq(seq)
                if len(seq)<a.min_length: continue
                n+=1; cid=f'{safe_id(sid)}::{n:08d}'
                if cid in seen: raise ValueError(cid)
                seen.add(cid); sha=seq_sha(seq); exact.setdefault(sha,[]).append(cid)
                write_record(fo,cid,seq)
                rows.append({'candidate_id':cid,'source_id':sid,'original_id':orig,'length':len(seq),'n_fraction':seq.count('N')/len(seq),'sha256':sha,**{k:v for k,v in sources[sid].items() if k!='source_id'}})
    if not rows: raise SystemExit('No candidate contigs')
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields: fields.append(k)
    with open(out/'candidates.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    (out/'PREPARE_SUMMARY.json').write_text(json.dumps({'n_candidates':len(rows),'n_sources':len(set(r['source_id'] for r in rows)),'n_exact_sequence_groups':len(exact),'n_exact_duplicates':sum(max(0,len(v)-1) for v in exact.values()),'min_length':a.min_length},indent=2)+'\n')
if __name__=='__main__': main()
