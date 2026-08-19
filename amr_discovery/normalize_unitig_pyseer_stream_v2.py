#!/usr/bin/env python3
"""Normalize unitig-caller sample IDs while preserving pyseer's `feature | calls` grammar."""
from __future__ import annotations
import argparse,gzip
from pathlib import Path
import pandas as pd


def args():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--manifest',required=True); p.add_argument('--output',required=True); return p.parse_args()

def op(path,mode): return gzip.open(path,mode,encoding='utf-8') if str(path).endswith('.gz') else open(path,mode,encoding='utf-8')

def clean_sample(token):
    token=token.strip()
    if ':' in token:
        sample,value=token.rsplit(':',1)
    else:
        sample,value=token,'1'
    sample=Path(sample).name
    for suffix in ('.fna.gz','.fasta.gz','.fa.gz','.fna','.fasta','.fa'):
        if sample.endswith(suffix): sample=sample[:-len(suffix)]; break
    return sample,value

def main():
    a=args(); m=pd.read_csv(a.manifest,dtype={'assembly_ID':str}); ids=m.assembly_ID.dropna().astype(str).drop_duplicates().tolist()
    exact={x:x for x in ids}; base={x.split('.')[0]:x for x in ids}
    if len(base)!=len(ids): raise SystemExit('Nonunique accession base IDs')
    features=calls=0; seen=set()
    with op(a.input,'rt') as src, open(a.output,'w',encoding='utf-8') as dst:
        for raw in src:
            line=raw.strip()
            if not line: continue
            if '|' in line:
                left,right=line.split('|',1); feature=left.strip(); tokens=right.split(); delimiter=True
            else:
                parts=line.split(); feature=parts[0]; tokens=parts[1:]; delimiter=False
            normalized=[]
            for token in tokens:
                if token=='|': continue
                sample,value=clean_sample(token)
                resolved=exact.get(sample,base.get(sample.split('.')[0]))
                if resolved is None: raise SystemExit(f'Unknown unitig sample ID: {sample}')
                normalized.append(f'{resolved}:{value}'); seen.add(resolved); calls+=1
            if delimiter:
                dst.write(feature+' |')
                if normalized: dst.write(' '+' '.join(normalized))
                dst.write('\n')
            else:
                dst.write(feature)
                if normalized: dst.write(' '+' '.join(normalized))
                dst.write('\n')
            features+=1
    if features==0: raise SystemExit('No unitig features normalized')
    Path(a.output+'.audit.txt').write_text(f'n_features={features}\nn_calls={calls}\nn_samples_seen={len(seen)}\nunknown_ids=0\ndelimiter_preserved=true\n')
if __name__=='__main__': main()
