#!/usr/bin/env python3
"""Normalize sample identifiers in unitig-caller pyseer streams.

Unitig-caller versions may emit absolute paths, basenames, or assembly-version
accessions. This script maps every emitted token onto the exact assembly IDs in a
manifest and aborts on ambiguity or unknown samples. It does not alter unitig
sequences or presence/absence calls.
"""
from __future__ import annotations
import argparse
import gzip
from pathlib import Path
import pandas as pd


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True)
    p.add_argument('--manifest',required=True)
    p.add_argument('--output',required=True)
    return p.parse_args()


def opener(path,mode):
    return gzip.open(path,mode,encoding='utf-8') if str(path).endswith('.gz') else open(path,mode,encoding='utf-8')


def token_candidates(token):
    token=str(token).strip()
    raw=token
    if ':' in token:
        sample,value=token.rsplit(':',1)
    else:
        sample,value=token,'1'
    sample=Path(sample).name
    for suffix in ['.fna.gz','.fa.gz','.fasta.gz','.fna','.fa','.fasta']:
        if sample.endswith(suffix): sample=sample[:-len(suffix)]; break
    return raw,sample,value


def main():
    a=parse_args()
    m=pd.read_csv(a.manifest,dtype={'assembly_ID':str})
    ids=m.assembly_ID.dropna().astype(str).drop_duplicates().tolist()
    exact={x:x for x in ids}
    base={x.split('.')[0]:x for x in ids}
    if len(base)!=len(ids): raise SystemExit('Assembly base IDs are not one-to-one')
    seen=set(); nlines=0; ncalls=0
    with opener(a.input,'rt') as src, open(a.output,'w',encoding='utf-8') as dst:
        for line in src:
            if not line.strip(): continue
            parts=line.rstrip('\n').split()
            if not parts: continue
            feature=parts[0]
            calls=[]
            for token in parts[1:]:
                _,sample,value=token_candidates(token)
                resolved=exact.get(sample,base.get(sample.split('.')[0]))
                if resolved is None:
                    raise SystemExit(f'Unknown unitig sample ID: {sample}')
                calls.append(f'{resolved}:{value}')
                seen.add(resolved); ncalls+=1
            dst.write(feature)
            if calls: dst.write(' ' + ' '.join(calls))
            dst.write('\n'); nlines+=1
    if not nlines: raise SystemExit('No unitig lines were written')
    Path(str(a.output)+'.audit.txt').write_text(
        f'n_features={nlines}\nn_calls={ncalls}\nn_samples_seen={len(seen)}\nunknown_ids=0\n',encoding='utf-8'
    )

if __name__=='__main__': main()
