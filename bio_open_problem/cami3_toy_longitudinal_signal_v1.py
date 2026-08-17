#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, tarfile
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd


def parse_profile(path: Path):
    sample_id=None
    rows=[]
    with path.open('r',encoding='utf-8',errors='replace') as f:
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('@') and 'SampleID:' in line:
                import re
                m=re.search(r'@SampleID:([^\s]+)',line)
                if m: sample_id=m.group(1)
            if not line or line.startswith('#') or line.startswith('@'): continue
            parts=line.split('\t')
            if len(parts)<5: continue
            taxid,rank,taxpath,taxpathsn,pct=parts[:5]
            try:p=float(pct)
            except:continue
            rows.append((taxid,rank,p))
    return sample_id,rows

def sample_index(name:str):
    import re
    m=re.search(r'(?:sample[_-]?)(\d+)',name,re.I)
    if not m: m=re.search(r'(^|[^0-9])(\d+)(?:[^0-9]|$)',name)
    return int(m.group(1 if m and m.lastindex==1 else 2)) if m else None

def bray(a,b):
    keys=set(a)|set(b)
    sa=sum(a.values());sb=sum(b.values())
    if sa+sb==0:return 0.0
    return 1.0-sum(abs(a.get(k,0)-b.get(k,0)) for k in keys)/(sa+sb)

def js_similarity(a,b):
    keys=sorted(set(a)|set(b));
    if not keys:return 0.0
    x=np.array([a.get(k,0.0) for k in keys],float);y=np.array([b.get(k,0.0) for k in keys],float)
    if x.sum()==0 or y.sum()==0:return 0.0
    x/=x.sum();y/=y.sum();m=.5*(x+y)
    def kl(p,q):
        z=p>0
        return float(np.sum(p[z]*np.log2(p[z]/q[z])))
    js=.5*kl(x,m)+.5*kl(y,m)
    return 1.0-js

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--profiles-tar',required=True);ap.add_argument('--out-dir',required=True);a=ap.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    extracted=out/'profiles';extracted.mkdir(exist_ok=True)
    with tarfile.open(a.profiles_tar,'r:*') as tf: tf.extractall(extracted)
    files=[p for p in extracted.rglob('*') if p.is_file()]
    profiles={};meta={}
    for p in files:
        sid,rows=parse_profile(p)
        if not rows:continue
        idx=sample_index(sid or p.name)
        if idx is None or not (0<=idx<=19): continue
        species={tax:pct for tax,rank,pct in rows if rank=='species'}
        genus={tax:pct for tax,rank,pct in rows if rank=='genus'}
        chosen=species if species else genus
        profiles[idx]=chosen;meta[idx]={'sample_id':sid,'file':str(p.relative_to(extracted)),'rank':'species' if species else 'genus','n_taxa':len(chosen)}
    if len(profiles)<20: raise RuntimeError(f'Expected 20 profiles, parsed {len(profiles)}: {meta}')
    pair_rows=[]
    for i in range(20):
        for j in range(i+1,20):
            same=(i//2)==(j//2)
            pair_rows.append({'i':i,'j':j,'same_individual':same,'bray_similarity':bray(profiles[i],profiles[j]),'js_similarity':js_similarity(profiles[i],profiles[j])})
    df=pd.DataFrame(pair_rows);df.to_csv(out/'pairwise_similarity.csv',index=False)
    rng=np.random.default_rng(20260817)
    observed={}
    for metric in ['bray_similarity','js_similarity']:
        same=df[df.same_individual][metric].to_numpy();other=df[~df.same_individual][metric].to_numpy()
        delta=float(same.mean()-other.mean())
        vals=df[metric].to_numpy();labels=df.same_individual.to_numpy();null=[]
        for _ in range(10000):
            sh=rng.permutation(labels);null.append(float(vals[sh].mean()-vals[~sh].mean()))
        p=(1+sum(x>=delta for x in null))/(1+len(null))
        observed[metric]={'same_individual_mean':float(same.mean()),'between_individual_mean':float(other.mean()),'delta':delta,'one_sided_permutation_p':float(p)}
    # Simple nearest-neighbor identification: for each sample, does the most similar other sample equal its paired longitudinal sample?
    nn={}
    for metric in ['bray_similarity','js_similarity']:
        wins=0;rows=[]
        for i in range(20):
            candidates=[]
            for j in range(20):
                if i==j:continue
                r=df[((df.i==min(i,j))&(df.j==max(i,j)))].iloc[0]
                candidates.append((float(r[metric]),j))
            candidates.sort(reverse=True);best=candidates[0][1];truth=i^1
            wins+=int(best==truth);rows.append({'sample':i,'nearest':best,'paired_truth':truth,'correct':best==truth,'similarity':candidates[0][0]})
        nn[metric]={'correct':wins,'total':20,'accuracy':wins/20,'rows':rows}
    result={'status':'TOY_LONGITUDINAL_SIGNAL_QUANTIFIED','n_samples':20,'n_individuals':10,'n_within_pairs':int(df.same_individual.sum()),'n_between_pairs':int((~df.same_individual).sum()),'metrics':observed,'nearest_neighbor_pair_identification':nn,
            'claim_boundary':'This uses gold-standard toy taxonomic profiles only. It tests whether longitudinal identity contains signal that could justify a prespecified prior. It is not a metagenomic profiler, not a blind CAMI III challenge result, and reference-based performance on this public-genome toy set can be inflated.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
