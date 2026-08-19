#!/usr/bin/env python3
"""Match sequence-audited residual resistant K. pneumoniae to susceptible controls.

Matching uses only pre-genomic metadata and Kleborate's prespecified known-mechanism flag.
It does not inspect candidate variants, unitigs, or GWAS results. A global one-to-one Hungarian
assignment minimizes a fixed cost over ST, BioProject, country, collection year, and isolation
source. The output is a discovery cohort input, not evidence of a resistance mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def args() -> argparse.Namespace:
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',required=True,help='QC-passed matrix manifest with BioProject metadata')
    p.add_argument('--kleborate-labels',required=True,help='labels_with_kleborate.csv from independent audit')
    p.add_argument('--out',required=True)
    return p.parse_args()


def text(v) -> str:
    return '' if pd.isna(v) else str(v)


def main() -> None:
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    m=pd.read_csv(a.manifest,dtype={'assembly_ID':str,'BioProject':str})
    k=pd.read_csv(a.kleborate_labels,dtype={'assembly_ID':str,'Kleborate_ST':str})
    keep=['assembly_ID','has_kleborate_known_colistin_evidence','has_mcr_family_call','has_kleborate_col_mutation','Kleborate_ST','Col_mutations','kleborate_result_present']
    keep=[c for c in keep if c in k.columns]
    x=m.merge(k[keep].drop_duplicates('assembly_ID'),on='assembly_ID',how='left',validate='one_to_one')
    known=x.get('has_kleborate_known_colistin_evidence',False)
    if not isinstance(known,pd.Series): known=pd.Series(False,index=x.index)
    x['has_sequence_known_colistin_evidence']=known.fillna(False).astype(bool)
    r=x[x.phenotype.astype(str).eq('R') & ~x.has_sequence_known_colistin_evidence].copy().sort_values('assembly_ID').reset_index(drop=True)
    s=x[x.phenotype.astype(str).eq('S') & ~x.has_sequence_known_colistin_evidence].copy().sort_values('assembly_ID').reset_index(drop=True)
    if len(s)<len(r): raise RuntimeError(f'Insufficient susceptible controls: R={len(r)} S={len(s)}')

    cost=np.zeros((len(r),len(s)),dtype=float)
    for i,rr in r.iterrows():
        for j,ss in s.iterrows():
            c=0.0
            rst,sst=text(rr.get('Kleborate_ST')),text(ss.get('Kleborate_ST'))
            if rst and sst: c += 12.0*(rst!=sst)
            elif rst or sst: c += 8.0
            c += 6.0*(text(rr.get('BioProject'))!=text(ss.get('BioProject')))
            c += 4.0*(text(rr.get('ISO_country_code'))!=text(ss.get('ISO_country_code')))
            c += 1.0*(text(rr.get('isolation_source_category'))!=text(ss.get('isolation_source_category')))
            yr,ys=rr.get('collection_year'),ss.get('collection_year')
            c += min(abs(float(yr)-float(ys)),10.0)/5.0 if pd.notna(yr) and pd.notna(ys) else 1.0
            cost[i,j]=c+j*1e-9
    ri,sj=linear_sum_assignment(cost)
    if len(ri)!=len(r): raise RuntimeError('Assignment did not match every residual resistant isolate')
    selected=s.iloc[sj].copy().reset_index(drop=True)
    pairs=[]
    for n,(ii,jj) in enumerate(zip(ri,sj),1):
        rr=r.iloc[ii]; ss=s.iloc[jj]
        pairs.append({'match_id':f'M{n:04d}','resistant_assembly':rr.assembly_ID,'susceptible_assembly':ss.assembly_ID,'cost':float(cost[ii,jj]),'same_ST':text(rr.get('Kleborate_ST'))!='' and text(rr.get('Kleborate_ST'))==text(ss.get('Kleborate_ST')),'same_BioProject':text(rr.get('BioProject'))==text(ss.get('BioProject')),'same_country':text(rr.get('ISO_country_code'))!='' and text(rr.get('ISO_country_code'))==text(ss.get('ISO_country_code')),'same_isolation_category':text(rr.get('isolation_source_category'))!='' and text(rr.get('isolation_source_category'))==text(ss.get('isolation_source_category'))})
    pair=pd.DataFrame(pairs)
    r=r.copy(); selected=selected.copy(); r['match_id']=pair.match_id.values; selected['match_id']=pair.match_id.values
    cohort=pd.concat([r,selected],ignore_index=True).sort_values(['match_id','phenotype'])
    if cohort.assembly_ID.duplicated().any(): raise RuntimeError('Duplicate assembly after matching')
    if cohort.phenotype.value_counts().to_dict()!={'R':len(r),'S':len(r)}: raise RuntimeError('Matched cohort is not balanced')
    cohort.to_csv(out/'DEEP_RESIDUAL_MATCHED_COHORT.csv',index=False); pair.to_csv(out/'MATCH_PAIRS.csv',index=False)
    summary={'n_qc_manifest':int(len(m)),'n_residual_resistant':int(len(r)),'n_eligible_susceptible_pool':int(len(s)),'n_matched_per_class':int(len(r)),'n_matched_total':int(len(cohort)),'mean_match_cost':float(pair.cost.mean()),'median_match_cost':float(pair.cost.median()),'max_match_cost':float(pair.cost.max()),'same_ST_pairs':int(pair.same_ST.sum()),'same_BioProject_pairs':int(pair.same_BioProject.sum()),'same_country_pairs':int(pair.same_country.sum()),'same_isolation_category_pairs':int(pair.same_isolation_category.sum()),'weights':{'ST_mismatch':12,'one_ST_missing':8,'BioProject_mismatch':6,'country_mismatch':4,'isolation_category_mismatch':1,'year_difference':'min(abs(year difference),10)/5'},'selection_boundary':'Matching used no candidate sequence feature or association result. Resistant isolates had no mcr-family call and no Kleborate-detected MgrB/PmrB inactivation; this does not exclude all known chromosomal mechanisms.'}
    (out/'MATCH_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']; (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
