#!/usr/bin/env python3
"""Phenotype-blind genetic clustering followed by discovery/validation CMH tests.

This prespecified sensitivity analysis asks whether a feature is associated with
colistin resistance repeatedly *within* genetic clusters rather than merely marking
one resistant lineage. Clusters are derived from Mash distances without phenotype
information. Candidate selection uses discovery samples only.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.cluster import KMeans


def args():
    p=argparse.ArgumentParser(); p.add_argument('--manifest',required=True); p.add_argument('--rtab',required=True); p.add_argument('--distance',required=True); p.add_argument('--feature-meta',required=True); p.add_argument('--out',required=True); p.add_argument('--clusters',default='8,12,16'); p.add_argument('--alpha',type=float,default=.05); p.add_argument('--min-af',type=float,default=.01); p.add_argument('--max-af',type=float,default=.99); p.add_argument('--min-stratum-class',type=int,default=3); return p.parse_args()

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p); v=p[ok]
    if len(v):
        o=np.argsort(v); r=v[o]; q=np.minimum.accumulate((r*len(r)/np.arange(1,len(r)+1))[::-1])[::-1]; out[np.flatnonzero(ok)[o]]=np.clip(q,0,1)
    return out

def pcoa(d,n=40):
    x=d.to_numpy(float); x=(x+x.T)/2; np.fill_diagonal(x,0); m=len(x); j=np.eye(m)-np.ones((m,m))/m; b=-.5*j@(x*x)@j
    e,v=np.linalg.eigh(b); o=np.argsort(e)[::-1]; e=e[o]; v=v[:,o]; keep=e>max(e[0]*1e-12,1e-15); e=e[keep]; v=v[:,keep]; k=min(n,len(e)); z=v[:,:k]*np.sqrt(e[:k]); sd=z.std(0); sd[sd==0]=1; return (z-z.mean(0))/sd

def cmh_for_feature(x,y,cluster,min_class):
    numer=0.; varsum=0.; or_num=0.; or_den=0.; informative=0; positive=0; details=[]
    for g in sorted(set(cluster)):
        idx=cluster==g; yy=y[idx]; xx=x[idx]
        r=int(yy.sum()); s=int(len(yy)-r)
        if r<min_class or s<min_class: continue
        a=int(((xx==1)&(yy==1)).sum()); b=int(((xx==1)&(yy==0)).sum()); c=int(((xx==0)&(yy==1)).sum()); d=int(((xx==0)&(yy==0)).sum()); n=a+b+c+d
        if n<=1 or (a+b)==0 or (c+d)==0: continue
        exp=(a+b)*(a+c)/n; var=(a+b)*(c+d)*(a+c)*(b+d)/(n*n*(n-1))
        if var<=0: continue
        numer+=a-exp; varsum+=var; or_num+=a*d/n; or_den+=b*c/n; informative+=1
        aa,bb,cc,dd=map(float,(a,b,c,d));
        if min(aa,bb,cc,dd)==0: aa+=.5;bb+=.5;cc+=.5;dd+=.5
        lor=math.log((aa*dd)/(bb*cc)); positive+=int(lor>0); details.append({'cluster':int(g),'a':a,'b':b,'c':c,'d':d,'log_or':lor})
    if informative==0 or varsum<=0: return {'p':np.nan,'z':np.nan,'or':np.nan,'n_strata':0,'positive_strata':0,'details':details}
    z=numer/math.sqrt(varsum); p=float(norm.sf(z)); common=(or_num/or_den) if or_den>0 else float('inf')
    return {'p':p,'z':z,'or':common,'n_strata':informative,'positive_strata':positive,'details':details}

def run_matrix(x,y,clusters,features,min_class):
    rows=[]; detail={}
    for j,f in enumerate(features):
        r=cmh_for_feature(x[:,j],y,clusters,min_class); detail[str(f)]=r.pop('details'); r['feature']=str(f); rows.append(r)
    z=pd.DataFrame(rows); z['q']=bh(z.p); return z,detail

def main():
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); ks=[int(x) for x in a.clusters.split(',')]
    m=pd.read_csv(a.manifest,dtype={'assembly_ID':str}).drop_duplicates('assembly_ID'); m['y']=m.phenotype.astype(str).eq('R').astype(int); ids=m.assembly_ID.tolist()
    r=pd.read_csv(a.rtab,sep='\t',index_col=0); r.index=r.index.astype(str); r.columns=r.columns.astype(str); r=r.apply(pd.to_numeric,errors='coerce').fillna(0).astype(np.uint8)
    d=pd.read_csv(a.distance,sep='\t',index_col=0); d.index=d.index.astype(str); d.columns=d.columns.astype(str)
    if set(ids)!=set(r.columns) or set(ids)!=set(d.index) or set(d.index)!=set(d.columns): raise SystemExit('Sample ID mismatch')
    r=r.loc[:,ids]; d=d.loc[ids,ids]; coords=pcoa(d,40); x=r.T.to_numpy(np.uint8); y=m.y.to_numpy(np.uint8); disc=m.split.eq('discovery').to_numpy(); val=m.split.eq('validation').to_numpy()
    af=x[disc].mean(0); eligible=(af>=a.min_af)&(af<=a.max_af); feats=r.index.to_numpy()[eligible]; xd=x[disc][:,eligible]; xv=x[val][:,eligible]; yd=y[disc]; yv=y[val]
    model_disc=[]; model_val=[]; all_details={}
    for k in ks:
        labels=KMeans(n_clusters=k,random_state=20260819,n_init=50).fit_predict(coords)
        m[f'genetic_cluster_k{k}']=labels
        zd,dd=run_matrix(xd,yd,labels[disc],feats,a.min_stratum_class); zv,dv=run_matrix(xv,yv,labels[val],feats,a.min_stratum_class)
        zd=zd.add_suffix(f'_k{k}').rename(columns={f'feature_k{k}':'feature'}); zv=zv.add_suffix(f'_k{k}').rename(columns={f'feature_k{k}':'feature'})
        zd.to_csv(out/f'DISCOVERY_CMH_K{k}.csv',index=False); zv.to_csv(out/f'VALIDATION_CMH_K{k}.csv',index=False)
        model_disc.append(zd); model_val.append(zv); all_details[f'k{k}']={'discovery':dd,'validation':dv}
    merged=model_disc[0]
    for z in model_disc[1:]: merged=merged.merge(z,on='feature',how='outer')
    pcols=[f'p_k{k}' for k in ks]; qcols=[f'q_k{k}' for k in ks]; ocols=[f'or_k{k}' for k in ks]; ncols=[f'n_strata_k{k}' for k in ks]
    merged['p_max_discovery']=merged[pcols].max(1); merged['q_max_discovery']=merged[qcols].max(1); merged['or_min_discovery']=merged[ocols].min(1); merged['min_informative_strata_discovery']=merged[ncols].min(1)
    merged['discovery_stable']=(merged.q_max_discovery<=a.alpha)&(merged.or_min_discovery>1)&(merged.min_informative_strata_discovery>=2)
    frozen=merged[merged.discovery_stable].copy().sort_values(['q_max_discovery','p_max_discovery']); frozen.to_csv(out/'FROZEN_LINEAGE_STRATIFIED_CANDIDATES.csv',index=False)
    vmerge=model_val[0]
    for z in model_val[1:]: vmerge=vmerge.merge(z,on='feature',how='outer')
    evidence=frozen.merge(vmerge,on='feature',how='left',suffixes=('','_validation'))
    vpcols=[f'p_k{k}' for k in ks]; vocols=[f'or_k{k}' for k in ks]; vncols=[f'n_strata_k{k}' for k in ks]
    if len(evidence):
        evidence['validation_p_max']=evidence[vpcols].max(1); evidence['validation_q_across_candidates']=bh(evidence.validation_p_max); evidence['validation_or_min']=evidence[vocols].min(1); evidence['validation_min_informative_strata']=evidence[vncols].min(1)
        evidence['heldout_lineage_replication']=(evidence.validation_q_across_candidates<=a.alpha)&(evidence.validation_or_min>1)&(evidence.validation_min_informative_strata>=2)
    else:
        evidence['validation_p_max']=pd.Series(dtype=float); evidence['validation_q_across_candidates']=pd.Series(dtype=float); evidence['validation_or_min']=pd.Series(dtype=float); evidence['validation_min_informative_strata']=pd.Series(dtype=float); evidence['heldout_lineage_replication']=pd.Series(dtype=bool)
    evidence.to_csv(out/'LINEAGE_STRATIFIED_DISCOVERY_VALIDATION_EVIDENCE.csv',index=False); strict=evidence[evidence.heldout_lineage_replication.eq(True)].copy(); strict.to_csv(out/'STRICT_LINEAGE_STRATIFIED_REPLICATES.csv',index=False)
    m.to_csv(out/'PHENOTYPE_BLIND_GENETIC_CLUSTERS.csv',index=False); (out/'CMH_DETAILS.json').write_text(json.dumps(all_details,ensure_ascii=False,default=str)+'\n')
    summary={'n_all':len(m),'n_discovery':int(disc.sum()),'n_validation':int(val.sum()),'n_eligible_features':len(feats),'cluster_counts':ks,'n_discovery_stable':len(frozen),'n_heldout_lineage_replicated':len(strict),'strict_features':strict.feature.astype(str).tolist() if len(strict) else [],'status':'LINEAGE_STRATIFIED_CANDIDATES_REQUIRE_KNOWN_MECHANISM_AND_UNITIG_CONCORDANCE' if len(strict) else 'NO_FEATURE_SURVIVED_LINEAGE_STRATIFIED_DISCOVERY_VALIDATION','boundary':'Within-lineage statistical replication does not establish novelty or causality. Known mechanisms, sequence context, independent tools, literature and biological validation remain required.'}
    (out/'LINEAGE_STRATIFIED_CMH_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    lines=['# Lineage-stratified CMH sensitivity analysis','',f"- Eligible features: **{len(feats):,}**",f"- Discovery-stable features: **{len(frozen):,}**",f"- Untouched validation replicates: **{len(strict):,}**",'',summary['boundary']]
    if len(evidence): lines += ['','## Evidence-ranked candidates','',evidence.head(50).to_markdown(index=False)]
    (out/'LINEAGE_STRATIFIED_CMH_REPORT.md').write_text('\n'.join(lines)+'\n')
    hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']; (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n'); print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__': main()
