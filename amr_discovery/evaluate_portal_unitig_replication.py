#!/usr/bin/env python3
"""Evaluate discovery-selected unitigs in an untouched source-held-out cohort.

This is a statistical replication gate only. It does not establish causality, biological
novelty, clinical validity, or resistance mechanism status. Candidates must subsequently be
intersected with the independent sequence-level known-mechanism audit and context review.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm


def args():
    p=argparse.ArgumentParser()
    p.add_argument('--selection',required=True); p.add_argument('--manifest',required=True)
    p.add_argument('--validation-rtab',required=True); p.add_argument('--all-rtab',required=True)
    p.add_argument('--whole-pyseer',nargs='+',required=True); p.add_argument('--out',required=True)
    p.add_argument('--alpha',type=float,default=0.05); p.add_argument('--min-validation-present',type=int,default=5)
    return p.parse_args()


def bh(s):
    p=pd.to_numeric(s,errors='coerce').to_numpy(float); out=np.full(len(p),np.nan); ok=np.isfinite(p); v=p[ok]
    if len(v):
        o=np.argsort(v); r=v[o]; q=np.minimum.accumulate((r*len(r)/np.arange(1,len(r)+1))[::-1])[::-1]
        out[np.flatnonzero(ok)[o]]=np.clip(q,0,1)
    return pd.Series(out,index=s.index)


def rtab(path):
    x=pd.read_csv(path,sep='\t',index_col=0); x.index=x.index.astype(str).str.upper(); x.columns=x.columns.astype(str)
    return x.apply(pd.to_numeric,errors='coerce').fillna(0).astype(int)


def cont(row,meta):
    z=meta[['assembly_ID','phenotype']].copy(); z['x']=z.assembly_ID.map(row).fillna(0).astype(int); r=z.phenotype.astype(str).eq('R')
    return int((r&(z.x==1)).sum()),int((~r&(z.x==1)).sum()),int((r&(z.x==0)).sum()),int((~r&(z.x==0)).sum())


def orci(a,b,c,d):
    aa,bb,cc,dd=map(float,(a,b,c,d))
    if min(aa,bb,cc,dd)==0: aa+=.5;bb+=.5;cc+=.5;dd+=.5
    l=math.log((aa*dd)/(bb*cc)); se=math.sqrt(1/aa+1/bb+1/cc+1/dd)
    return math.exp(l),math.exp(l-1.96*se),math.exp(l+1.96*se)


def random_meta(row,meta,col):
    eff=[]
    for g,sub in meta.groupby(col,dropna=False):
        a,b,c,d=cont(row,sub)
        if (a+c)==0 or (b+d)==0 or (a+b)==0 or (c+d)==0: continue
        aa,bb,cc,dd=map(float,(a,b,c,d))
        if min(aa,bb,cc,dd)==0: aa+=.5;bb+=.5;cc+=.5;dd+=.5
        l=math.log((aa*dd)/(bb*cc)); var=1/aa+1/bb+1/cc+1/dd; eff.append((str(g),l,var,a,b,c,d))
    if not eff: return {'n':0}
    y=np.array([e[1] for e in eff]); w=1/np.array([e[2] for e in eff]); mu=float(np.sum(w*y)/np.sum(w))
    q=float(np.sum(w*(y-mu)**2)); df=len(y)-1; cval=float(np.sum(w)-np.sum(w*w)/np.sum(w)); tau=max(0,(q-df)/cval) if df>0 and cval>0 else 0
    wr=1/(np.array([e[2] for e in eff])+tau); mur=float(np.sum(wr*y)/np.sum(wr)); se=math.sqrt(1/float(np.sum(wr)))
    return {'n':len(eff),'or':math.exp(mur),'lo':math.exp(mur-1.96*se),'hi':math.exp(mur+1.96*se),'p':float(norm.sf(mur/se)) if se else 1.0,'I2':max(0,(q-df)/q*100) if q>0 and df>0 else 0,'positive':int(sum(e[1]>0 for e in eff)),'details':[{'group':e[0],'log_or':e[1],'var':e[2],'a':e[3],'b':e[4],'c':e[5],'d':e[6]} for e in eff]}


def pyseer(path,i):
    x=pd.read_csv(path,sep='\t',dtype={'variant':str}); x.variant=x.variant.astype(str).str.upper()
    x['beta']=pd.to_numeric(x['beta'],errors='coerce'); x['lrt-pvalue']=pd.to_numeric(x['lrt-pvalue'],errors='coerce'); x['q']=bh(x['lrt-pvalue'])
    return x[['variant','beta','lrt-pvalue','q']].drop_duplicates('variant').rename(columns={'beta':f'b{i}','lrt-pvalue':f'p{i}','q':f'q{i}'})


def main():
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    sel=pd.read_csv(a.selection,dtype={'canonical_sequence':str,'candidate_id':str}); sel.canonical_sequence=sel.canonical_sequence.str.upper()
    meta=pd.read_csv(a.manifest,dtype={'assembly_ID':str}); val=meta[meta.split.eq('validation')].copy(); vr=rtab(a.validation_rtab); ar=rtab(a.all_rtab)
    whole=None
    for i,p in enumerate(a.whole_pyseer): whole=pyseer(p,i) if whole is None else whole.merge(pyseer(p,i),on='variant',how='outer')
    bcols=[c for c in whole if c.startswith('b')]; qcols=[c for c in whole if c.startswith('q')]
    whole['models']=whole[bcols].notna().sum(axis=1); whole['beta_min']=whole[bcols].min(axis=1); whole['q_max']=whole[qcols].max(axis=1)
    rows=[]; detail={}
    for _,s in sel.iterrows():
        seq=s.canonical_sequence
        if seq not in vr.index or seq not in ar.index: continue
        va=cont(vr.loc[seq],val); aa=cont(ar.loc[seq],meta); vo,vlo,vhi=orci(*va); ao,alo,ahi=orci(*aa)
        vp=float(fisher_exact([[va[0],va[1]],[va[2],va[3]]],alternative='greater').pvalue)
        src=random_meta(ar.loc[seq],meta,'source_group'); ctr=random_meta(ar.loc[seq],meta,'ISO_country_code')
        detail[str(s.candidate_id)]={'source':src,'country':ctr}
        d=s.to_dict(); d.update({'variant':seq,'validation_R_present':va[0],'validation_S_present':va[1],'validation_R_absent':va[2],'validation_S_absent':va[3],'validation_or':vo,'validation_ci_low':vlo,'validation_ci_high':vhi,'validation_p':vp,'all_or':ao,'all_ci_low':alo,'all_ci_high':ahi,'source_n':src.get('n',0),'source_or':src.get('or'),'source_ci_low':src.get('lo'),'source_ci_high':src.get('hi'),'source_p':src.get('p'),'source_I2':src.get('I2'),'country_n':ctr.get('n',0),'country_or':ctr.get('or'),'country_ci_low':ctr.get('lo'),'country_ci_high':ctr.get('hi'),'country_p':ctr.get('p')}); rows.append(d)
    res=pd.DataFrame(rows)
    if res.empty: raise RuntimeError('No selected unitig overlapped validation and all matrices')
    res['validation_q']=bh(res.validation_p); res=res.merge(whole,on='variant',how='left'); nmodels=len(a.whole_pyseer)
    res['validation_replication']=(res.validation_R_present+res.validation_S_present>=a.min_validation_present)&(res.validation_or>1)&(res.validation_ci_low>1)&(res.validation_q<=a.alpha)
    res['whole_adjusted_stable']=(res.models==nmodels)&(res.beta_min>0)&(res.q_max<=a.alpha)
    res['source_replication']=(res.source_n>=3)&(pd.to_numeric(res.source_ci_low,errors='coerce')>1)&(pd.to_numeric(res.source_p,errors='coerce')<=a.alpha)
    res['country_replication']=(res.country_n>=3)&(pd.to_numeric(res.country_ci_low,errors='coerce')>1)&(pd.to_numeric(res.country_p,errors='coerce')<=a.alpha)
    res['strict_statistical_replication']=res.validation_replication&res.whole_adjusted_stable&res.source_replication&res.country_replication
    res['score']=4*res.validation_replication.astype(int)+3*res.whole_adjusted_stable.astype(int)+4*res.source_replication.astype(int)+3*res.country_replication.astype(int)
    res=res.sort_values(['strict_statistical_replication','score','validation_q','q_max'],ascending=[False,False,True,True])
    res.to_csv(out/'ALL_UNITIG_REPLICATION_EVIDENCE.csv',index=False); strict=res[res.strict_statistical_replication].copy(); strict.to_csv(out/'STRICT_STATISTICALLY_REPLICATED_UNITIGS.csv',index=False)
    (out/'UNITIG_META_ANALYSIS_DETAILS.json').write_text(json.dumps(detail,indent=2,ensure_ascii=False,default=str)+'\n')
    summary={'n_selected':int(len(sel)),'n_evaluated':int(len(res)),'n_validation_replicated':int(res.validation_replication.sum()),'n_whole_adjusted_stable':int(res.whole_adjusted_stable.sum()),'n_source_replicated':int(res.source_replication.sum()),'n_country_replicated':int(res.country_replication.sum()),'n_strict_statistically_replicated':int(res.strict_statistical_replication.sum()),'strict_candidate_ids':strict.candidate_id.astype(str).tolist(),'status':'PORTAL_COHORT_UNITIGS_REQUIRE_KNOWN_MECHANISM_INTERSECTION' if len(strict) else 'NO_UNITIG_SURVIVED_COMPLETE_PORTAL_COHORT_GATE','boundary':'Statistical replication in the portal-residual cohort is not a novel resistance determinant. Candidates must survive the independent sequence-level known-mechanism filter, context/database/literature review, and biological validation.'}
    (out/'UNITIG_REPLICATION_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    cols=[c for c in ['candidate_id','sequence_length','p_max','q_max','validation_or','validation_ci_low','validation_ci_high','validation_q','beta_min','source_n','source_or','source_ci_low','source_ci_high','country_n','country_or','country_ci_low','country_ci_high','strict_statistical_replication'] if c in res]
    report=['# Portal-residual whole-genome unitig replication audit','',f"- Discovery-selected unitigs: **{len(sel):,}**",f"- Evaluated in untouched validation: **{len(res):,}**",f"- Held-out replicates: **{int(res.validation_replication.sum()):,}**",f"- Stable across adjusted models: **{int(res.whole_adjusted_stable.sum()):,}**",f"- Cross-source replicates: **{int(res.source_replication.sum()):,}**",f"- Cross-country replicates: **{int(res.country_replication.sum()):,}**",f"- Complete statistical gate: **{int(res.strict_statistical_replication.sum()):,}**",'', '## Claim boundary','',summary['boundary'],'','## Top evidence-ranked unitigs','',res.head(30)[cols].to_markdown(index=False)]
    (out/'UNITIG_REPLICATION_REPORT.md').write_text('\n'.join(report)+'\n')
    hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']; (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
