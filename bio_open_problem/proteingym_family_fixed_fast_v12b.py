#!/usr/bin/env python3
"""Vectorized family-held-out fixed ensemble audit for ProteinGym clinical variants."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import proteingym_family_gating_v12 as core

SPECS=[
 core.FixedSpec('mean3','mean',3),core.FixedSpec('mean5','mean',5),
 core.FixedSpec('weighted3_p1','weighted',3,1.0),core.FixedSpec('weighted5_p1','weighted',5,1.0)
]

def fast_cluster_ci(tab:pd.DataFrame,n=20000,seed=20260817):
    g=tab.groupby('cluster_id')['diff'].agg(['sum','count'])
    sums=g['sum'].to_numpy(float);counts=g['count'].to_numpy(float);k=len(g)
    rng=np.random.default_rng(seed);vals=[];batch=500
    for start in range(0,n,batch):
        b=min(batch,n-start);idx=rng.integers(0,k,size=(b,k))
        vals.append(sums[idx].sum(axis=1)/counts[idx].sum(axis=1))
    v=np.concatenate(vals);return [float(np.quantile(v,.025)),float(np.quantile(v,.975))]

def choose_fixed(tr):
    old=core.FIXED_SPECS;core.FIXED_SPECS=SPECS
    try:return core.choose_spec(tr,'fixed',3)
    finally:core.FIXED_SPECS=old

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--clinical-score-zip',required=True);ap.add_argument('--sequence-inventory',required=True);ap.add_argument('--edge',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d,inv,diag=core.load_common_scores(Path(a.clinical_score_zip));seq=pd.read_csv(a.sequence_inventory)
    cmap,cdiag=core.parse_cluster_edges(seq['header'].astype(str).tolist(),Path(a.edge));d['cluster_id']=d['header'].map(cmap)
    omap=core.assign_cluster_folds(d,'cluster_id',5);d['outer_fold']=d['cluster_id'].map(omap)
    parts=[];folds=[];tunes=[]
    for f in range(5):
        tr=d[d.outer_fold!=f].copy();te=d[d.outer_fold==f].copy();spec,tab,_=choose_fixed(tr)
        fixed=core.fit_fixed(tr,spec);best,signs=core.fit_best_individual(tr)
        pf=core.per_protein_auc(te,fixed.predict(te),'fixed');pb=core.per_protein_auc(te,core.predict_best_individual(te,best,signs),'best')
        p=pf.merge(pb[['protein_file','best']],on='protein_file');p['outer_fold']=f;parts.append(p);tunes.append(tab.assign(outer_fold=f))
        folds.append({'fold':f,'spec':spec.name,'best':best,'proteins':len(p),'fixed':float(p.fixed.mean()),'best_auc':float(p.best.mean()),'gain':float((p.fixed-p.best).mean())})
        print(json.dumps(folds[-1]),flush=True)
    allp=pd.concat(parts,ignore_index=True);allp['diff']=allp.fixed-allp.best;ci=fast_cluster_ci(allp)
    cm=allp.groupby('cluster_id').diff.mean().to_numpy(float)
    from scipy.stats import wilcoxon
    pval=float(wilcoxon(cm,alternative='greater').pvalue) if len(cm)>=5 and np.any(cm!=0) else None
    comp={'n_proteins':len(allp),'n_clusters':int(allp.cluster_id.nunique()),'fixed_mean':float(allp.fixed.mean()),'best_mean':float(allp.best.mean()),'gain':float(allp['diff'].mean()),'cluster_bootstrap_95ci':ci,'fraction_improved':float((allp['diff']>1e-12).mean()),'fraction_tied':float((allp['diff'].abs()<=1e-12).mean()),'cluster_wilcoxon_p':pval,'confirmed':bool(ci[0]>0)}
    result={'load':diag,'cluster':cdiag,'folds':folds,'comparison':comp,'family_heldout_fixed_signal_confirmed':comp['confirmed'],'grand_problem_fully_solved':False}
    allp.to_csv(out/'protein_auc.csv',index=False);pd.concat(tunes).to_csv(out/'tuning.csv',index=False);(out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8');(out/'REPORT.md').write_text(f"# Family-held-out fixed audit v12b\n\n- Gain: {comp['gain']:+.6f}\n- Cluster-bootstrap 95% CI: {ci}\n- Confirmed: {comp['confirmed']}\n",encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
