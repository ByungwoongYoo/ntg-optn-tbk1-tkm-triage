#!/usr/bin/env python3
"""Fastest decisive audit: nested family-held-out fixed ensemble versus best individual."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import proteingym_family_gating_v12 as core

core.FIXED_SPECS = [
    core.FixedSpec("mean3", "mean", 3),
    core.FixedSpec("mean5", "mean", 5),
    core.FixedSpec("weighted3_p1", "weighted", 3, 1.0),
    core.FixedSpec("weighted5_p1", "weighted", 5, 1.0),
]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--clinical-score-zip',required=True);ap.add_argument('--sequence-inventory',required=True);ap.add_argument('--edge',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d,inv,diag=core.load_common_scores(Path(a.clinical_score_zip));seq=pd.read_csv(a.sequence_inventory)
    cmap,cdiag=core.parse_cluster_edges(seq['header'].astype(str).tolist(),Path(a.edge));d['cluster_id']=d['header'].map(cmap)
    omap=core.assign_cluster_folds(d,'cluster_id',5);d['outer_fold']=d['cluster_id'].map(omap)
    parts=[];folds=[];tunes=[]
    for f in range(5):
        tr=d[d.outer_fold!=f].copy();te=d[d.outer_fold==f].copy()
        spec,tab,_=core.choose_spec(tr,'fixed',3);fixed=core.fit_fixed(tr,spec);best,signs=core.fit_best_individual(tr)
        pf=core.per_protein_auc(te,fixed.predict(te),'fixed');pb=core.per_protein_auc(te,core.predict_best_individual(te,best,signs),'best')
        p=pf.merge(pb[['protein_file','best']],on='protein_file');p['outer_fold']=f;parts.append(p)
        folds.append({'fold':f,'spec':spec.name,'best':best,'proteins':len(p),'fixed':float(p.fixed.mean()),'best_auc':float(p.best.mean()),'gain':float((p.fixed-p.best).mean())})
        tunes.append(tab.assign(outer_fold=f))
        print(json.dumps(folds[-1]),flush=True)
    allp=pd.concat(parts,ignore_index=True);summ=core.paired_summary(allp,'fixed','best')
    result={'load':diag,'cluster':cdiag,'folds':folds,'comparison':summ,'family_heldout_fixed_signal_confirmed':summ['confirmed'],'grand_problem_fully_solved':False}
    allp.to_csv(out/'protein_auc.csv',index=False);pd.concat(tunes).to_csv(out/'tuning.csv',index=False);(out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (out/'REPORT.md').write_text(f"# Family-held-out fixed-only audit\n\n- Fixed ensemble gain: {summ['mean_difference']:+.6f}\n- Cluster-bootstrap 95% CI: {summ['cluster_bootstrap_95ci']}\n- Confirmed: {summ['confirmed']}\n",encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
