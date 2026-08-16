#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20260817
META = {'Unnamed: 0','protein','protein_sequence','mutant','mutated_sequence','DMS_bin_score','DMS_score_bin','label'}


def label01(s: pd.Series) -> pd.Series:
    x = s.astype(str).str.strip().str.lower()
    return x.map({'pathogenic':1,'likely pathogenic':1,'benign':0,'likely benign':0})


def rank01(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors='coerce').rank(method='average', pct=True)


def stable_bucket(text: str, modulo: int = 10) -> int:
    h = hashlib.sha256(str(text).encode()).hexdigest()
    return int(h[:12],16) % modulo


def safe_auc(y, p) -> float:
    y=np.asarray(y,dtype=float); p=np.asarray(p,dtype=float)
    m=np.isfinite(y)&np.isfinite(p)
    if m.sum()<4 or len(np.unique(y[m]))!=2: return float('nan')
    return float(roc_auc_score(y[m].astype(int),p[m]))


def bootstrap_ci(v: np.ndarray, n: int=20000) -> list[float|None]:
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if len(v)<3:return [None,None]
    rng=np.random.default_rng(SEED); idx=rng.integers(0,len(v),size=(n,len(v)))
    means=v[idx].mean(axis=1)
    return [float(np.quantile(means,.025)),float(np.quantile(means,.975))]


def family(name: str) -> str:
    n=str(name).lower()
    rules=[
      (r'polyphen','PolyPhen'),(r'sift','SIFT'),(r'varity','VARITY'),(r'bayesdel','BayesDel'),
      (r'trancepteve','TranceptEVE'),(r'^eve$','EVE'),(r'poet','PoET'),(r'gemme','GEMME'),
      (r'esm1b','ESM1b'),(r'mutationtaster','MutationTaster'),(r'mutationassessor','MutationAssessor'),
      (r'clinpred','ClinPred'),(r'dann','DANN'),(r'cadd','CADD'),(r'fathmm','FATHMM'),
      (r'list','LIST'),(r'primate','PrimateAI'),(r'mutpred','MutPred'),(r'metar','MetaRNN'),
      (r'revel','REVEL'),(r'vest','VEST'),(r'proven','PROVEAN'),(r'deogen','DEOGEN'),
      (r'gMVP'.lower(),'gMVP'),(r'mpc','MPC'),(r'lrt','LRT')]
    for pat,f in rules:
      if re.search(pat,n):return f
    return re.sub(r'[^a-z0-9]+','',n)


def zipmap(z): return {Path(n).name:n for n in z.namelist() if n.lower().endswith('.csv')}


def load_clinical(score_zip: Path) -> tuple[pd.DataFrame,list[str],dict[str,Any]]:
    frames=[]; failures=[]; headers=[]
    with zipfile.ZipFile(score_zip) as z:
      members=zipmap(z)
      for i,(fn,member) in enumerate(sorted(members.items()),1):
        try:
          with z.open(member) as f:d=pd.read_csv(f,low_memory=False)
          if not {'mutant','DMS_bin_score','protein_sequence'}.issubset(d.columns):
            failures.append({'file':fn,'reason':'missing required columns'}); continue
          model_cols=[c for c in d.columns if c not in META]
          headers.append(set(model_cols))
          d['label']=label01(d['DMS_bin_score'])
          d=d[d['label'].isin([0,1])].drop_duplicates('mutant')
          if len(d)<4 or d['label'].nunique()!=2: continue
          d['protein_file']=fn
          d['protein_group']=fn.split('.')[0]
          keep=['protein_file','protein_group','protein_sequence','mutant','label']+model_cols
          frames.append(d[keep].copy())
        except Exception as exc: failures.append({'file':fn,'reason':repr(exc)})
        if i%500==0:print(f'load {i}/{len(members)} accepted={len(frames)}',flush=True)
    if not frames:raise RuntimeError('No clinical data loaded')
    common=set.intersection(*headers) if headers else set()
    # Require column to appear in at least 95% of accepted files, rather than exact intersection only.
    counts={}
    for hs in headers:
      for c in hs:counts[c]=counts.get(c,0)+1
    models=sorted([c for c,n in counts.items() if n>=math.ceil(.95*len(headers))])
    data=pd.concat(frames,ignore_index=True)
    for c in models:
      if c not in data.columns:data[c]=np.nan
    data=data[['protein_file','protein_group','protein_sequence','mutant','label']+models]
    for c in models:
      data[c]=data.groupby('protein_file')[c].transform(rank01)
    return data,models,{'files':len(frames),'variants':len(data),'models':models,'failures':failures}


def protein_aucs(df:pd.DataFrame,pred:np.ndarray,name='prediction')->pd.DataFrame:
    t=df[['protein_file','label']].copy();t[name]=pred
    rows=[]
    for p,g in t.groupby('protein_file',sort=False):
      rows.append({'protein_file':p,'n':len(g),'n_pathogenic':int(g.label.sum()),name:safe_auc(g.label,g[name])})
    return pd.DataFrame(rows)


def train_orientation_perf(df:pd.DataFrame,models:list[str])->tuple[dict[str,int],dict[str,float]]:
    signs={};perf={}
    for m in models:
      vals=[]
      for _,g in df.groupby('protein_file'):
        vals.append(safe_auc(g.label,g[m]))
      mean=float(np.nanmean(vals)); signs[m]=1 if not np.isfinite(mean) or mean>=.5 else -1
      perf[m]=mean if signs[m]>0 else 1-mean
    return signs,perf


def select_diverse(perf:dict[str,float],k:int)->list[str]:
    best={}
    for m,p in perf.items():
      f=family(m)
      if f not in best or p>best[f][1]:best[f]=(m,p)
    return [m for m,_ in sorted(best.values(),key=lambda x:(-np.nan_to_num(x[1],nan=-9),x[0]))[:k]]


def oriented(df:pd.DataFrame,models:list[str],signs:dict[str,int])->np.ndarray:
    X=df[models].to_numpy(float)
    for j,m in enumerate(models):
      if signs[m]<0:X[:,j]=1-X[:,j]
    return X


def weights_per_protein_class(df:pd.DataFrame)->np.ndarray:
    w=np.zeros(len(df),float)
    for _,idx in df.groupby('protein_file').groups.items():
      idx=np.asarray(list(idx),int); y=df.loc[idx,'label'].to_numpy(int)
      for cls in [0,1]:
        ii=idx[y==cls]
        if len(ii):w[ii]=.5/len(ii)
    w*=len(w)/w.sum()
    return w


def features(df,models,signs):
    X=oriented(df,models,signs);miss=(~np.isfinite(X)).astype(float)
    return np.column_stack([X,miss,np.nanmean(X,axis=1),np.nanmedian(X,axis=1),np.nanstd(X,axis=1),np.isfinite(X).sum(axis=1)/max(1,X.shape[1])])

@dataclass(frozen=True)
class Spec:
    name:str;kind:str;k:int=0;power:float=1.;C:float=1.

SPECS=[
 Spec('mean5','mean',5),Spec('mean10','mean',10),Spec('mean20','mean',20),Spec('mean_all','mean',99),
 Spec('median5','median',5),Spec('median10','median',10),Spec('median20','median',20),Spec('median_all','median',99),
 Spec('weighted10_p1','weighted',10,1),Spec('weighted20_p1','weighted',20,1),Spec('weighted_all_p1','weighted',99,1),
 Spec('weighted10_p2','weighted',10,2),Spec('weighted20_p2','weighted',20,2),Spec('weighted_all_p2','weighted',99,2),
 Spec('logit10_c01','logit',10,C=.1),Spec('logit20_c01','logit',20,C=.1),Spec('logit_all_c01','logit',99,C=.1),
 Spec('logit10_c1','logit',10,C=1),Spec('logit20_c1','logit',20,C=1),Spec('logit_all_c1','logit',99,C=1),
 Spec('hgb20','hgb',20),Spec('hgb_all','hgb',99),Spec('extra20','extra',20),Spec('extra_all','extra',99)]

class Fitted:
  def __init__(self,spec,models,signs,perf,model=None):self.spec=spec;self.models=models;self.signs=signs;self.perf=perf;self.model=model
  def predict(self,df):
    X=oriented(df,self.models,self.signs)
    if self.spec.kind=='mean':return np.nanmean(X,axis=1)
    if self.spec.kind=='median':return np.nanmedian(X,axis=1)
    if self.spec.kind=='weighted':
      w=np.array([max(self.perf[m]-.5,.001)**self.spec.power for m in self.models]);w/=w.sum();ok=np.isfinite(X);W=ok*w;den=W.sum(1);p=np.full(len(df),np.nan);v=den>0;p[v]=np.nansum(X[v]*w,axis=1)/den[v];return p
    return self.model.predict_proba(features(df,self.models,self.signs))[:,1]

def fit_spec(train,all_models,spec):
    signs,perf=train_orientation_perf(train,all_models);sel=select_diverse(perf,min(spec.k,len(set(family(m) for m in all_models))))
    model=None
    if spec.kind in {'logit','hgb','extra'}:
      F=features(train,sel,signs);y=train.label.to_numpy(int);w=weights_per_protein_class(train)
      if spec.kind=='logit':model=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),LogisticRegression(C=spec.C,max_iter=5000,class_weight='balanced'))
      elif spec.kind=='hgb':model=make_pipeline(SimpleImputer(strategy='median'),HistGradientBoostingClassifier(max_iter=250,learning_rate=.05,max_leaf_nodes=15,min_samples_leaf=50,l2_regularization=2,random_state=SEED))
      else:model=make_pipeline(SimpleImputer(strategy='median'),ExtraTreesClassifier(n_estimators=400,min_samples_leaf=20,max_features=.7,n_jobs=-1,class_weight='balanced',random_state=SEED))
      model.fit(F,y,**({'logisticregression__sample_weight':w} if spec.kind=='logit' else {'histgradientboostingclassifier__sample_weight':w} if spec.kind=='hgb' else {'extratreesclassifier__sample_weight':w}))
    return Fitted(spec,sel,signs,perf,model)


def eval_mean(df,fitted,name):
    p=fitted.predict(df);tab=protein_aucs(df,p,name);return tab,float(tab[name].mean())


def compare(tab,a,b):
    d=(tab[a]-tab[b]).to_numpy(float);d=d[np.isfinite(d)]
    try:p=float(wilcoxon(d,alternative='greater').pvalue) if len(d)>=5 and np.any(d!=0) else None
    except Exception:p=None
    return {'n_proteins':len(d),'ensemble_mean_auc':float(tab[a].mean()),'comparator_mean_auc':float(tab[b].mean()),'mean_gain':float(d.mean()),'gain_95ci':bootstrap_ci(d),'fraction_improved':float((d>0).mean()),'wilcoxon_one_sided_p':p}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--clinical-score-zip',required=True);ap.add_argument('--out-dir',required=True);args=ap.parse_args()
    out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    data,models,diag=load_clinical(Path(args.clinical_score_zip))
    proteins=sorted(data.protein_file.unique());bucket={p:stable_bucket(p) for p in proteins}
    data['split']=data.protein_file.map(lambda x:'train' if bucket[x]<=5 else ('tune' if bucket[x]<=7 else 'final'))
    counts=data[['protein_file','split']].drop_duplicates().groupby('split').size().to_dict()
    train=data[data.split.eq('train')].copy();tune=data[data.split.eq('tune')].copy();final=data[data.split.eq('final')].copy()
    tune_rows=[]
    for s in SPECS:
      f=fit_spec(train,models,s);tab,mean=eval_mean(tune,f,s.name);tune_rows.append({'method':s.name,'kind':s.kind,'k':s.k,'tune_mean_protein_auc':mean,'tune_median_auc':float(tab[s.name].median()),'models':';'.join(f.models)});print(s.name,mean,flush=True)
    tuning=pd.DataFrame(tune_rows).sort_values(['tune_mean_protein_auc','method'],ascending=[False,True]);chosen_name=str(tuning.iloc[0].method);chosen=next(s for s in SPECS if s.name==chosen_name)
    train_tune=data[data.split.isin(['train','tune'])].copy();fit=fit_spec(train_tune,models,chosen);pred=fit.predict(final);ens=protein_aucs(final,pred,'ensemble')
    signs,perf=train_orientation_perf(train_tune,models);best=max(perf,key=lambda m:np.nan_to_num(perf[m],nan=-9));bp=final[best].to_numpy(float);bp=bp if signs[best]>0 else 1-bp;base=protein_aucs(final,bp,'best_selected')
    comp=ens.merge(base[['protein_file','best_selected']],on='protein_file')
    indiv=[]
    for m in models:
      p=final[m].to_numpy(float);p=p if signs[m]>0 else 1-p;t=protein_aucs(final,p,'auc');indiv.append({'model':m,'final_mean_auc':float(t.auc.mean()),'n_proteins':int(t.auc.notna().sum())})
    indiv=pd.DataFrame(indiv).sort_values('final_mean_auc',ascending=False);oracle=str(indiv.iloc[0].model);op=final[oracle].to_numpy(float);op=op if signs[oracle]>0 else 1-op;ora=protein_aucs(final,op,'oracle_best');comp=comp.merge(ora[['protein_file','oracle_best']],on='protein_file')
    vs_selected=compare(comp,'ensemble','best_selected');vs_oracle=compare(comp,'ensemble','oracle_best')
    poet=None
    if 'PoET' in models:
      pp=final.PoET.to_numpy(float);pp=pp if signs['PoET']>0 else 1-pp;pt=protein_aucs(final,pp,'PoET');comp=comp.merge(pt[['protein_file','PoET']],on='protein_file');poet=compare(comp,'ensemble','PoET')
    result={'question':'Can a protein-held-out stack of official ProteinGym clinical scores outperform the best individual predictor?','diagnostics':diag,'split_proteins':counts,'chosen_method':chosen_name,'chosen_models':fit.models,'sealed_final_proteins':int(len(ens)),'sealed_final_variants':int(len(final)),'best_individual_selected_without_final_labels':best,'posthoc_best_final_individual':oracle,'vs_selected_best':vs_selected,'vs_posthoc_oracle':vs_oracle,'vs_PoET':poet,'benchmark_advance_confirmed':bool(vs_selected['gain_95ci'][0] is not None and vs_selected['gain_95ci'][0]>0),'beats_posthoc_oracle_confirmed':bool(vs_oracle['gain_95ci'][0] is not None and vs_oracle['gain_95ci'][0]>0),'open_problem_fully_solved':False}
    tuning.to_csv(out/'tuning_results.csv',index=False);indiv.to_csv(out/'final_individual_models.csv',index=False);comp.to_csv(out/'sealed_final_protein_auc.csv',index=False)
    pd.DataFrame({'protein_file':final.protein_file,'mutant':final.mutant,'label':final.label,'ensemble_prediction':pred}).to_csv(out/'sealed_final_predictions.csv.gz',index=False,compression='gzip')
    (out/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'REPORT.md').write_text(f'''# ProteinGym clinical stack v7\n\nProteins were split before method selection into train/tune/sealed-final sets. The final labels were not used to choose the algorithm, predictors, orientation, or comparator.\n\n- Clinical proteins: {diag['files']}\n- Clinical variants: {diag['variants']:,}\n- Predictors: {len(models)}\n- Split: {counts}\n- Chosen method: `{chosen_name}`\n- Sealed-final proteins: {len(ens)}\n- Ensemble mean protein AUC: {vs_selected['ensemble_mean_auc']:.4f}\n- Best individual selected on train+tune: `{best}` = {vs_selected['comparator_mean_auc']:.4f}\n- Gain: {vs_selected['mean_gain']:+.4f}, 95% CI {vs_selected['gain_95ci']}\n- Post-hoc best final individual: `{oracle}`\n- Gain versus post-hoc oracle: {vs_oracle['mean_gain']:+.4f}, 95% CI {vs_oracle['gain_95ci']}\n- Benchmark advance confirmed: {result['benchmark_advance_confirmed']}\n\nThis is a protein-held-out benchmark experiment. It does not solve all missense pathogenicity, calibration, penetrance, or patient-level diagnosis.\n''',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
