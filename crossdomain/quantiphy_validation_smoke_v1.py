#!/usr/bin/env python3
from __future__ import annotations
import json,math,os,re
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from datasets import load_dataset
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

TH=np.array([.10,.20,.30,.40,.50,.60,.70,.80,.90,.95])
NUM=re.compile(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?')

def mra(y,p):
    y=np.asarray(y,float);p=np.asarray(p,float);valid=np.isfinite(y)&np.isfinite(p)&(np.abs(y)>0)
    if not valid.any():return 0.0
    rel=np.abs(p[valid]-y[valid])/np.abs(y[valid]);return float(np.mean(rel[:,None] < (1-TH)[None,:]))
def get_y(df):
    for c in ['answer','ground_truth','gt','value','target','target_value','parsed_value']:
        if c in df.columns:
            v=pd.to_numeric(df[c],errors='coerce')
            if v.notna().mean()>.8:return c,v.to_numpy(float)
    raise RuntimeError('No numeric ground truth column found: '+repr(list(df.columns)))
def numbers(s):
    return [float(x) for x in NUM.findall(str(s or ''))]
def category(df):
    # Official category convention is inference type static/dynamic x dimension 2D/3D; preserve raw values for diagnostics.
    cols=[c for c in ['inference_type','video_type','depth_info'] if c in df.columns]
    return df[cols].astype(str).agg('|'.join,axis=1) if cols else pd.Series(['all']*len(df))
def eval_pred(df,y,p):
    cats=category(df);vals={}
    for c in sorted(cats.unique()):
        ix=(cats==c).to_numpy(); vals[str(c)]=mra(y[ix],p[ix])
    return {'overall_mra':mra(y,p),'category_mra':vals,'finite_fraction':float(np.isfinite(p).mean())}
def main():
    out=Path(os.environ.get('OUT_DIR','artifact/quantiphy_v1'));out.mkdir(parents=True,exist_ok=True)
    ds=load_dataset('PaulineLi/QuantiPhy-validation',split='validation');df=ds.to_pandas();yc,y=get_y(df)
    # Save schema and first rows without video binary objects.
    safe=df.copy()
    for c in safe.columns:
        if safe[c].map(lambda x:isinstance(x,(bytes,bytearray,dict,list))).any():safe[c]=safe[c].map(lambda x:str(x)[:1000])
    safe.head(30).to_csv(out/'sample30.csv',index=False)
    res={'n':len(df),'columns':list(df.columns),'ground_truth_column':yc,'baselines':{}}
    # H0: constant geometric median-ish positive target (LOO median by full data is optimistic but a scale sanity bound, marked exploratory).
    pos=np.abs(y[np.isfinite(y)&(y!=0)]);global_med=float(np.exp(np.median(np.log(pos)))) if len(pos) else 1.0
    res['baselines']['global_geometric_median_exploratory']=eval_pred(df,y,np.full(len(df),global_med))
    # H1/H2: numbers explicitly supplied in question or prior, without video.
    for field in ['prior','question']:
        if field not in df.columns:continue
        ns=[numbers(x) for x in df[field]]
        for which in ('first','last'):
            p=np.array([(z[0] if which=='first' and z else z[-1] if z else np.nan) for z in ns],float)
            res['baselines'][f'{field}_{which}_number']=eval_pred(df,y,p)
    # H3: leave-one-video-out nearest-neighbor retrieval on question+prior text. Uses validation labels intentionally; measures repeated textual structure, not blind generalization.
    text=(df.get('question',pd.Series(['']*len(df))).astype(str)+' [PRIOR] '+df.get('prior',pd.Series(['']*len(df))).astype(str)).tolist()
    V=TfidfVectorizer(ngram_range=(1,2),min_df=1,sublinear_tf=True).fit_transform(text);sim=cosine_similarity(V)
    np.fill_diagonal(sim,-1)
    vids=df.get('video_id',pd.Series(np.arange(len(df)))).astype(str).to_numpy();pred=np.full(len(df),np.nan);nn=[]
    for i in range(len(df)):
        order=np.argsort(-sim[i]);j=next((j for j in order if vids[j]!=vids[i] and np.isfinite(y[j])),None)
        if j is not None:pred[i]=y[j];nn.append({'i':i,'j':int(j),'similarity':float(sim[i,j]),'y':float(y[i]),'pred':float(y[j])})
    res['baselines']['leave_one_video_out_text_nearest_neighbor_exploratory']=eval_pred(df,y,pred)
    (out/'nn_pairs.json').write_text(json.dumps(nn,indent=2),encoding='utf-8')
    # H4: category-specific geometric median leave-one-out group prior. Again exploratory diagnostic only.
    cats=category(df);p=np.full(len(df),np.nan)
    for i in range(len(df)):
        vals=np.abs(y[(cats==cats.iloc[i]).to_numpy() & np.isfinite(y) & (y!=0)])
        if len(vals):p[i]=float(np.exp(np.median(np.log(vals))))
    res['baselines']['category_geometric_median_exploratory']=eval_pred(df,y,p)
    res['best_baseline']=max(res['baselines'].items(),key=lambda kv:kv[1]['overall_mra'])
    res['success_gate']='Only advance this track if a no-paid-compute method can materially approach or exceed the published ~53.1 frontier-model MRA on validation without target leakage. Label-retrieval baselines are diagnostic and cannot satisfy the gate.'
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
