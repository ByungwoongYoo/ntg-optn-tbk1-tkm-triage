#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import make_column_transformer
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RNG = np.random.default_rng(20260817)
META_COLS = {
    "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin", "DMS_score_bin_manual",
    "fitness", "score", "target_seq", "sequence", "wildtype", "wild_type", "wt"
}


def rank01(x: pd.Series) -> pd.Series:
    return x.rank(method="average", pct=True)


def safe_spearman(a, b) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 20:
        return float("nan")
    if np.nanstd(np.asarray(a)[m]) == 0 or np.nanstd(np.asarray(b)[m]) == 0:
        return float("nan")
    return float(spearmanr(np.asarray(a)[m], np.asarray(b)[m]).statistic)


def list_csv_members(z: zipfile.ZipFile) -> list[str]:
    return [n for n in z.namelist() if n.lower().endswith('.csv') and not n.endswith('/')]


def detect_score_column(df: pd.DataFrame) -> str | None:
    numeric = []
    for c in df.columns:
        if str(c).lower() in META_COLS or str(c).lower().startswith('dms_'):
            continue
        s = pd.to_numeric(df[c], errors='coerce')
        n = int(s.notna().sum())
        if n:
            numeric.append((n, float(s.nunique(dropna=True)), str(c)))
    if not numeric:
        return None
    numeric.sort(reverse=True)
    return numeric[0][2]


def normalize_model_name(path: str) -> str:
    parts = Path(path).parts
    if len(parts) < 2:
        return "unknown"
    # The score archive convention is root/model/assay.csv. Use the parent directory.
    return parts[-2]


def choose_models(model_coverage: dict[str, int], n_assays: int, max_models: int = 12) -> list[str]:
    desired = [
        r"AIDO.*RAG", r"VenusREM", r"ProSST.*2048", r"S3F.*MSA", r"PoET.*200",
        r"ESM3", r"VespaG", r"SaProt.*650", r"TranceptEVE.*L", r"GEMME",
        r"ProteinMPNN", r"ESM.?1v", r"Tranception.*L", r"MSA.*Transformer"
    ]
    selected: list[str] = []
    names = sorted(model_coverage, key=lambda x: (-model_coverage[x], x.lower()))
    for pat in desired:
        cand = [x for x in names if re.search(pat, x, re.I) and model_coverage[x] >= max(10, int(0.55*n_assays))]
        if cand:
            selected.append(cand[0])
    # Fill with high-coverage models, avoiding near-duplicate names.
    for x in names:
        if len(selected) >= max_models:
            break
        if model_coverage[x] < max(10, int(0.75*n_assays)):
            continue
        stem = re.sub(r"\([^)]*\)|[_\- ]?(large|medium|small|base|ensemble|single)$", "", x.lower()).strip()
        if any(stem and stem in y.lower() for y in selected):
            continue
        selected.append(x)
    return selected[:max_models]


def build_member_maps(score_zip: zipfile.ZipFile, assay_files: set[str]):
    members = list_csv_members(score_zip)
    coverage: dict[str, set[str]] = defaultdict(set)
    member_map: dict[tuple[str, str], str] = {}
    for n in members:
        base = Path(n).name
        if base not in assay_files:
            continue
        model = normalize_model_name(n)
        coverage[model].add(base)
        member_map[(model, base)] = n
    return {m: len(v) for m, v in coverage.items()}, member_map


def read_csv_from_zip(z: zipfile.ZipFile, member: str) -> pd.DataFrame:
    with z.open(member) as f:
        return pd.read_csv(f, low_memory=False)


def assay_member_map(dms_zip: zipfile.ZipFile) -> dict[str, str]:
    out = {}
    for n in list_csv_members(dms_zip):
        out[Path(n).name] = n
    return out


def load_assay(
    meta: pd.Series,
    dms_zip: zipfile.ZipFile,
    dms_member: str,
    score_zip: zipfile.ZipFile,
    score_member_map: dict[tuple[str, str], str],
    models: list[str],
) -> pd.DataFrame | None:
    d = read_csv_from_zip(dms_zip, dms_member)
    if 'mutant' not in d.columns or 'DMS_score' not in d.columns:
        return None
    d = d[['mutant', 'DMS_score']].copy()
    d['mutant'] = d['mutant'].astype(str)
    d = d[~d['mutant'].str.contains(':', regex=False)]
    d['DMS_score'] = pd.to_numeric(d['DMS_score'], errors='coerce')
    d = d.dropna(subset=['DMS_score']).drop_duplicates('mutant')
    if len(d) < 100:
        return None
    for model in models:
        member = score_member_map.get((model, meta['DMS_filename']))
        if member is None:
            d[model] = np.nan
            continue
        try:
            s = read_csv_from_zip(score_zip, member)
            if 'mutant' not in s.columns:
                d[model] = np.nan
                continue
            col = detect_score_column(s)
            if col is None:
                d[model] = np.nan
                continue
            ss = s[['mutant', col]].copy()
            ss['mutant'] = ss['mutant'].astype(str)
            ss[col] = pd.to_numeric(ss[col], errors='coerce')
            ss = ss.drop_duplicates('mutant').rename(columns={col: model})
            d = d.merge(ss, on='mutant', how='left')
        except Exception:
            d[model] = np.nan
    d['DMS_id'] = meta['DMS_id']
    d['UniProt_ID'] = meta['UniProt_ID']
    d['taxon'] = meta['taxon']
    d['selection_type'] = meta.get('coarse_selection_type', meta.get('selection_type', 'unknown'))
    d['msa_depth'] = meta.get('MSA_Neff_L_category', 'unknown')
    d['target_rank'] = rank01(d['DMS_score'])
    for model in models:
        d[model + '__rank'] = rank01(d[model])
    return d


def bootstrap_ci(values: np.ndarray, n=10000, seed=20260817):
    values = values[np.isfinite(values)]
    if len(values) < 3:
        return [None, None]
    rng = np.random.default_rng(seed)
    b = np.empty(n)
    for i in range(n):
        b[i] = np.mean(rng.choice(values, size=len(values), replace=True))
    return [float(np.quantile(b, .025)), float(np.quantile(b, .975))]


def orient_features(train: pd.DataFrame, models: list[str]) -> dict[str, int]:
    signs = {}
    for m in models:
        vals = []
        for _, g in train.groupby('DMS_id'):
            vals.append(safe_spearman(g[m+'__rank'].to_numpy(), g['target_rank'].to_numpy()))
        med = np.nanmedian(vals)
        signs[m] = 1 if not np.isfinite(med) or med >= 0 else -1
    return signs


def make_feature_frame(df: pd.DataFrame, models: list[str], signs: dict[str, int]) -> pd.DataFrame:
    x = pd.DataFrame(index=df.index)
    for m in models:
        r = df[m+'__rank'].astype(float)
        x[m] = r if signs[m] == 1 else (1-r)
        x[m+'__missing'] = r.isna().astype(float)
    x['model_mean'] = x[models].mean(axis=1)
    x['model_median'] = x[models].median(axis=1)
    x['model_sd'] = x[models].std(axis=1)
    x['n_models'] = x[models].notna().sum(axis=1)
    return x


def assay_balanced_sample(df: pd.DataFrame, max_per_assay=1500) -> pd.DataFrame:
    chunks=[]
    for _,g in df.groupby('DMS_id'):
        if len(g)>max_per_assay:
            chunks.append(g.sample(max_per_assay, random_state=20260817))
        else:
            chunks.append(g)
    return pd.concat(chunks, ignore_index=True)


def score_predictions(df: pd.DataFrame, pred_cols: list[str], models: list[str]):
    rows=[]
    for assay,g in df.groupby('DMS_id'):
        row={'DMS_id':assay,'UniProt_ID':g['UniProt_ID'].iloc[0],'n_mutants':len(g),'taxon':g['taxon'].iloc[0],'selection_type':g['selection_type'].iloc[0],'msa_depth':g['msa_depth'].iloc[0]}
        y=g['target_rank'].to_numpy()
        for c in pred_cols:
            row[c]=safe_spearman(g[c].to_numpy(),y)
        for m in models:
            row['base__'+m]=safe_spearman(g[m+'__pred'].to_numpy(),y)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--dms-zip',required=True)
    ap.add_argument('--score-zip',required=True)
    ap.add_argument('--reference',required=True)
    ap.add_argument('--out-dir',required=True)
    ap.add_argument('--human-only',action='store_true')
    ap.add_argument('--max-assays',type=int,default=0)
    ap.add_argument('--max-models',type=int,default=12)
    args=ap.parse_args()
    out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
    ref=pd.read_csv(args.reference)
    if args.human_only:
        ref=ref[ref['taxon'].astype(str).str.lower().eq('human')].copy()
    ref=ref[(ref['DMS_number_single_mutants']>=100)&(~ref['includes_multiple_mutants'].astype(bool))].copy()
    if args.max_assays and len(ref)>args.max_assays:
        # stratified deterministic sample by selection type and MSA depth
        ref=ref.sort_values(['coarse_selection_type','MSA_Neff_L_category','DMS_id']).groupby(['coarse_selection_type','MSA_Neff_L_category'],group_keys=False).head(max(1,args.max_assays//max(1,ref.groupby(['coarse_selection_type','MSA_Neff_L_category']).ngroups)))
        ref=ref.head(args.max_assays)
    with zipfile.ZipFile(args.dms_zip) as dz, zipfile.ZipFile(args.score_zip) as sz:
        dmap=assay_member_map(dz)
        ref=ref[ref['DMS_filename'].isin(dmap)].copy()
        coverage,smap=build_member_maps(sz,set(ref['DMS_filename']))
        models=choose_models(coverage,len(ref),args.max_models)
        (out/'model_inventory.json').write_text(json.dumps({'coverage':coverage,'selected_models':models,'n_reference_assays':len(ref)},indent=2),encoding='utf-8')
        assays=[]; failures=[]
        for i,(_,meta) in enumerate(ref.iterrows(),1):
            try:
                d=load_assay(meta,dz,dmap[meta['DMS_filename']],sz,smap,models)
                if d is not None and d[[m+'__rank' for m in models]].notna().sum().sum()>0:
                    assays.append(d)
                else:
                    failures.append({'DMS_id':meta['DMS_id'],'reason':'no usable merged data'})
            except Exception as e:
                failures.append({'DMS_id':meta['DMS_id'],'reason':repr(e)})
            print(f'loaded {i}/{len(ref)} usable={len(assays)}',flush=True)
    if len(assays)<15:
        raise RuntimeError(f'Only {len(assays)} usable assays; failures={failures[:5]}')
    data=pd.concat(assays,ignore_index=True)
    # Require a reasonable number of model scores per mutation.
    rank_cols=[m+'__rank' for m in models]
    data=data[data[rank_cols].notna().sum(axis=1)>=max(3,len(models)//2)].copy()
    groups=data[['DMS_id','UniProt_ID']].drop_duplicates()
    n_splits=min(5,groups['UniProt_ID'].nunique())
    gkf=GroupKFold(n_splits=n_splits)
    assay_index=groups.reset_index(drop=True)
    fold_assignment={}
    for fold,(_,te) in enumerate(gkf.split(assay_index,groups=assay_index['UniProt_ID']),1):
        for a in assay_index.iloc[te]['DMS_id']:
            fold_assignment[a]=fold
    data['fold']=data['DMS_id'].map(fold_assignment)
    data['pred_mean']=np.nan; data['pred_median']=np.nan; data['pred_ridge']=np.nan; data['pred_hgb']=np.nan; data['pred_extratrees']=np.nan
    for m in models: data[m+'__pred']=np.nan
    fold_info=[]
    for fold in sorted(data['fold'].dropna().unique()):
        tr=data[data['fold']!=fold].copy(); te=data[data['fold']==fold].copy()
        signs=orient_features(tr,models)
        Xtr=make_feature_frame(tr,models,signs); Xte=make_feature_frame(te,models,signs)
        # Save oriented base predictions for fair held-out evaluation.
        for m in models:
            data.loc[te.index,m+'__pred']=Xte[m]
        data.loc[te.index,'pred_mean']=Xte[models].mean(axis=1)
        data.loc[te.index,'pred_median']=Xte[models].median(axis=1)
        train_sample=assay_balanced_sample(pd.concat([tr[['DMS_id','target_rank']],Xtr],axis=1),max_per_assay=1200)
        ytr=train_sample['target_rank'].to_numpy()
        feats=[c for c in Xtr.columns]
        Xs=train_sample[feats]
        ridge=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),Ridge(alpha=25.0))
        ridge.fit(Xs,ytr)
        data.loc[te.index,'pred_ridge']=ridge.predict(Xte[feats])
        hgb=make_pipeline(SimpleImputer(strategy='median'),HistGradientBoostingRegressor(max_iter=200,learning_rate=.05,max_leaf_nodes=15,l2_regularization=2.0,random_state=20260817))
        hgb.fit(Xs,ytr)
        data.loc[te.index,'pred_hgb']=hgb.predict(Xte[feats])
        et=make_pipeline(SimpleImputer(strategy='median'),ExtraTreesRegressor(n_estimators=200,min_samples_leaf=25,max_features=.8,n_jobs=-1,random_state=20260817))
        et.fit(Xs,ytr)
        data.loc[te.index,'pred_extratrees']=et.predict(Xte[feats])
        fold_info.append({'fold':int(fold),'train_assays':int(tr['DMS_id'].nunique()),'test_assays':int(te['DMS_id'].nunique()),'signs':signs,'train_mutants_sampled':len(train_sample)})
    pred_cols=['pred_mean','pred_median','pred_ridge','pred_hgb','pred_extratrees']
    assay_scores=score_predictions(data,pred_cols,models)
    # Primary summary: equal weight per assay. Secondary: equal weight per UniProt.
    model_cols=pred_cols+['base__'+m for m in models]
    summary=[]
    for c in model_cols:
        vals=assay_scores[c].to_numpy(dtype=float)
        summary.append({'model':c,'mean_assay_spearman':float(np.nanmean(vals)),'median_assay_spearman':float(np.nanmedian(vals)),'n_assays':int(np.isfinite(vals).sum()),'mean_assay_95ci':bootstrap_ci(vals)})
    summary_df=pd.DataFrame(summary).sort_values('mean_assay_spearman',ascending=False)
    # Protein-level averaging avoids overweighting proteins with multiple assays.
    protein_scores=assay_scores.groupby('UniProt_ID')[model_cols].mean(numeric_only=True)
    protein_summary=[]
    for c in model_cols:
        v=protein_scores[c].to_numpy(float)
        protein_summary.append({'model':c,'mean_uniprot_spearman':float(np.nanmean(v)),'median_uniprot_spearman':float(np.nanmedian(v)),'n_uniprot':int(np.isfinite(v).sum()),'mean_uniprot_95ci':bootstrap_ci(v)})
    protein_summary_df=pd.DataFrame(protein_summary).sort_values('mean_uniprot_spearman',ascending=False)
    # Paired comparison of best deployable ensemble against best individual base model.
    ens_candidates=summary_df[summary_df['model'].isin(pred_cols)]
    base_candidates=summary_df[summary_df['model'].str.startswith('base__')]
    best_ens=ens_candidates.iloc[0]['model']; best_base=base_candidates.iloc[0]['model']
    paired=(assay_scores[best_ens]-assay_scores[best_base]).to_numpy(float)
    paired=paired[np.isfinite(paired)]
    result={
        'scope':'human-only' if args.human_only else 'all taxa',
        'n_assays':int(assay_scores['DMS_id'].nunique()),
        'n_uniprot':int(assay_scores['UniProt_ID'].nunique()),
        'n_mutants':int(len(data)),
        'selected_models':models,
        'best_ensemble':best_ens,
        'best_ensemble_mean_assay_spearman':float(summary_df.set_index('model').loc[best_ens,'mean_assay_spearman']),
        'best_base':best_base,
        'best_base_mean_assay_spearman':float(summary_df.set_index('model').loc[best_base,'mean_assay_spearman']),
        'paired_mean_difference':float(np.mean(paired)),
        'paired_difference_95ci':bootstrap_ci(paired),
        'fraction_assays_ensemble_better':float(np.mean(paired>0)),
        'folds':fold_info,
        'failures':failures,
    }
    data[['DMS_id','UniProt_ID','mutant','DMS_score','target_rank','fold']+pred_cols+[m+'__pred' for m in models]].to_csv(out/'heldout_predictions.csv.gz',index=False,compression='gzip')
    assay_scores.to_csv(out/'assay_scores.csv',index=False)
    summary_df.to_csv(out/'assay_summary.csv',index=False)
    protein_summary_df.to_csv(out/'uniprot_summary.csv',index=False)
    (out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    report=f"""# ProteinGym cross-protein missense-effect ensemble experiment

## Leakage control

All mutation-level predictions were evaluated out of fold. Folds were grouped by UniProt ID, so no assay from a held-out protein contributed labels to the meta-model or score orientation used for that protein.

## Scope

- Assays: {result['n_assays']}
- Unique UniProt proteins: {result['n_uniprot']}
- Single amino-acid substitutions: {result['n_mutants']:,}
- Public zero-shot models combined: {len(models)}

## Primary result

- Best ensemble: `{best_ens}`
- Mean assay Spearman: {result['best_ensemble_mean_assay_spearman']:.4f}
- Best individual public score in this merged subset: `{best_base}`
- Mean assay Spearman: {result['best_base_mean_assay_spearman']:.4f}
- Paired mean difference: {result['paired_mean_difference']:+.4f}
- Paired bootstrap 95% CI: [{result['paired_difference_95ci'][0]:+.4f}, {result['paired_difference_95ci'][1]:+.4f}]
- Fraction of assays improved: {100*result['fraction_assays_ensemble_better']:.1f}%

## Interpretation boundary

This experiment tests whether a cross-protein ensemble of existing public predictors improves held-out DMS ranking. It does not solve all human missense effects, establish clinical pathogenicity, or create new experimental labels. A positive result is a benchmark advance; a null result shows that simple score stacking is not the missing solution.
"""
    (out/'REPORT.md').write_text(report,encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    main()
