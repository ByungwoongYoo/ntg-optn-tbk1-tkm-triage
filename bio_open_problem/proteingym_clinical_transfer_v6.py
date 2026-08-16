#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

from proteingym_wide_v2 import bootstrap_mean_ci, load_data, model_train_stats, rank01, stable_bucket

SEED = 20260817
META = {
    'Unnamed: 0', 'protein', 'protein_sequence', 'mutant', 'mutated_sequence',
    'DMS_bin_score', 'DMS_score_bin', 'label'
}


def normalize_name(name: str) -> str:
    x = re.sub(r'[^a-z0-9]+', '', str(name).lower())
    aliases = {
        'eveensemble': 'eve', 'evesingle': 'eve', 'eve': 'eve',
        'tranceptevelarge': 'tranceptevel', 'tranceptevel': 'tranceptevel',
        'poet200m': 'poet', 'poet': 'poet', 'gemme': 'gemme', 'esm1b': 'esm1b',
    }
    return aliases.get(x, x)


def label01(s: pd.Series) -> pd.Series:
    x = s.astype(str).str.strip().str.lower()
    return x.map({'pathogenic': 1, 'likely pathogenic': 1, 'benign': 0, 'likely benign': 0})


def safe_auc(y, p) -> float:
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    m = np.isfinite(y) & np.isfinite(p)
    if m.sum() < 4 or len(np.unique(y[m])) != 2:
        return float('nan')
    return float(roc_auc_score(y[m].astype(int), p[m]))


def sha_seq(seq: str) -> str:
    return hashlib.sha256(str(seq).strip().upper().encode()).hexdigest()


def zipmap(z: zipfile.ZipFile) -> dict[str, str]:
    return {Path(n).name: n for n in z.namelist() if n.lower().endswith('.csv')}


def common_model_pairs(dms_models: list[str], clinical_columns: list[str], perf: dict[str, float]) -> list[tuple[str, str]]:
    dm: dict[str, list[str]] = {}
    cm: dict[str, list[str]] = {}
    for m in dms_models:
        dm.setdefault(normalize_name(m), []).append(m)
    for c in clinical_columns:
        if c not in META:
            cm.setdefault(normalize_name(c), []).append(c)
    pairs = []
    for key in sorted(set(dm).intersection(cm)):
        d = max(dm[key], key=lambda m: np.nan_to_num(perf.get(m, np.nan), nan=-9))
        pairs.append((d, cm[key][0]))
    return sorted(pairs, key=lambda x: -np.nan_to_num(perf.get(x[0], np.nan), nan=-9))


def load_clinical(score_zip_path: Path, pairs: list[tuple[str, str]], excluded_seq_hashes: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    frames = []
    failures = []
    exact_overlap_files = []
    with zipfile.ZipFile(score_zip_path) as z:
        members = zipmap(z)
        for i, (fn, member) in enumerate(sorted(members.items()), 1):
            try:
                with z.open(member) as f:
                    d = pd.read_csv(f, low_memory=False)
                needed = ['mutant', 'DMS_bin_score', 'protein_sequence'] + [c for _, c in pairs if c in d.columns]
                if not {'mutant', 'DMS_bin_score', 'protein_sequence'}.issubset(d.columns):
                    failures.append({'file': fn, 'reason': 'missing label/mutant/sequence'})
                    continue
                if len(needed) < 6:
                    failures.append({'file': fn, 'reason': f'only {len(needed)-3} common predictors'})
                    continue
                d = d[needed].copy()
                d['label'] = label01(d['DMS_bin_score'])
                d = d[d['label'].isin([0, 1])].drop_duplicates('mutant')
                if len(d) < 4 or d['label'].nunique() != 2:
                    continue
                seq_hash = sha_seq(d['protein_sequence'].iloc[0])
                d['protein_file'] = fn
                d['sequence_hash'] = seq_hash
                d['strict_external'] = int(seq_hash not in excluded_seq_hashes)
                if seq_hash in excluded_seq_hashes:
                    exact_overlap_files.append(fn)
                frames.append(d)
            except Exception as exc:
                failures.append({'file': fn, 'reason': repr(exc)})
            if i % 500 == 0:
                print(f'clinical {i}/{len(members)} accepted={len(frames)}', flush=True)
    if not frames:
        raise RuntimeError(f'No clinical proteins loaded: {failures[:10]}')
    return pd.concat(frames, ignore_index=True), {
        'accepted_files': len(frames), 'failures': failures,
        'exact_sequence_overlap_files': exact_overlap_files,
    }


def auc_table(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for protein, g in df.groupby('protein_file', sort=False):
        row = {'protein_file': protein, 'n': len(g), 'n_pathogenic': int(g['label'].sum()), 'strict_external': int(g['strict_external'].iloc[0])}
        for c in cols:
            row[c] = safe_auc(g['label'], g[c])
        rows.append(row)
    return pd.DataFrame(rows)


def compare(tab: pd.DataFrame, a: str, b: str) -> dict[str, Any]:
    d = (tab[a] - tab[b]).to_numpy(float)
    d = d[np.isfinite(d)]
    try:
        p = float(wilcoxon(d, alternative='greater').pvalue) if len(d) >= 5 and np.any(d != 0) else None
    except Exception:
        p = None
    return {
        'n_proteins': int(len(d)),
        'ensemble_mean_auc': float(tab[a].mean()),
        'comparator_mean_auc': float(tab[b].mean()),
        'mean_gain': float(np.mean(d)),
        'gain_95ci': bootstrap_mean_ci(d),
        'fraction_improved': float(np.mean(d > 0)),
        'wilcoxon_one_sided_p': p,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--dms-score-zip', required=True)
    ap.add_argument('--dms-reference', required=True)
    ap.add_argument('--clinical-score-zip', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    dms, dms_models, dms_diag = load_data(Path(args.dms_score_zip), Path(args.dms_reference), human_only=False)
    human = dms[dms['taxon'].astype(str).str.lower().eq('human')].copy()
    proteins = sorted(human['UniProt_ID'].unique())
    buckets = {p: stable_bucket(p) for p in proteins}
    human['split'] = human['UniProt_ID'].map(lambda x: 'train' if buckets[x] <= 5 else ('tune' if buckets[x] <= 7 else 'final'))
    counts = human[['DMS_id','UniProt_ID','split']].drop_duplicates().groupby('split').agg(assays=('DMS_id','nunique'), proteins=('UniProt_ID','nunique')).to_dict('index')
    if any(counts.get(s,{}).get('proteins',0) < 6 for s in ['train','tune','final']):
        assignments = {}
        for i,p in enumerate(proteins):
            q = i/max(1,len(proteins)); assignments[p] = 'train' if q < .6 else ('tune' if q < .8 else 'final')
        human['split'] = human['UniProt_ID'].map(assignments)
    train_tune = human[human['split'].isin(['train','tune'])].copy()
    signs, perf = model_train_stats(train_tune, dms_models)

    ref = pd.read_csv(args.dms_reference)
    ids = set(train_tune['DMS_id'].astype(str))
    train_seq_hashes = {sha_seq(x) for x in ref.loc[ref['DMS_id'].astype(str).isin(ids), 'target_seq'].dropna()}

    with zipfile.ZipFile(args.clinical_score_zip) as z:
        first = next(n for n in z.namelist() if n.lower().endswith('.csv'))
        with z.open(first) as f:
            header = pd.read_csv(f, nrows=0).columns.tolist()
    pairs = common_model_pairs(dms_models, header, perf)
    if len(pairs) < 3:
        raise RuntimeError(f'Only {len(pairs)} common predictors: {pairs}')
    print('common predictors', pairs, flush=True)

    clinical, clinical_diag = load_clinical(Path(args.clinical_score_zip), pairs, train_seq_hashes)
    dms_names = [d for d,_ in pairs]
    weights = np.array([max(perf[d], .001)**2 for d in dms_names], dtype=float)
    weights /= weights.sum()

    for dms_name, clinical_name in pairs:
        fit_col = 'fit__' + dms_name
        clinical[fit_col] = clinical.groupby('protein_file')[clinical_name].transform(rank01)
        if signs[dms_name] < 0:
            clinical[fit_col] = 1 - clinical[fit_col]
        clinical['risk__' + dms_name] = 1 - clinical[fit_col]
    fit_cols = ['fit__' + d for d in dms_names]
    X = clinical[fit_cols].to_numpy(float)
    ok = np.isfinite(X); W = ok * weights.reshape(1,-1); denom = W.sum(axis=1)
    fit = np.full(len(clinical), np.nan); valid = denom > 0
    fit[valid] = np.nansum(X[valid]*weights.reshape(1,-1), axis=1)/denom[valid]
    clinical['ensemble_risk'] = 1-fit

    pred_cols = ['ensemble_risk'] + ['risk__'+d for d in dms_names]
    tab = auc_table(clinical, pred_cols)
    strict = tab[tab['strict_external'].eq(1)].copy()
    summaries = []
    for c in pred_cols:
        summaries.append({'model': c, 'all_mean_auc': float(tab[c].mean()), 'strict_mean_auc': float(strict[c].mean()), 'n_all': int(tab[c].notna().sum()), 'n_strict': int(strict[c].notna().sum())})
    summary = pd.DataFrame(summaries).sort_values('strict_mean_auc', ascending=False)
    frozen_best = max(dms_names, key=lambda d: np.nan_to_num(perf[d], nan=-9))
    frozen_col = 'risk__' + frozen_best
    posthoc_col = str(summary[summary['model']!='ensemble_risk'].iloc[0]['model'])

    result = {
        'question': 'Does a human-DMS-trained frozen rank ensemble improve clinical pathogenicity discrimination?',
        'common_predictor_pairs': pairs,
        'dms_weights': {d: float(w) for d,w in zip(dms_names,weights)},
        'frozen_best_individual_from_dms': frozen_best,
        'posthoc_best_clinical_individual': posthoc_col,
        'all_clinical': {
            'proteins': int(len(tab)), 'variants': int(len(clinical)),
            'vs_frozen_best': compare(tab,'ensemble_risk',frozen_col),
            'vs_posthoc_best': compare(tab,'ensemble_risk',posthoc_col),
        },
        'strict_no_exact_sequence_overlap': {
            'proteins': int(len(strict)), 'variants': int(clinical.loc[clinical['strict_external'].eq(1)].shape[0]),
            'vs_frozen_best': compare(strict,'ensemble_risk',frozen_col),
            'vs_posthoc_best': compare(strict,'ensemble_risk',posthoc_col),
        },
        'clinical_diagnostics': clinical_diag,
        'dms_diagnostics': dms_diag,
        'clinical_transfer_advance_confirmed': False,
        'open_problem_fully_solved': False,
    }
    ci = result['strict_no_exact_sequence_overlap']['vs_posthoc_best']['gain_95ci']
    result['clinical_transfer_advance_confirmed'] = bool(ci[0] is not None and ci[0] > 0)

    tab.to_csv(out/'clinical_protein_auc.csv', index=False)
    summary.to_csv(out/'clinical_model_summary.csv', index=False)
    clinical[['protein_file','mutant','label','strict_external','ensemble_risk'] + ['risk__'+d for d in dms_names]].to_csv(out/'clinical_predictions.csv.gz', index=False, compression='gzip')
    (out/'result.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    s = result['strict_no_exact_sequence_overlap']['vs_posthoc_best']
    (out/'REPORT.md').write_text(f'''# ProteinGym DMS-to-clinical transfer v6\n\nA rank ensemble was frozen from human DMS train+tuning assays and transferred to clinical pathogenic/benign variants. Clinical proteins with an exact sequence match to any DMS training protein were excluded from the strict analysis.\n\n- Common predictors: {', '.join(dms_names)}\n- Strict clinical proteins: {result['strict_no_exact_sequence_overlap']['proteins']}\n- Strict clinical variants: {result['strict_no_exact_sequence_overlap']['variants']:,}\n- Ensemble mean protein AUC: {s['ensemble_mean_auc']:.4f}\n- Post-hoc best common individual mean AUC: {s['comparator_mean_auc']:.4f}\n- Paired gain: {s['mean_gain']:+.4f}\n- 95% CI: {s['gain_95ci']}\n- Advance confirmed: {result['clinical_transfer_advance_confirmed']}\n\nThis is a transfer benchmark, not a complete solution to clinical variant interpretation.\n''', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)

if __name__ == '__main__':
    main()
