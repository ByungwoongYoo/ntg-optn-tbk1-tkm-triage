#!/usr/bin/env python3
"""Compare the frozen DMS-trained ensemble with every available ProteinGym clinical baseline.

The comparison uses the exact 700-gene conservative unseen-protein set from the
prelabel-frozen transfer experiment. For each baseline, AUC is computed on the
intersection of variants with the frozen evaluated set. Models must cover at least
90% of the 700 genes and 90% of frozen variants to enter the primary zero-shot
comparison. Clinically supervised predictors are reported separately and are not
used to judge the DMS-transfer hypothesis.
"""
from __future__ import annotations

import argparse
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SEED = 20260817


def member_map(archive: zipfile.ZipFile):
    return {Path(name).name: name for name in archive.namelist() if name.lower().endswith('.csv')}


def direction_configs(config: dict):
    out = {}
    for group_name in [
        'model_list_zero_shot_substitutions_clinical',
        'model_list_supervised_substitutions_clinical',
    ]:
        category = 'zero_shot' if 'zero_shot' in group_name else 'clinically_supervised'
        for model, spec in config.get(group_name, {}).items():
            out[model] = {
                'direction': float(spec.get('directionality', 1)),
                'category': category,
            }
    return out


def rank01(series):
    return pd.to_numeric(series, errors='coerce').rank(method='average', pct=True)


def gene_auc(frame, score_col):
    rows = []
    for gene, group in frame.groupby('DMS_id', sort=False):
        valid = group[['benign_label', score_col]].dropna()
        n_benign = int((valid['benign_label'] == 1).sum())
        n_pathogenic = int((valid['benign_label'] == 0).sum())
        if min(n_benign, n_pathogenic) < 5:
            continue
        rows.append({
            'DMS_id': gene,
            'n_variants': len(valid),
            'n_benign': n_benign,
            'n_pathogenic': n_pathogenic,
            'auc': float(roc_auc_score(valid['benign_label'], valid[score_col])),
        })
    return pd.DataFrame(rows)


def bootstrap(values, n=50000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(SEED)
    draws = np.empty(n)
    for i in range(n):
        draws[i] = np.mean(rng.choice(values, size=len(values), replace=True))
    return [float(np.quantile(draws, .025)), float(np.quantile(draws, .975))]


def sign_flip(values, n=200000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    obs = float(np.mean(values))
    rng = np.random.default_rng(SEED)
    count = 0
    for _ in range(n):
        perm = float(np.mean(values * rng.choice([-1.0, 1.0], size=len(values))))
        count += int(perm >= obs)
    return float((count + 1) / (n + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--clinical-scores', required=True)
    ap.add_argument('--clinical-reference', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--frozen-variants', required=True)
    ap.add_argument('--corrected-overlap', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    config = json.loads(Path(args.config).read_text())
    specs = direction_configs(config)
    clinical_ref = pd.read_csv(args.clinical_reference)
    ref_member = dict(zip(clinical_ref['DMS_id'].astype(str), clinical_ref['DMS_filename'].astype(str)))
    frozen = pd.read_csv(args.frozen_variants, compression='gzip')
    overlap = pd.read_csv(args.corrected_overlap)
    unseen = set(overlap.loc[~overlap['seen_conservative_union'], 'DMS_id'].astype(str))
    frozen = frozen[frozen['DMS_id'].astype(str).isin(unseen)].copy()
    frozen['DMS_id'] = frozen['DMS_id'].astype(str)
    frozen['mutant'] = frozen['mutant'].astype(str)
    ensemble_auc = gene_auc(frozen, 'prediction_ridge100').rename(columns={'auc':'ensemble_auc'})
    if len(ensemble_auc) != 700:
        raise RuntimeError(f'Expected 700 ensemble genes, got {len(ensemble_auc)}')

    summaries = []
    pairwise_rows = []
    model_failures = []
    with zipfile.ZipFile(args.clinical_scores) as archive:
        members = member_map(archive)
        for model, spec in specs.items():
            gene_parts = []
            variants_scored = 0
            genes_with_file = 0
            for gene in sorted(unseen):
                filename = ref_member.get(gene, f'{gene}.csv')
                member = members.get(filename) or members.get(f'{gene}.csv')
                if member is None:
                    continue
                target = frozen[frozen['DMS_id'] == gene][['DMS_id','mutant','benign_label','prediction_ridge100']]
                if target.empty:
                    continue
                try:
                    with archive.open(member) as handle:
                        header = pd.read_csv(handle, nrows=0).columns.tolist()
                    if model not in header or 'mutant' not in header:
                        continue
                    with archive.open(member) as handle:
                        scores = pd.read_csv(handle, usecols=['mutant', model], low_memory=False)
                    scores['mutant'] = scores['mutant'].astype(str)
                    scores[model] = pd.to_numeric(scores[model], errors='coerce') * spec['direction']
                    scores = scores.drop_duplicates('mutant')
                    merged = target.merge(scores, on='mutant', how='left')
                    merged['baseline_score'] = rank01(merged[model])
                    gene_parts.append(merged)
                    variants_scored += int(merged['baseline_score'].notna().sum())
                    genes_with_file += 1
                except Exception as exc:
                    model_failures.append({'model':model,'DMS_id':gene,'error':repr(exc)})
            if not gene_parts:
                continue
            model_frame = pd.concat(gene_parts, ignore_index=True)
            baseline_auc = gene_auc(model_frame, 'baseline_score').rename(columns={'auc':'baseline_auc'})
            paired = ensemble_auc[['DMS_id','ensemble_auc']].merge(
                baseline_auc[['DMS_id','baseline_auc']], on='DMS_id', how='inner'
            )
            paired['difference'] = paired['ensemble_auc'] - paired['baseline_auc']
            gene_coverage = len(paired) / 700
            variant_coverage = variants_scored / len(frozen)
            summaries.append({
                'model': model,
                'category': spec['category'],
                'n_paired_genes': len(paired),
                'gene_coverage': gene_coverage,
                'variant_coverage': variant_coverage,
                'ensemble_mean_auc_on_common_genes': float(paired['ensemble_auc'].mean()),
                'baseline_mean_auc': float(paired['baseline_auc'].mean()),
                'ensemble_minus_baseline': float(paired['difference'].mean()),
                'median_difference': float(paired['difference'].median()),
                'fraction_genes_ensemble_better': float((paired['difference'] > 0).mean()),
            })
            for _, row in paired.iterrows():
                pairwise_rows.append({'model':model,'category':spec['category'],**row.to_dict()})
            print(model, spec['category'], len(paired), variant_coverage, flush=True)

    summary = pd.DataFrame(summaries)
    pairwise = pd.DataFrame(pairwise_rows)
    eligible_zero = summary[
        (summary['category']=='zero_shot') &
        (summary['gene_coverage']>=0.90) &
        (summary['variant_coverage']>=0.90)
    ].copy()
    if eligible_zero.empty:
        raise RuntimeError('No zero-shot baseline met 90% coverage gate')
    strongest = eligible_zero.sort_values(
        ['baseline_mean_auc','gene_coverage'], ascending=False
    ).iloc[0]
    strongest_name = strongest['model']
    p = pairwise[pairwise['model']==strongest_name].copy()
    differences = p['difference'].to_numpy(float)
    comparison = {
        'strongest_eligible_zero_shot_baseline': strongest_name,
        'n_common_genes': int(len(p)),
        'ensemble_mean_auc': float(p['ensemble_auc'].mean()),
        'baseline_mean_auc': float(p['baseline_auc'].mean()),
        'mean_difference': float(differences.mean()),
        'paired_bootstrap_95ci': bootstrap(differences),
        'one_sided_sign_flip_p': sign_flip(differences),
        'fraction_genes_ensemble_better': float((differences>0).mean()),
    }
    comparison['decision'] = (
        'ENSEMBLE_BEATS_STRONGEST_AVAILABLE_ZERO_SHOT_BASELINE'
        if comparison['paired_bootstrap_95ci'][0] > 0 and comparison['one_sided_sign_flip_p'] < 0.05
        else 'NO_ROBUST_ADVANCE_OVER_STRONGEST_AVAILABLE_ZERO_SHOT_BASELINE'
    )
    result = {
        'scope':'700 conservative DMS-unseen clinical genes from frozen transfer experiment',
        'frozen_variant_rows':int(len(frozen)),
        'models_evaluated':int(len(summary)),
        'eligible_zero_shot_models_90pct_coverage':eligible_zero['model'].tolist(),
        'primary_comparison':comparison,
        'important_boundary':(
            'Clinically supervised predictors are reported only descriptively because they use clinical labels or related features. '
            'This analysis is post-label benchmarking of fixed predictions; it does not change the prediction freeze.'
        ),
        'failures_n':len(model_failures),
    }
    summary.sort_values(['category','baseline_mean_auc'],ascending=[True,False]).to_csv(out/'all_baseline_summary.csv',index=False)
    pairwise.to_csv(out/'all_baseline_gene_pairwise.csv.gz',index=False,compression='gzip')
    pd.DataFrame(model_failures).to_csv(out/'baseline_failures.csv',index=False)
    (out/'result.json').write_text(json.dumps(result,indent=2))
    c = comparison
    (out/'REPORT.md').write_text(f'''# Frozen clinical ensemble versus all available ProteinGym baselines

## Scope

The exact 700-gene conservative DMS-unseen set and frozen ensemble predictions are unchanged.
All available clinical score columns were evaluated on the same frozen variant set. A zero-shot
baseline had to cover at least 90% of genes and variants to enter the primary comparison.

## Primary result

- Strongest eligible zero-shot baseline: **{c['strongest_eligible_zero_shot_baseline']}**
- Common genes: **{c['n_common_genes']}**
- Frozen DMS-trained ensemble mean AUC: **{c['ensemble_mean_auc']:.4f}**
- Strongest baseline mean AUC: **{c['baseline_mean_auc']:.4f}**
- Difference: **{c['mean_difference']:+.4f}**
- Paired bootstrap 95% CI: **[{c['paired_bootstrap_95ci'][0]:+.4f}, {c['paired_bootstrap_95ci'][1]:+.4f}]**
- One-sided sign-flip p: **{c['one_sided_sign_flip_p']:.6f}**
- Genes improved: **{100*c['fraction_genes_ensemble_better']:.1f}%**
- Decision: **`{c['decision']}`**

Clinically supervised predictors are secondary descriptive comparators, not fair zero-shot baselines.
''')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    main()
