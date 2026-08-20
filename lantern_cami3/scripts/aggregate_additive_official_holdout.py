#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, math, statistics
from pathlib import Path

EXPECTED = [
    'official_02_07','official_03_05','official_04_08','official_06_19',
    'official_09_12','official_10_15','official_11_13','official_14_16','official_17_18'
]
TIERS = ['exact_target_withheld','species_neighbour_withheld','genus_neighbour_withheld','deep_neighbour_withheld']

def args():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);return p.parse_args()

def mean(vals):
    vals=[float(v) for v in vals if v is not None and math.isfinite(float(v))]
    return statistics.fmean(vals) if vals else None

def median(vals):
    vals=[float(v) for v in vals if v is not None and math.isfinite(float(v))]
    return statistics.median(vals) if vals else None

def main():
    a=args(); root=Path(a.root); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    files=sorted(root.rglob('PAIR_DECISION.json'))
    rows=[]
    for f in files:
        x=json.loads(f.read_text()); x['_source']=str(f); rows.append(x)
    by={r['pair_id']:r for r in rows}
    missing=sorted(set(EXPECTED)-set(by)); extra=sorted(set(by)-set(EXPECTED))
    ordered=[by[p] for p in EXPECTED if p in by]
    complete=len(ordered)==len(EXPECTED) and not missing and not extra

    gf=[r['genome_fraction_gain_vs_strongest_pp'] for r in ordered]
    mr=[r['mean_recovery_gain_vs_strongest_pp'] for r in ordered]
    ab=[r['longitudinal_ablation_drop_pp'] for r in ordered]
    chim=[r['relative_chimera_change_vs_strongest'] for r in ordered]
    align=[r['alignment_ratio_change_vs_strongest'] for r in ordered]
    pseudo={t:[r['pseudo_novel_gain_pp'].get(t) for r in ordered] for t in TIERS}

    directional={
        'all_pairs_complete': complete,
        'mean_genome_fraction_gain_positive': (mean(gf) or 0)>0,
        'mean_recovery_gain_positive': (mean(mr) or 0)>0,
        'at_least_six_positive_pairs': sum(float(v)>0 for v in gf)>=6,
        'mean_longitudinal_ablation_positive': (mean(ab) or 0)>0,
        'maximum_relative_chimera_increase_le_0_10': bool(chim) and max(float(v) for v in chim)<=0.10,
        'maximum_alignment_ratio_increase_le_0_10': bool(align) and max(float(v) for v in align)<=0.10,
        'pseudo_novel_mean_positive_all_tiers': all((mean(pseudo[t]) or 0)>0 for t in TIERS),
        'pseudo_novel_positive_at_least_six_pairs_all_tiers': all(sum(v is not None and float(v)>0 for v in pseudo[t])>=6 for t in TIERS),
    }
    strict={
        'all_pairs_complete': complete,
        'mean_genome_fraction_gain_ge_0_5_pp': (mean(gf) or 0)>=0.5,
        'mean_recovery_gain_ge_1_0_pp': (mean(mr) or 0)>=1.0,
        'mean_longitudinal_ablation_ge_0_5_pp': (mean(ab) or 0)>=0.5,
        'maximum_relative_chimera_increase_le_0_10': bool(chim) and max(float(v) for v in chim)<=0.10,
        'maximum_alignment_ratio_increase_le_0_10': bool(align) and max(float(v) for v in align)<=0.10,
        'pseudo_novel_mean_positive_all_tiers': all((mean(pseudo[t]) or 0)>0 for t in TIERS),
    }
    summary={
        'status':'COMPLETE' if complete else 'INCOMPLETE',
        'verdict':'STRICT_PASS' if all(strict.values()) else ('DIRECTIONAL_REPLICATION' if all(directional.values()) else 'FAIL'),
        'n_expected_pairs':len(EXPECTED),'n_completed_pairs':len(ordered),'missing_pairs':missing,'extra_pairs':extra,
        'mean_genome_fraction_gain_vs_strongest_pp':mean(gf),
        'median_genome_fraction_gain_vs_strongest_pp':median(gf),
        'positive_genome_fraction_pairs':sum(float(v)>0 for v in gf),
        'mean_mean_recovery_gain_vs_strongest_pp':mean(mr),
        'median_mean_recovery_gain_vs_strongest_pp':median(mr),
        'positive_mean_recovery_pairs':sum(float(v)>0 for v in mr),
        'mean_longitudinal_ablation_drop_pp':mean(ab),
        'positive_longitudinal_contribution_pairs':sum(float(v)>0 for v in ab),
        'maximum_relative_chimera_change_vs_strongest':max(chim) if chim else None,
        'maximum_alignment_ratio_change_vs_strongest':max(align) if align else None,
        'pseudo_novel':{t:{'mean_gain_pp':mean(pseudo[t]),'median_gain_pp':median(pseudo[t]),'positive_pairs':sum(v is not None and float(v)>0 for v in pseudo[t]),'informative_pairs':sum(v is not None for v in pseudo[t])} for t in TIERS},
        'directional_gate':directional,'strict_gate':strict,'per_pair':ordered,
        'boundary':('The frozen additive rule was selected on public Toy samples 0/1, then reapplied without pair-specific tuning to pre-truth inputs preserved from all nine officially corrected public Toy holdout pairs at 1% depth. '
                    'These Toy pairs had been evaluated previously by the broader LANTERN project, so this is a broad frozen-rule reapplication audit, not a prospectively untouched or hidden CAMI challenge validation. No CAMI rank or novel-genome discovery is established.')
    }
    (out/'ADDITIVE_OFFICIAL_HOLDOUT_AGGREGATE.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    flat=[]
    for r in ordered:
        flat.append({k:v for k,v in r.items() if not isinstance(v,(dict,list)) and not k.startswith('_')})
    keys=sorted({k for r in flat for k in r})
    with open(out/'ADDITIVE_OFFICIAL_HOLDOUT_RESULTS.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(flat)
    lines=['# Frozen additive-rescue reapplication across corrected CAMI III Toy pairs','',f"- Verdict: **{summary['verdict']}**",f"- Complete pairs: **{len(ordered)}/{len(EXPECTED)}**",f"- Mean genome-fraction gain vs strongest baseline: **{summary['mean_genome_fraction_gain_vs_strongest_pp']} pp**",f"- Positive genome-fraction pairs: **{summary['positive_genome_fraction_pairs']}/{len(EXPECTED)}**",f"- Mean per-genome recovery gain: **{summary['mean_mean_recovery_gain_vs_strongest_pp']} pp**",f"- Mean longitudinal-ablation drop: **{summary['mean_longitudinal_ablation_drop_pp']} pp**",f"- Maximum relative chimera change: **{summary['maximum_relative_chimera_change_vs_strongest']}**",'',summary['boundary']]
    (out/'ADDITIVE_OFFICIAL_HOLDOUT_AGGREGATE.md').write_text('\n'.join(lines)+'\n')
    hashes=[]
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name!='SHA256SUMS.txt': hashes.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='per_pair'},indent=2,sort_keys=True))
    if not complete: raise SystemExit(1)
if __name__=='__main__': main()
