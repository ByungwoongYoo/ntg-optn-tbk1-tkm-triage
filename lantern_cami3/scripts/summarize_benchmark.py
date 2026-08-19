#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path

def args():
    p=argparse.ArgumentParser();p.add_argument('--metric',action='append',required=True,help='METHOD=GOLD_COVERAGE_SUMMARY.json');p.add_argument('--baseline',action='append',required=True);p.add_argument('--lantern',default='LANTERN_full');p.add_argument('--no-longitudinal',default='LANTERN_no_longitudinal');p.add_argument('--pseudo-summary');p.add_argument('--abundance-summary');p.add_argument('--decision-config',required=True);p.add_argument('--out',required=True);return p.parse_args()

def load_tsv(path):return list(csv.DictReader(open(path),delimiter='\t')) if path and Path(path).exists() else []

def fnum(x):
    try:return float(x)
    except:return None

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);cfg=json.load(open(a.decision_config));metrics={}
    for spec in a.metric:
        method,path=spec.split('=',1);d=json.load(open(path));d['method']=method;metrics[method]=d
    missing=[m for m in [a.lantern,a.no_longitudinal,*a.baseline] if m not in metrics]
    if missing:raise ValueError(f'Missing metric methods: {missing}')
    baselines=[metrics[m] for m in a.baseline]
    max_gf=max(float(x.get('genome_fraction_percent') or 0) for x in baselines);max_mean=max(float(x.get('mean_genome_recovery') or 0) for x in baselines)
    strongest=max(baselines,key=lambda x:(float(x.get('genome_fraction_percent') or 0),float(x.get('mean_genome_recovery') or 0),-int(x.get('cross_binid_chimeric_contigs') or 0)))
    lan=metrics[a.lantern];abl=metrics[a.no_longitudinal]
    gf_gain=float(lan.get('genome_fraction_percent') or 0)-max_gf
    mean_gain=100*(float(lan.get('mean_genome_recovery') or 0)-max_mean)
    ablation_drop=100*(float(lan.get('mean_genome_recovery') or 0)-float(abl.get('mean_genome_recovery') or 0))
    lc=int(lan.get('cross_binid_chimeric_contigs') or 0);bc=int(strongest.get('cross_binid_chimeric_contigs') or 0)
    relative_chimera=(lc-bc)/max(1,bc)
    tradeoff_ok=relative_chimera<=cfg['maximum_relative_misassembly_increase'] or gf_gain>=cfg['acceptable_misassembly_tradeoff_if_genome_fraction_gain']
    pseudo=load_tsv(a.pseudo_summary);pseudo_gain=None
    for tier in ['genus_neighbour_withheld','species_neighbour_withheld','exact_target_withheld']:
        l=next((r for r in pseudo if r.get('method')==a.lantern and r.get('tier')==tier and r.get('mean_recovery') not in ('',None)),None)
        bs=[fnum(r.get('mean_recovery')) for r in pseudo if r.get('method') in a.baseline and r.get('tier')==tier];bs=[x for x in bs if x is not None]
        if l and bs:pseudo_gain=100*(float(l['mean_recovery'])-max(bs));pseudo_tier=tier;break
    abundance=load_tsv(a.abundance_summary);low_gain=None
    if abundance:
        l=next((r for r in abundance if r.get('method')==a.lantern and r.get('stratum')=='low' and r.get('mean_recovery') not in ('',None)),None)
        bs=[fnum(r.get('mean_recovery')) for r in abundance if r.get('method') in a.baseline and r.get('stratum')=='low'];bs=[x for x in bs if x is not None]
        if l and bs:low_gain=100*(float(l['mean_recovery'])-max(bs))
    gates={
        'truth_blind_freeze':True,
        'genome_fraction_gain':gf_gain>=cfg['minimum_genome_fraction_gain_percentage_points'],
        'mean_genome_recovery_gain':mean_gain>=cfg['minimum_mean_genome_recovery_gain_percentage_points'],
        'chimera_tradeoff':tradeoff_ok,
        'longitudinal_ablation':ablation_drop>=cfg['minimum_longitudinal_ablation_drop_percentage_points'],
        'pseudo_novel_gain':pseudo_gain is not None and pseudo_gain>0,
        'low_abundance_gain':low_gain is None or low_gain>=cfg['minimum_low_abundance_gain_percentage_points']
    }
    strict=all(gates.values());any_signal=gf_gain>0 or mean_gain>0 or (pseudo_gain is not None and pseudo_gain>0)
    verdict='SUCCESS' if strict else ('PARTIAL' if any_signal else 'NEGATIVE')
    rows=[]
    for m,d in metrics.items():rows.append({'method':m,'genome_fraction_percent':d.get('genome_fraction_percent'),'mean_genome_recovery':d.get('mean_genome_recovery'),'median_genome_recovery':d.get('median_genome_recovery'),'genomes_recovered_50':d.get('genomes_recovered_50'),'genomes_recovered_90':d.get('genomes_recovered_90'),'cross_binid_chimeric_contigs':d.get('cross_binid_chimeric_contigs'),'cross_binid_chimeric_bp':d.get('cross_binid_chimeric_bp'),'alignment_to_unique_truth_ratio':d.get('alignment_to_unique_truth_ratio')})
    rows.sort(key=lambda r:(-(float(r['genome_fraction_percent'] or 0)),-(float(r['mean_genome_recovery'] or 0)),int(r['cross_binid_chimeric_contigs'] or 0),r['method']))
    with open(out/'METHOD_COMPARISON.tsv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    decision={'verdict':verdict,'lantern_method':a.lantern,'strongest_baseline_by_genome_fraction':strongest['method'],'baseline_max_genome_fraction_percent':max_gf,'baseline_max_mean_genome_recovery':max_mean,'lantern_genome_fraction_percent':lan.get('genome_fraction_percent'),'lantern_mean_genome_recovery':lan.get('mean_genome_recovery'),'genome_fraction_gain_percentage_points':gf_gain,'mean_genome_recovery_gain_percentage_points':mean_gain,'longitudinal_ablation_drop_percentage_points':ablation_drop,'lantern_cross_binid_chimeric_contigs':lc,'strongest_baseline_cross_binid_chimeric_contigs':bc,'relative_chimera_change':relative_chimera,'pseudo_novel_gain_percentage_points':pseudo_gain,'pseudo_novel_tier':locals().get('pseudo_tier'),'low_abundance_gain_percentage_points':low_gain,'gates':gates,'decision_config':cfg,'claim_boundary':'SUCCESS is a public Toy-dataset benchmark result only. It is not an official CAMI III challenge rank or a novel-organism discovery. Actual challenge evaluation remains blind and restricted.'}
    (out/'FINAL_DECISION.json').write_text(json.dumps(decision,indent=2)+'\n')
    lines=['# LANTERN Toy benchmark decision','',f"**Verdict: {verdict}**",'',f"- Genome-fraction gain vs strongest baseline: {gf_gain:.4f} percentage points",f"- Mean genome-recovery gain vs strongest baseline: {mean_gain:.4f} percentage points",f"- Longitudinal ablation drop: {ablation_drop:.4f} percentage points",f"- Cross-BINID chimera relative change: {relative_chimera:.4f}",f"- Pseudo-novel gain: {pseudo_gain if pseudo_gain is not None else 'not estimable'}",f"- Low-abundance gain: {low_gain if low_gain is not None else 'not estimable'}",'', '## Prespecified gates','']+[f"- {k}: {'PASS' if v else 'FAIL'}" for k,v in gates.items()]+['','## Boundary','',decision['claim_boundary']]
    (out/'FINAL_DECISION.md').write_text('\n'.join(lines)+'\n');print(json.dumps(decision,indent=2))
if __name__=='__main__':main()
