#!/usr/bin/env python3
"""Summarize truth-side CAMI Toy assembly evaluations under a frozen decision gate.

The script does not select or tune contigs. It opens only already-frozen assemblies and
gold-side evaluation tables, ranks predeclared baselines by the primary genome-fraction
metric, calculates paired per-genome differences, and applies the prespecified gate.
"""
from __future__ import annotations
import argparse,csv,json,math,random,statistics
from collections import defaultdict
from pathlib import Path


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--method',action='append',required=True,help='METHOD=EVALUATION_DIR')
    p.add_argument('--manifest',required=True,help='TSV with method,role,mode,scope,ablation')
    p.add_argument('--decision-config',required=True)
    p.add_argument('--full-method',default='lantern_full')
    p.add_argument('--no-longitudinal-method',default='lantern_no_longitudinal')
    p.add_argument('--pseudo-splits',default=None,help='optional TSV genome_id,tier')
    p.add_argument('--seed',type=int,default=20260819)
    p.add_argument('--bootstrap',type=int,default=5000)
    p.add_argument('--out',required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_csv(path,rows):
    path=Path(path);fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def domain(g):
    if g.startswith('hASV'):return 'host'
    if g.startswith('hvASV') or g.startswith('vASV'):return 'virus'
    if g.startswith('pASV'):return 'plasmid'
    if g.startswith('fASV'):return 'fungus'
    return 'bacteria_archaea'

def percentile(x,p):
    if not x:return None
    y=sorted(x);i=(len(y)-1)*p;lo=math.floor(i);hi=math.ceil(i)
    return y[lo] if lo==hi else y[lo]*(hi-i)+y[hi]*(i-lo)

def bootstrap_mean_ci(diffs,n,seed):
    if not diffs:return (None,None,None)
    rng=random.Random(seed);m=len(diffs);vals=[]
    for _ in range(n):vals.append(sum(diffs[rng.randrange(m)] for _ in range(m))/m)
    return (sum(diffs)/m,percentile(vals,.025),percentile(vals,.975))

def main():
    a=parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    cfg=json.load(open(a.decision_config));manifest={r['method']:r for r in read_tsv(a.manifest)};method_dirs={}
    for spec in a.method:
        name,path=spec.split('=',1)
        if name in method_dirs:raise ValueError(f'duplicate method {name}')
        if name not in manifest:raise ValueError(f'{name} absent from manifest')
        method_dirs[name]=Path(path)
    summaries={};per_genome={}
    for name,root in method_dirs.items():
        summary_file=root/'GOLD_COVERAGE_SUMMARY.json';genome_file=root/'per_genome_recovery.tsv'
        if not summary_file.exists() or not genome_file.exists():raise FileNotFoundError(f'{name}: missing evaluation output')
        s=json.load(open(summary_file));s.update(manifest[name]);s['method']=name;summaries[name]=s
        per_genome[name]={r['genome_id']:float(r['recovery_fraction']) for r in read_tsv(genome_file)}
    baselines=[s for s in summaries.values() if s.get('role')=='baseline' and str(s.get('eligible_primary_baseline','true')).lower()=='true']
    if not baselines:raise SystemExit('No eligible baseline in manifest')
    strongest=max(baselines,key=lambda s:(float(s['genome_fraction_percent']),float(s['mean_genome_recovery']),-float(s.get('cross_binid_chimeric_bp_fraction') or 0),s['method']))
    full=summaries.get(a.full_method)
    if full is None:raise SystemExit(f'Missing full method {a.full_method}')
    ab=summaries.get(a.no_longitudinal_method)
    common=sorted(set(per_genome[a.full_method])&set(per_genome[strongest['method']]))
    diffs=[100*(per_genome[a.full_method][g]-per_genome[strongest['method']][g]) for g in common]
    paired_mean,ci_lo,ci_hi=bootstrap_mean_ci(diffs,a.bootstrap,a.seed)
    gf_gain=float(full['genome_fraction_percent'])-float(strongest['genome_fraction_percent'])
    mean_gain=100*(float(full['mean_genome_recovery'])-float(strongest['mean_genome_recovery']))
    full_ch=float(full.get('cross_binid_chimeric_bp_fraction') or 0);base_ch=float(strongest.get('cross_binid_chimeric_bp_fraction') or 0)
    relative_chimera=(full_ch-base_ch)/base_ch if base_ch>0 else (0 if full_ch==0 else None)
    ablation_drop=float(full['genome_fraction_percent'])-float(ab['genome_fraction_percent']) if ab is not None else None
    pseudo={};pseudo_positive=False
    if a.pseudo_splits:
        tiers=defaultdict(list)
        for r in read_tsv(a.pseudo_splits):tiers[r['tier']].append(r['genome_id'])
        for tier,gs in tiers.items():
            ds=[100*(per_genome[a.full_method][g]-per_genome[strongest['method']][g]) for g in gs if g in per_genome[a.full_method] and g in per_genome[strongest['method']]]
            pseudo[tier]={'n':len(ds),'mean_recovery_gain_percentage_points':sum(ds)/len(ds) if ds else None,'median_gain_percentage_points':statistics.median(ds) if ds else None,'positive_genomes':sum(x>0 for x in ds),'negative_genomes':sum(x<0 for x in ds)}
            if ds and sum(ds)/len(ds)>0:pseudo_positive=True
    gate_gf=gf_gain>=float(cfg['minimum_genome_fraction_gain_percentage_points'])
    gate_mean=mean_gain>=float(cfg['minimum_mean_genome_recovery_gain_percentage_points'])
    gate_chimera=((relative_chimera is not None and relative_chimera<=float(cfg['maximum_relative_misassembly_increase'])) or (gf_gain>=float(cfg['acceptable_misassembly_tradeoff_if_genome_fraction_gain'])))
    gate_ablation=ablation_drop is not None and ablation_drop>=float(cfg['minimum_longitudinal_ablation_drop_percentage_points'])
    gate_pseudo=bool(pseudo) and pseudo_positive
    gates={'genome_fraction_gain':gate_gf,'mean_genome_recovery_gain':gate_mean,'chimera_tradeoff':gate_chimera,'longitudinal_ablation_contribution':gate_ablation,'pseudo_novel_or_isolation_gain':gate_pseudo}
    status='SUCCESS' if all(gates.values()) else ('PARTIAL' if gf_gain>0 or mean_gain>0 or pseudo_positive else 'NEGATIVE')
    decision={'status':status,'full_method':a.full_method,'strongest_baseline':strongest['method'],'genome_fraction_gain_percentage_points':gf_gain,'mean_genome_recovery_gain_percentage_points':mean_gain,'paired_per_genome_mean_gain_percentage_points':paired_mean,'paired_bootstrap_95_ci':[ci_lo,ci_hi],'n_paired_genomes':len(common),'full_chimera_bp_fraction':full_ch,'baseline_chimera_bp_fraction':base_ch,'relative_chimera_change':relative_chimera,'longitudinal_ablation_drop_percentage_points':ablation_drop,'gates':gates,'pseudo_novel_or_isolation_tiers':pseudo,'decision_config':cfg,'claim_boundary':'Toy-dataset performance is development evidence only. It is not a CAMI III challenge result, rank, or novel-organism discovery.'}
    (out/'FINAL_DECISION.json').write_text(json.dumps(decision,indent=2,allow_nan=False)+'\n')
    rows=[]
    for name,s in sorted(summaries.items(),key=lambda kv:(kv[1].get('role')!='baseline',-float(kv[1]['genome_fraction_percent']),kv[0])):
        rows.append({'method':name,'role':s.get('role',''),'mode':s.get('mode',''),'scope':s.get('scope',''),'ablation':s.get('ablation',''),'genome_fraction_percent':s['genome_fraction_percent'],'mean_genome_recovery_percent':100*float(s['mean_genome_recovery']),'median_genome_recovery_percent':100*float(s['median_genome_recovery']),'genomes_recovered_50':s['genomes_recovered_50'],'genomes_recovered_90':s['genomes_recovered_90'],'assembly_contigs_total':s['assembly_contigs_total'],'assembly_total_bp':s['assembly_total_bp'],'cross_binid_chimeric_contigs':s['cross_binid_chimeric_contigs'],'cross_binid_chimeric_bp_fraction':s['cross_binid_chimeric_bp_fraction'],'alignment_to_unique_truth_ratio':s['alignment_to_unique_truth_ratio'],'strongest_baseline':str(name==strongest['method']).lower()})
    write_csv(out/'METHOD_COMPARISON.csv',rows)
    domain_rows=[];genomes=sorted(set().union(*[set(v) for v in per_genome.values()]))
    for name in sorted(per_genome):
        by=defaultdict(list)
        for g in genomes:
            if g in per_genome[name]:by[domain(g)].append(per_genome[name][g])
        for d,vals in sorted(by.items()):domain_rows.append({'method':name,'domain':d,'n_genomes':len(vals),'mean_recovery_percent':100*sum(vals)/len(vals),'median_recovery_percent':100*statistics.median(vals),'recovered_50':sum(v>=.5 for v in vals),'recovered_90':sum(v>=.9 for v in vals)})
    write_csv(out/'DOMAIN_RECOVERY.csv',domain_rows)
    paired_rows=[{'genome_id':g,'domain':domain(g),'full_recovery_percent':100*per_genome[a.full_method][g],'baseline_recovery_percent':100*per_genome[strongest['method']][g],'difference_percentage_points':100*(per_genome[a.full_method][g]-per_genome[strongest['method']][g])} for g in common]
    write_csv(out/'PAIRED_PER_GENOME_DIFFERENCES.csv',paired_rows)
    lines=['# LANTERN frozen decision-gate result','',f"- Status: **{status}**",f"- Full method: `{a.full_method}`",f"- Strongest prespecified baseline: `{strongest['method']}`",f"- Genome-fraction gain: **{gf_gain:.4f} percentage points**",f"- Mean per-genome recovery gain: **{mean_gain:.4f} percentage points**",f"- Paired bootstrap mean gain (95% CI): **{paired_mean:.4f} ({ci_lo:.4f}, {ci_hi:.4f}) percentage points**" if paired_mean is not None else '- Paired bootstrap: not evaluable',f"- Longitudinal ablation drop: **{ablation_drop:.4f} percentage points**" if ablation_drop is not None else '- Longitudinal ablation: not evaluable',f"- Relative conservative chimera-proxy change: **{relative_chimera}**",'', '## Gates','']
    for k,v in gates.items():lines.append(f"- {k}: **{'PASS' if v else 'FAIL'}**")
    lines+=['','## Claim boundary','',decision['claim_boundary'],'']
    (out/'FINAL_DECISION.md').write_text('\n'.join(lines)+'\n');print(json.dumps(decision,indent=2,allow_nan=False))
if __name__=='__main__':main()
