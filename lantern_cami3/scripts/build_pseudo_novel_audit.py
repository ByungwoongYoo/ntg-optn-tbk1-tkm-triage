#!/usr/bin/env python3
"""Build evaluation-only pseudo-novel tiers from truth genomes.

The LANTERN assembly path is reference-free and already frozen before this script runs.
This post hoc *evaluation* withholds each target and defines progressively harder target
sets by the nearest remaining Toy truth-genome Mash distance. It never feeds truth or
Mash distances back into contig selection.
"""
from __future__ import annotations
import argparse,csv,json,statistics
from collections import defaultdict
from pathlib import Path

def args():
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--mash',required=True);p.add_argument('--recovery',action='append',required=True,help='METHOD=per_genome_recovery.tsv');p.add_argument('--min-truth-bp',type=int,default=50000);p.add_argument('--species-distance',type=float,default=.05);p.add_argument('--genus-distance',type=float,default=.15);p.add_argument('--deep-distance',type=float,default=.25);p.add_argument('--out',required=True);return p.parse_args()

def norm(x):
    p=Path(str(x));return str(p),p.name,p.stem

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    manifest=list(csv.DictReader(open(a.manifest),delimiter='\t'));eligible={r['genome_id']:r for r in manifest if int(r['truth_bp'])>=a.min_truth_bp}
    alias={}
    for gid,r in eligible.items():
        for key in norm(r['reference_path']):alias[key]=gid
    dist=defaultdict(list);unresolved=set()
    with open(a.mash) as f:
        for line in f:
            x=line.rstrip().split('\t')
            if len(x)<3:continue
            ra,qa,d=x[:3];g1=next((alias[k] for k in norm(ra) if k in alias),None);g2=next((alias[k] for k in norm(qa) if k in alias),None)
            if g1 is None:unresolved.add(ra)
            if g2 is None:unresolved.add(qa)
            if g1 is None or g2 is None or g1==g2:continue
            try:dv=float(d)
            except ValueError:continue
            dist[g1].append((dv,g2))
    rows=[]
    for gid,r in sorted(eligible.items()):
        nearest=min(dist.get(gid,[(1.0,'NONE')]))
        nd,ng=nearest
        rows.append({'genome_id':gid,'truth_bp':r['truth_bp'],'samples':r['samples'],'nearest_remaining_genome':ng,'nearest_mash_distance':nd,'exact_target_withheld':'true','species_neighbour_withheld':str(nd>=a.species_distance).lower(),'genus_neighbour_withheld':str(nd>=a.genus_distance).lower(),'deep_neighbour_withheld':str(nd>=a.deep_distance).lower()})
    with open(out/'PSEUDO_NOVEL_TARGETS.tsv','w',newline='') as f:
        fields=list(rows[0]) if rows else ['genome_id','truth_bp','samples','nearest_remaining_genome','nearest_mash_distance','exact_target_withheld','species_neighbour_withheld','genus_neighbour_withheld','deep_neighbour_withheld'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    target={r['genome_id']:r for r in rows};detail=[];summary=[]
    tiers=['exact_target_withheld','species_neighbour_withheld','genus_neighbour_withheld','deep_neighbour_withheld']
    for spec in a.recovery:
        method,path=spec.split('=',1);rec={r['genome_id']:r for r in csv.DictReader(open(path),delimiter='\t')}
        for gid,t in target.items():
            if gid not in rec:continue
            detail.append({'method':method,**t,'recovery_fraction':float(rec[gid]['recovery_fraction']),'recovered_50':rec[gid]['recovered_50'],'recovered_90':rec[gid]['recovered_90']})
        for tier in tiers:
            vals=[r for r in detail if r['method']==method and r[tier]=='true']
            x=[r['recovery_fraction'] for r in vals]
            summary.append({'method':method,'tier':tier,'n_targets':len(vals),'mean_recovery':sum(x)/len(x) if x else None,'median_recovery':statistics.median(x) if x else None,'minimum_recovery':min(x) if x else None,'recovered_50':sum(r['recovered_50']=='true' for r in vals),'recovered_90':sum(r['recovered_90']=='true' for r in vals)})
    for name,data in [('PSEUDO_NOVEL_RECOVERY_DETAIL.tsv',detail),('PSEUDO_NOVEL_RECOVERY_SUMMARY.tsv',summary)]:
        fields=list(data[0]) if data else ['method','tier','n_targets','mean_recovery'];
        with open(out/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
    report={'n_manifest_genomes':len(manifest),'n_eligible_targets':len(rows),'min_truth_bp':a.min_truth_bp,'thresholds':{'species':a.species_distance,'genus':a.genus_distance,'deep':a.deep_distance},'n_unresolved_mash_ids':len(unresolved),'unresolved_examples':sorted(unresolved)[:10],'boundary':'Pseudo-novel tiers are evaluation-only target-withholding stress tests within the Toy community. They are not claims that a genome is taxonomically novel in RefSeq, and no truth information entered LANTERN selection.'}
    (out/'PSEUDO_NOVEL_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
