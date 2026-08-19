#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def args():
 p=argparse.ArgumentParser();p.add_argument('--paf',required=True);p.add_argument('--truth-mapping',required=True);p.add_argument('--min-identity',type=float,default=.90);p.add_argument('--min-alignment',type=int,default=500);p.add_argument('--out',required=True);return p.parse_args()

def union_len(intervals):
 if not intervals:return 0
 xs=sorted(intervals);s,e=xs[0];n=0
 for a,b in xs[1:]:
  if a<=e:e=max(e,b)
  else:n+=e-s;s,e=a,b
 return n+e-s

def main():
 a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);truth={};glen=defaultdict(int)
 with open(a.truth_mapping,newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):truth[r['sequence_id']]=r;glen[r['genome_id']]+=int(r['length'])
 iv=defaultdict(list);best_identity=defaultdict(float);queries=defaultdict(set)
 with open(a.paf) as f:
  for line in f:
   x=line.rstrip().split('\t')
   if len(x)<12:continue
   q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12];nm=int(nm);al=int(al);ts=int(ts);te=int(te)
   ident=nm/al if al else 0
   if t not in truth or ident<a.min_identity or al<a.min_alignment:continue
   iv[t].append((ts,te));best_identity[t]=max(best_identity[t],ident);queries[t].add(q)
 cont=[];gcov=defaultdict(int);gqueries=defaultdict(set)
 for sid,r in truth.items():
  cov=union_len(iv[sid]);gcov[r['genome_id']]+=cov;gqueries[r['genome_id']]|=queries[sid];cont.append({'sequence_id':sid,'genome_id':r['genome_id'],'length':r['length'],'covered_bp':cov,'covered_fraction':cov/int(r['length']),'best_identity':best_identity[sid],'assembly_contigs':len(queries[sid])})
 genomes=[]
 for g in sorted(glen):genomes.append({'genome_id':g,'length':glen[g],'covered_bp':gcov[g],'recovery_fraction':gcov[g]/glen[g] if glen[g] else 0,'assembly_contigs':len(gqueries[g]),'recovered_50':str(gcov[g]>=.5*glen[g]).lower(),'recovered_90':str(gcov[g]>=.9*glen[g]).lower()})
 for name,rows in [('per_truth_contig.tsv',cont),('per_genome_recovery.tsv',genomes)]:
  with open(out/name,'w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 vals=[r['recovery_fraction'] for r in genomes];summary={'n_genomes':len(genomes),'mean_genome_recovery':sum(vals)/len(vals) if vals else 0,'median_genome_recovery':sorted(vals)[len(vals)//2] if vals else 0,'genomes_recovered_50':sum(r['recovered_50']=='true' for r in genomes),'genomes_recovered_90':sum(r['recovered_90']=='true' for r in genomes),'min_identity':a.min_identity,'truth_only_evaluation':True}
 (out/'GOLD_COVERAGE_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
