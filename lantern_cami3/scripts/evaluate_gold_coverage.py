#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def args():
 p=argparse.ArgumentParser();p.add_argument('--paf',required=True);p.add_argument('--truth-mapping',required=True);p.add_argument('--min-identity',type=float,default=.90);p.add_argument('--min-alignment',type=int,default=500);p.add_argument('--chimera-min-bp',type=int,default=1000);p.add_argument('--chimera-min-query-fraction',type=float,default=.10);p.add_argument('--out',required=True);return p.parse_args()

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
 iv=defaultdict(list);best_identity=defaultdict(float);queries=defaultdict(set);query_hits=defaultdict(list);query_lengths={};total_accepted_alignment_bp=0
 with open(a.paf) as f:
  for line in f:
   x=line.rstrip().split('\t')
   if len(x)<12:continue
   q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12];ql=int(ql);qs=int(qs);qe=int(qe);nm=int(nm);al=int(al);ts=int(ts);te=int(te)
   ident=nm/al if al else 0
   if t not in truth or ident<a.min_identity or al<a.min_alignment:continue
   gid=truth[t]['genome_id'];iv[t].append((ts,te));best_identity[t]=max(best_identity[t],ident);queries[t].add(q);query_lengths[q]=ql;query_hits[q].append((qs,qe,gid,ident,t));total_accepted_alignment_bp+=al
 cont=[];gcov=defaultdict(int);gqueries=defaultdict(set)
 for sid,r in truth.items():
  cov=union_len(iv[sid]);gcov[r['genome_id']]+=cov;gqueries[r['genome_id']]|=queries[sid];cont.append({'sequence_id':sid,'genome_id':r['genome_id'],'length':r['length'],'covered_bp':cov,'covered_fraction':cov/int(r['length']),'best_identity':best_identity[sid],'assembly_contigs':len(queries[sid])})
 genomes=[]
 for g in sorted(glen):genomes.append({'genome_id':g,'length':glen[g],'covered_bp':gcov[g],'recovery_fraction':gcov[g]/glen[g] if glen[g] else 0,'assembly_contigs':len(gqueries[g]),'recovered_50':str(gcov[g]>=.5*glen[g]).lower(),'recovered_90':str(gcov[g]>=.9*glen[g]).lower()})
 chimera=[]
 for q,hits in sorted(query_hits.items()):
  by=defaultdict(list)
  for qs,qe,gid,ident,t in hits:by[gid].append((qs,qe))
  support=sorted(((union_len(v),g) for g,v in by.items()),reverse=True)
  second=support[1][0] if len(support)>1 else 0;ql=query_lengths[q];threshold=max(a.chimera_min_bp,int(a.chimera_min_query_fraction*ql));flag=len(support)>1 and second>=threshold
  chimera.append({'assembly_contig':q,'length':ql,'n_genome_ids':len(support),'primary_genome':support[0][1] if support else '','primary_aligned_bp':support[0][0] if support else 0,'secondary_genome':support[1][1] if len(support)>1 else '','secondary_aligned_bp':second,'chimera_threshold_bp':threshold,'cross_binid_chimera':str(flag).lower()})
 for name,rows,fields in [
  ('per_truth_contig.tsv',cont,list(cont[0]) if cont else ['sequence_id','genome_id','length','covered_bp','covered_fraction','best_identity','assembly_contigs']),
  ('per_genome_recovery.tsv',genomes,list(genomes[0]) if genomes else ['genome_id','length','covered_bp','recovery_fraction','assembly_contigs','recovered_50','recovered_90']),
  ('cross_binid_chimera_audit.tsv',chimera,list(chimera[0]) if chimera else ['assembly_contig','length','n_genome_ids','primary_genome','primary_aligned_bp','secondary_genome','secondary_aligned_bp','chimera_threshold_bp','cross_binid_chimera'])]:
  with open(out/name,'w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 vals=[r['recovery_fraction'] for r in genomes];unique_covered=sum(gcov.values());truth_bp=sum(glen.values());chimeric=[r for r in chimera if r['cross_binid_chimera']=='true'];summary={'n_genomes':len(genomes),'truth_bp_total':truth_bp,'unique_truth_bp_covered':unique_covered,'genome_fraction_percent':100*unique_covered/truth_bp if truth_bp else 0,'mean_genome_recovery':sum(vals)/len(vals) if vals else 0,'median_genome_recovery':sorted(vals)[len(vals)//2] if vals else 0,'genomes_recovered_50':sum(r['recovered_50']=='true' for r in genomes),'genomes_recovered_90':sum(r['recovered_90']=='true' for r in genomes),'accepted_alignment_bp':total_accepted_alignment_bp,'alignment_to_unique_truth_ratio':total_accepted_alignment_bp/unique_covered if unique_covered else None,'assembly_contigs_with_accepted_alignment':len(chimera),'cross_binid_chimeric_contigs':len(chimeric),'cross_binid_chimeric_bp':sum(r['length'] for r in chimeric),'min_identity':a.min_identity,'min_alignment':a.min_alignment,'truth_only_evaluation':True,'boundary':'Cross-BINID mappings are a conservative chimera proxy against the CAMI gold-standard assembly; they are not identical to the official MetaQUAST misassembly metric.'}
 (out/'GOLD_COVERAGE_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
