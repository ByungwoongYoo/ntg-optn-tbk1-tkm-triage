#!/usr/bin/env python3
"""Conservative two-contig scaffolding from long-read end bridges."""
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
from fasta_utils import read_fasta,write_record,revcomp

def args():
 p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--paf',action='append',required=True);p.add_argument('--config',required=True);p.add_argument('--out',required=True);return p.parse_args()

def canon_edge(a,ao,b,bo):
 x=(a,ao,b,bo); y=(b,'-' if bo=='+' else '+',a,'-' if ao=='+' else '+')
 return min(x,y)

def main():
 a=args();cfg=json.load(open(a.config));out=Path(a.out);out.mkdir(parents=True,exist_ok=True);seqs=dict(read_fasta(a.fasta));support=defaultdict(set);gaps=defaultdict(list)
 for path in a.paf:
  byread=defaultdict(list)
  with open(path) as f:
   for line in f:
    x=line.rstrip().split('\t')
    if len(x)<12:continue
    q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12];qs=int(qs);qe=int(qe);tl=int(tl);ts=int(ts);te=int(te);al=int(al);mq=int(mq)
    if t not in seqs or mq<cfg['bridge_min_mapq'] or al<500:continue
    left=ts<=cfg['bridge_end_margin'];right=te>=tl-cfg['bridge_end_margin']
    if not (left or right):continue
    byread[q].append({'qstart':qs,'qend':qe,'strand':st,'target':t,'left':left,'right':right})
  for read,hits in byread.items():
   best={}
   for h in hits:
    span=h['qend']-h['qstart']
    if h['target'] not in best or span>best[h['target']]['qend']-best[h['target']]['qstart']:best[h['target']]=h
   hs=sorted(best.values(),key=lambda h:h['qstart'])
   if len(hs)!=2:continue
   u,v=hs
   if u['qend']>v['qstart']+200:continue
   u_exit=u['right'] if u['strand']=='+' else u['left'];v_entry=v['left'] if v['strand']=='+' else v['right']
   if not (u_exit and v_entry):continue
   e=canon_edge(u['target'],u['strand'],v['target'],v['strand']);support[e].add(read);gaps[e].append(max(10,v['qstart']-u['qend']))
 candidates=[]
 for e,reads in support.items():
  if len(reads)>=cfg['bridge_min_reads']:candidates.append((len(reads),e,sorted(gaps[e])))
 end_edges=defaultdict(list)
 for n,e,gs in candidates:
  a1,o1,b1,o2=e;end_edges[(a1,'R' if o1=='+' else 'L')].append((n,e));end_edges[(b1,'L' if o2=='+' else 'R')].append((n,e))
 accepted=[]
 for n,e,gs in sorted(candidates,reverse=True):
  a1,o1,b1,o2=e;ends=[(a1,'R' if o1=='+' else 'L'),(b1,'L' if o2=='+' else 'R')];conflict=False
  for end in ends:
   rivals=sorted(end_edges[end],reverse=True)
   if rivals and rivals[0][1]!=e and rivals[0][0]>cfg['bridge_max_conflict_reads']:conflict=True
   if len(rivals)>1 and rivals[0][1]==e and rivals[1][0]>cfg['bridge_max_conflict_reads']:conflict=True
  if not conflict:accepted.append((n,e,gs))
 used=set();chosen=[]
 for n,e,gs in sorted(accepted,reverse=True):
  if e[0] in used or e[2] in used:continue
  used|={e[0],e[2]};chosen.append((n,e,gs))
 with open(out/'LANTERN_SCAFFOLDED.fasta','w') as f:
  for i,(n,(x,xo,y,yo),gs) in enumerate(chosen,1):
   sx=seqs[x] if xo=='+' else revcomp(seqs[x]);sy=seqs[y] if yo=='+' else revcomp(seqs[y]);gap=min(10000,max(10,round(sum(gs)/len(gs)))) if gs else cfg['bridge_gap_n_default']
   write_record(f,f'LANTERN_SCAFFOLD_{i:06d}_support_{n}',sx+'N'*gap+sy)
  for x,s in seqs.items():
   if x not in used:write_record(f,x,s)
 rows=[]
 for n,e,gs in candidates:rows.append({'contig_a':e[0],'orientation_a':e[1],'contig_b':e[2],'orientation_b':e[3],'supporting_reads':n,'median_gap':gs[len(gs)//2] if gs else '', 'accepted':str(any(e==z[1] for z in chosen)).lower()})
 with open(out/'BRIDGE_AUDIT.tsv','w',newline='') as f:
  fields=['contig_a','orientation_a','contig_b','orientation_b','supporting_reads','median_gap','accepted'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 (out/'SCAFFOLD_SUMMARY.json').write_text(json.dumps({'n_input':len(seqs),'n_candidate_edges':len(candidates),'n_chosen_edges':len(chosen),'n_output':len(seqs)-len(chosen),'boundary':'Conservative long-read scaffolds remain hypotheses until gold-standard evaluation; the unscaffolded assembly is retained as an ablation.'},indent=2)+'\n')
if __name__=='__main__':main()
