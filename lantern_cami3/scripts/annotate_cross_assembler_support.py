#!/usr/bin/env python3
"""Annotate strict cluster representatives with looser independent-assembler support."""
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();p.add_argument('--representatives',required=True);p.add_argument('--members',required=True);p.add_argument('--candidate-metadata',required=True);p.add_argument('--paf',required=True);p.add_argument('--min-identity',type=float,default=.95);p.add_argument('--min-shorter-coverage',type=float,default=.50);p.add_argument('--out',required=True);return p.parse_args()

def main():
    a=parse_args();out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);reps=list(csv.DictReader(open(a.representatives),delimiter='\t'));members=list(csv.DictReader(open(a.members),delimiter='\t'));meta={r['candidate_id']:r for r in csv.DictReader(open(a.candidate_metadata),delimiter='\t')};c2r={r['candidate_id']:r['representative_id'] for r in members};support_asm=defaultdict(set);support_src=defaultdict(set);support_reps=defaultdict(set)
    for r in members:
        m=meta[r['candidate_id']];support_asm[r['representative_id']].add(m.get('assembler',''));support_src[r['representative_id']].add(m.get('source_id',''))
    accepted=0
    with open(a.paf) as f:
        for line in f:
            x=line.rstrip().split('\t')
            if len(x)<12:continue
            q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12]
            if q==t or q not in c2r or t not in c2r:continue
            ql=int(ql);tl=int(tl);nm=int(nm);al=int(al);ident=nm/al if al else 0;cov=al/min(ql,tl) if min(ql,tl) else 0
            if ident<a.min_identity or cov<a.min_shorter_coverage:continue
            qr,tr=c2r[q],c2r[t]
            if qr==tr:continue
            mt=meta[t];mqm=meta[q];support_asm[qr].add(mt.get('assembler',''));support_src[qr].add(mt.get('source_id',''));support_reps[qr].add(tr);support_asm[tr].add(mqm.get('assembler',''));support_src[tr].add(mqm.get('source_id',''));support_reps[tr].add(qr);accepted+=1
    rows=[]
    for r in reps:
        rid=r['representative_id'];asm=sorted(x for x in support_asm[rid] if x);src=sorted(x for x in support_src[rid] if x);x=dict(r);x['cluster_assembler_count']=r.get('assembler_count','');x['support_assembler_count']=len(asm);x['support_source_count']=len(src);x['supporting_representative_count']=len(support_reps[rid]);x['support_assemblers']=','.join(asm);x['support_sources']=','.join(src);rows.append(x)
    with open(out,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    report={'n_representatives':len(rows),'accepted_cross_cluster_alignment_edges':accepted,'min_identity':a.min_identity,'min_shorter_coverage':a.min_shorter_coverage,'truth_blind':True,'boundary':'Support is derived only from independent assembler sequences and never from gold truth or taxonomic references.'};out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
