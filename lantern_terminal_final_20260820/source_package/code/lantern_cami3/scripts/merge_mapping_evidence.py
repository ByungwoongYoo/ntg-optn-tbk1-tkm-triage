#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def args():
    p=argparse.ArgumentParser(); p.add_argument('--short-coverage',action='append',default=[],help='SAMPLE=FILE'); p.add_argument('--long-coverage',action='append',default=[]); p.add_argument('--long-paf',action='append',default=[]); p.add_argument('--end-margin',type=int,default=250); p.add_argument('--min-mapq',type=int,default=20); p.add_argument('--min-align',type=int,default=500); p.add_argument('--out',required=True); return p.parse_args()

def covfiles(specs,kind,rows):
    for spec in specs:
        sample,path=spec.split('=',1)
        with open(path) as f:
            rd=csv.DictReader(f,delimiter='\t')
            for r in rd:
                cid=r.get('#rname') or r.get('rname') or r.get('chrom')
                if not cid: continue
                key=(cid,sample); d=rows[key]; d['contig_id']=cid; d['sample_id']=sample
                d[f'{kind}_reads']=int(float(r.get('numreads',0) or 0)); d[f'{kind}_breadth']=float(r.get('coverage',0) or 0)/100.0; d[f'{kind}_depth']=float(r.get('meandepth',0) or 0)

def pafs(specs,rows,a):
    for spec in specs:
        sample,path=spec.split('=',1); seen=defaultdict(set)
        with open(path) as f:
            for line in f:
                x=line.rstrip().split('\t')
                if len(x)<12: continue
                q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12]; tl=int(tl); ts=int(ts); te=int(te); al=int(al); mq=int(mq)
                if mq<a.min_mapq or al<a.min_align: continue
                span=(ts<=a.end_margin and te>=tl-a.end_margin) or al>=.8*tl
                if span: seen[t].add(q)
        for cid,reads in seen.items():
            d=rows[(cid,sample)]; d['contig_id']=cid; d['sample_id']=sample; d['long_spanning_reads']=len(reads)

def main():
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); rows=defaultdict(dict)
    covfiles(a.short_coverage,'short',rows); covfiles(a.long_coverage,'long',rows); pafs(a.long_paf,rows,a)
    fields=['contig_id','sample_id','short_reads','short_breadth','short_depth','long_reads','long_breadth','long_depth','long_spanning_reads']
    vals=[{k:d.get(k,0) for k in fields} for _,d in sorted(rows.items())]
    with open(out/'mapping_evidence.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(vals)
    (out/'EVIDENCE_SUMMARY.json').write_text(json.dumps({'n_rows':len(vals),'n_contigs':len({r['contig_id'] for r in vals}),'n_samples':len({r['sample_id'] for r in vals})},indent=2)+'\n')
if __name__=='__main__':main()
