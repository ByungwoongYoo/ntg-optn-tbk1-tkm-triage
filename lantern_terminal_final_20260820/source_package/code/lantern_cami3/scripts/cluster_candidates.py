#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
from fasta_utils import read_fasta,write_record,seq_sha

class UF:
    def __init__(self,xs): self.p={x:x for x in xs}; self.r={x:0 for x in xs}
    def find(self,x):
        while self.p[x]!=x: self.p[x]=self.p[self.p[x]]; x=self.p[x]
        return x
    def union(self,a,b):
        a=self.find(a); b=self.find(b)
        if a==b:return
        if self.r[a]<self.r[b]:a,b=b,a
        self.p[b]=a
        if self.r[a]==self.r[b]:self.r[a]+=1

def args():
    p=argparse.ArgumentParser(); p.add_argument('--fasta',required=True); p.add_argument('--metadata',required=True); p.add_argument('--paf',required=True)
    p.add_argument('--min-identity',type=float,default=.97); p.add_argument('--min-shorter-coverage',type=float,default=.85); p.add_argument('--out',required=True); return p.parse_args()

def main():
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    seqs=dict(read_fasta(a.fasta)); meta={}
    with open(a.metadata,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'): meta[r['candidate_id']]=r
    if set(seqs)!=set(meta): raise ValueError(f'FASTA/metadata mismatch {len(seqs)} {len(meta)}')
    uf=UF(seqs); accepted=0
    with open(a.paf) as f:
        for line in f:
            x=line.rstrip().split('\t')
            if len(x)<12: continue
            q,ql,qs,qe,st,t,tl,ts,te,nm,al,mq=x[:12]
            if q==t or q not in seqs or t not in seqs: continue
            ql=int(ql); tl=int(tl); nm=int(nm); al=int(al)
            ident=nm/al if al else 0; cov=al/min(ql,tl) if min(ql,tl) else 0
            if ident>=a.min_identity and cov>=a.min_shorter_coverage:
                uf.union(q,t); accepted+=1
    groups=defaultdict(list)
    for x in seqs: groups[uf.find(x)].append(x)
    ordered=sorted(groups.values(),key=lambda g:(-max(len(seqs[x]) for x in g),sorted(g)[0]))
    member_rows=[]; rep_rows=[]
    with open(out/'cluster_representatives.fasta','w') as fo:
        for g in ordered:
            best=max(g,key=lambda x:(len(seqs[x]),-seqs[x].count('N'),x))
            stable='LANTERN_C'+seq_sha(seqs[best])[:16].upper()
            write_record(fo,stable,seqs[best])
            source_ids=sorted({meta[x]['source_id'] for x in g}); assemblers=sorted({meta[x].get('assembler','') for x in g if meta[x].get('assembler','')}); scopes=sorted({meta[x].get('scope','') for x in g if meta[x].get('scope','')}); modes=sorted({meta[x].get('mode','') for x in g if meta[x].get('mode','')})
            rep_rows.append({'representative_id':stable,'selected_candidate_id':best,'length':len(seqs[best]),'n_fraction':seqs[best].count('N')/len(seqs[best]),'cluster_size':len(g),'source_count':len(source_ids),'assembler_count':len(assemblers),'source_ids':','.join(source_ids),'assemblers':','.join(assemblers),'scopes':','.join(scopes),'modes':','.join(modes),'sequence_sha256':seq_sha(seqs[best])})
            for x in sorted(g): member_rows.append({'representative_id':stable,'candidate_id':x,'source_id':meta[x]['source_id'],'original_id':meta[x]['original_id'],'length':len(seqs[x]),'is_selected':str(x==best).lower()})
    for name,rows in [('cluster_members.tsv',member_rows),('representative_metadata.tsv',rep_rows)]:
        with open(out/name,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    (out/'CLUSTER_SUMMARY.json').write_text(json.dumps({'n_input':len(seqs),'n_clusters':len(ordered),'n_edges_accepted':accepted,'min_identity':a.min_identity,'min_shorter_coverage':a.min_shorter_coverage},indent=2)+'\n')
if __name__=='__main__':main()
