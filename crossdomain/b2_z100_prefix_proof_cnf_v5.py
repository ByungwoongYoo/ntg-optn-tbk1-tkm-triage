#!/usr/bin/env python3
"""Independent proof-CNF generator with mathematically proven propagation cuts.

The canonical prefix semantics and original distance constraints are identical
to `b2_z100_prefix_proof_cnf_v3.py`. Additional parity and quotient-capacity
constraints are consequences of the original B2[2] condition; see
`B2_Z100_REDUNDANT_CUTS_PROOF.md`.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

N=100
K=14
MODS=(2,4,5,10,20,25,50)

class CNF:
    def __init__(self):
        self.nvars=0;self.clauses=[];self.names={}
    def var(self,name):
        self.nvars+=1;self.names[name]=self.nvars;return self.nvars
    def add(self,*lits):self.clauses.append(tuple(int(x) for x in lits))
    def atmost(self,lits,k,label):
        lits=list(lits);n=len(lits)
        if k<0:self.add();return
        if k>=n:return
        if k==0:
            for x in lits:self.add(-x)
            return
        s=[[self.var(f'{label}_s_{i}_{j}') for j in range(k)] for i in range(n-1)]
        self.add(-lits[0],s[0][0])
        for j in range(1,k):self.add(-s[0][j])
        for i in range(1,n-1):
            self.add(-lits[i],s[i][0]);self.add(-s[i-1][0],s[i][0])
            for j in range(1,k):
                self.add(-lits[i],-s[i-1][j-1],s[i][j]);self.add(-s[i-1][j],s[i][j])
            self.add(-lits[i],-s[i-1][k-1])
        self.add(-lits[-1],-s[n-2][k-1])
    def equals(self,lits,k,label):
        self.atmost(lits,k,label+'_le');self.atmost([-x for x in lits],len(lits)-k,label+'_ge')
    def write(self,path):
        with path.open('w',encoding='ascii',newline='\n') as f:
            f.write(f'p cnf {self.nvars} {len(self.clauses)}\n')
            for c in self.clauses:f.write(' '.join(map(str,c))+' 0\n')

def fsha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def validate_prefix(prefix):
    if len(prefix)<3 or len(prefix)>K:raise ValueError('prefix length must be 3..14')
    if prefix[:2]!=[0,1]:raise ValueError('prefix must begin [0,1]')
    if any(type(x) is not int for x in prefix):raise ValueError('prefix entries must be integers')
    if any(x<0 or x>=N for x in prefix):raise ValueError('prefix outside representatives 0..99')
    if any(prefix[i]>=prefix[i+1] for i in range(len(prefix)-1)):raise ValueError('prefix must increase')
    t=prefix[2]
    if not 2<=t<=45:raise ValueError('third must be 2..45')
    max_selected=101-t
    if prefix[-1]>max_selected-(K-len(prefix)):raise ValueError('prefix cannot be completed under reflection bound')
    return t,max_selected

def build(prefix):
    t,max_selected=validate_prefix(prefix)
    c=CNF();x=[c.var(f'x_{i}') for i in range(N)];c.equals(x,K,'card14')
    c.add(x[0]);c.add(x[1]);prev=1
    for p in prefix[2:]:
        c.add(x[p])
        for i in range(prev+1,p):c.add(-x[i])
        prev=p
    for i in range(max_selected+1,N):c.add(-x[i])

    # Every valid 14-set has 5..9 even selected residues.
    evens=[x[i] for i in range(0,N,2)]
    c.atmost(evens,9,'even_atmost9')
    c.atmost([-z for z in evens],45,'even_atleast5')

    dist={d:[] for d in range(1,51)};pairs={}
    for i in range(N):
        for j in range(i+1,N):
            y=c.var(f'y_{i}_{j}');pairs[(i,j)]=y
            c.add(-y,x[i]);c.add(-y,x[j]);c.add(y,-x[i],-x[j])
            d=min(j-i,N-(j-i));dist[d].append(y)
    for d in range(1,50):c.atmost(dist[d],2,f'dist_{d}')
    c.atmost(dist[50],1,'dist_50')

    quotient_meta={}
    for m in MODS:
        qb={q:[] for q in range(0,m//2+1)}
        for (i,j),y in pairs.items():
            r=(j-i)%m;q=min(r,(m-r)%m);qb[q].append(y)
        bounds={0:100//m-1}
        c.atmost(qb[0],bounds[0],f'q{m}_0')
        for q in range(1,m//2+1):
            bound=100//m if (m%2==0 and q==m//2) else 2*(100//m)
            bounds[q]=bound;c.atmost(qb[q],bound,f'q{m}_{q}')
        quotient_meta[str(m)]={'bucket_sizes':{str(q):len(v) for q,v in qb.items()},'bounds':{str(q):b for q,b in bounds.items()}}
    return c,dist,max_selected,quotient_meta

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prefix-json',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--metadata',type=Path)
    a=ap.parse_args();prefix=json.loads(a.prefix_json);c,dist,max_selected,qm=build(prefix)
    a.output.parent.mkdir(parents=True,exist_ok=True);c.write(a.output)
    remaining=K-len(prefix);nr=None
    if remaining:nr=[prefix[-1]+1,max_selected-(remaining-1)]
    meta={
      'schema':'b2-z100-independent-prefix-sat-v5','modulus':N,'target_size':K,
      'normalization':[0,1],'prefix':prefix,'third':prefix[2],
      'reflection_max_selected':max_selected,'remaining_points':remaining,
      'next_selected_range':nr,'variables':c.nvars,'clauses':len(c.clauses),
      'cnf_bytes':a.output.stat().st_size,'cnf_sha256':fsha(a.output),
      'distance_bucket_sizes':{str(d):len(v) for d,v in dist.items()},
      'redundant_cuts':{'even_selected_range':[5,9],'quotient_moduli':list(MODS),'proof_document':'crossdomain/B2_Z100_REDUNDANT_CUTS_PROOF.md','quotient_data':qm},
      'semantics':'UNSAT excludes exactly this canonical prefix branch because all added cuts are proven consequences; verified SAT yields a 14-set counterexample.',
      'generator':'standalone Python; no optimized exact-search code and no PySAT'}
    m=a.metadata or a.output.with_suffix(a.output.suffix+'.json');m.write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
    print(json.dumps(meta,indent=2,sort_keys=True))
if __name__=='__main__':main()
