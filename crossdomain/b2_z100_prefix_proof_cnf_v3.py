#!/usr/bin/env python3
"""Standalone proof-CNF generator for a canonical prefix branch of B2[2] in Z/100Z.

The branch always contains the affine normalization 0,1. A prefix
[0,1,t,u,v,...] means these are exactly the first selected residues in
increasing order. The reflection orientation is fixed by
max(A) <= 101-t. The generated CNF is independent of the optimized DFS.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

N=100
K=14

class CNF:
    def __init__(self):
        self.nvars=0
        self.clauses=[]
        self.names={}
    def var(self,name):
        self.nvars+=1
        self.names[name]=self.nvars
        return self.nvars
    def add(self,*lits):
        self.clauses.append(tuple(int(x) for x in lits))
    def atmost(self,lits,k,label):
        lits=list(lits)
        n=len(lits)
        if k<0:
            self.add()
            return
        if k>=n:
            return
        if k==0:
            for x in lits:
                self.add(-x)
            return
        s=[[self.var(f"{label}_s_{i}_{j}") for j in range(k)] for i in range(n-1)]
        self.add(-lits[0],s[0][0])
        for j in range(1,k):
            self.add(-s[0][j])
        for i in range(1,n-1):
            self.add(-lits[i],s[i][0])
            self.add(-s[i-1][0],s[i][0])
            for j in range(1,k):
                self.add(-lits[i],-s[i-1][j-1],s[i][j])
                self.add(-s[i-1][j],s[i][j])
            self.add(-lits[i],-s[i-1][k-1])
        self.add(-lits[-1],-s[n-2][k-1])
    def equals(self,lits,k,label):
        self.atmost(lits,k,label+"_le")
        self.atmost([-x for x in lits],len(lits)-k,label+"_ge")
    def write(self,path):
        with path.open('w',encoding='ascii',newline='\n') as f:
            f.write(f"p cnf {self.nvars} {len(self.clauses)}\n")
            for c in self.clauses:
                f.write(" ".join(map(str,c))+" 0\n")

def file_sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

def validate_prefix(prefix):
    if len(prefix)<3 or len(prefix)>K:
        raise ValueError('prefix length must be 3..14')
    if prefix[:2] != [0,1]:
        raise ValueError('prefix must begin [0,1]')
    if any(not isinstance(x,int) for x in prefix):
        raise ValueError('prefix entries must be integers')
    if any(x<0 or x>=N for x in prefix):
        raise ValueError('prefix entries outside Z_100 representatives')
    if any(prefix[i]>=prefix[i+1] for i in range(len(prefix)-1)):
        raise ValueError('prefix must be strictly increasing')
    t=prefix[2]
    if not 2<=t<=45:
        raise ValueError('third element must be 2..45')
    max_allowed=101-t
    if prefix[-1]>max_allowed-(K-len(prefix)):
        raise ValueError('prefix leaves too few positions under reflection bound')
    return t,max_allowed

def build(prefix):
    t,max_allowed=validate_prefix(prefix)
    cnf=CNF()
    x=[cnf.var(f"x_{i}") for i in range(N)]
    cnf.equals(x,K,'card14')
    cnf.add(x[0])
    cnf.add(x[1])
    prev=1
    for p in prefix[2:]:
        cnf.add(x[p])
        for i in range(prev+1,p):
            cnf.add(-x[i])
        prev=p
    for i in range(max_allowed+1,N):
        cnf.add(-x[i])
    buckets={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            y=cnf.var(f"y_{i}_{j}")
            cnf.add(-y,x[i])
            cnf.add(-y,x[j])
            cnf.add(y,-x[i],-x[j])
            d=min(j-i,N-(j-i))
            buckets[d].append(y)
    for d in range(1,50):
        cnf.atmost(buckets[d],2,f"dist_{d}")
    cnf.atmost(buckets[50],1,'dist_50')
    return cnf,buckets,max_allowed

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--prefix-json',required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--metadata',type=Path)
    a=ap.parse_args()
    prefix=json.loads(a.prefix_json)
    cnf,buckets,max_allowed=build(prefix)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    cnf.write(a.output)
    remaining=K-len(prefix)
    next_range=None
    if remaining:
        lo=prefix[-1]+1
        hi=max_allowed-(remaining-1)
        next_range=[lo,hi]
    meta={
      'schema':'b2-z100-independent-prefix-sat-v3',
      'modulus':N,
      'target_size':K,
      'normalization':[0,1],
      'prefix':prefix,
      'third':prefix[2],
      'reflection_max_selected':max_allowed,
      'remaining_points':remaining,
      'next_selected_range':next_range,
      'variables':cnf.nvars,
      'clauses':len(cnf.clauses),
      'cnf_bytes':a.output.stat().st_size,
      'cnf_sha256':file_sha(a.output),
      'distance_bucket_sizes':{str(d):len(v) for d,v in buckets.items()},
      'semantics':'UNSAT excludes exactly this canonical prefix branch; verified SAT yields a 14-set counterexample.',
      'generator':'standalone Python; no optimized exact-search code and no PySAT'
    }
    m=a.metadata or a.output.with_suffix(a.output.suffix+'.json')
    m.write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps(meta,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
