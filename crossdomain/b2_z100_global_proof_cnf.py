#!/usr/bin/env python3
"""Independent DIMACS generator for proof-producing SAT runs.

The default instance is global after the elementary affine normalization
{0,1} subset A.  Optional --third/--fourth arguments create the same disjoint
reflection-reduced branches used by the exact search, but no production search
code or PySAT dependency is used.

Variables:
  x_i       selected residue i
  y_i_j     conjunction x_i and x_j for each unordered pair
  sequential-counter auxiliaries

Constraints:
  |A| = 14;
  x_0 = x_1 = 1;
  each circular distance 1..49 occurs among unordered pairs at most twice;
  distance 50 occurs at most once.
These are exactly equivalent to r_A(t) <= 2 for all nonzero ordered t.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

N=100; K=14

class CNF:
    def __init__(self): self.nvars=0; self.clauses=[]; self.names={}
    def var(self,name):
        self.nvars+=1; self.names[name]=self.nvars; return self.nvars
    def add(self,*lits): self.clauses.append(tuple(int(x) for x in lits))
    def atmost(self,lits,k,label):
        lits=list(lits); n=len(lits)
        if k<0: self.add(); return
        if k>=n: return
        if k==0:
            for x in lits:self.add(-x)
            return
        # Sinz sequential counter. s[i,j] means at least j+1 true among lits[:i+1].
        s=[[self.var(f"{label}_s_{i}_{j}") for j in range(k)] for i in range(n-1)]
        self.add(-lits[0],s[0][0])
        for j in range(1,k): self.add(-s[0][j])
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
            for c in self.clauses: f.write(" ".join(map(str,c))+" 0\n")

def sha(path):
    h=hashlib.sha256();
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def build(third=None,fourth=None):
    if fourth is not None and third is None: raise ValueError('--fourth requires --third')
    if third is not None and not (2<=third<=45): raise ValueError('third must be 2..45')
    if fourth is not None and not (third+1<=fourth<=91-third): raise ValueError('fourth outside exact branch range')
    cnf=CNF(); x=[cnf.var(f"x_{i}") for i in range(N)]
    cnf.equals(x,K,'card14')
    cnf.add(x[0]);cnf.add(x[1])
    if third is not None:
        cnf.add(x[third])
        for i in range(2,third):cnf.add(-x[i])
        # Chosen reflection orientation: wrap >= third-1, hence largest <= 101-third.
        for i in range(102-third,N):cnf.add(-x[i])
    if fourth is not None:
        cnf.add(x[fourth])
        for i in range(third+1,fourth):cnf.add(-x[i])
    buckets={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            y=cnf.var(f"y_{i}_{j}")
            cnf.add(-y,x[i]);cnf.add(-y,x[j]);cnf.add(y,-x[i],-x[j])
            d=min(j-i,N-(j-i));buckets[d].append(y)
    for d in range(1,50):cnf.atmost(buckets[d],2,f"dist_{d}")
    cnf.atmost(buckets[50],1,"dist_50")
    return cnf,x,buckets

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--metadata',type=Path);ap.add_argument('--third',type=int);ap.add_argument('--fourth',type=int)
    a=ap.parse_args();cnf,x,b=build(a.third,a.fourth);a.output.parent.mkdir(parents=True,exist_ok=True);cnf.write(a.output)
    meta={"schema":"b2-z100-independent-sat-v1","modulus":N,"target_size":K,"normalization":[0,1],"third":a.third,"fourth":a.fourth,"variables":cnf.nvars,"clauses":len(cnf.clauses),"cnf_bytes":a.output.stat().st_size,"cnf_sha256":sha(a.output),"distance_bucket_sizes":{str(d):len(v) for d,v in b.items()},"semantics":"UNSAT proves no normalized 14-set in this scope; a verified SAT assignment gives a counterexample.","generator":"standalone Python; no project exact-search code and no PySAT"}
    m=a.metadata or a.output.with_suffix(a.output.suffix+'.json');m.write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
    print(json.dumps(meta,indent=2,sort_keys=True))
if __name__=='__main__':main()
