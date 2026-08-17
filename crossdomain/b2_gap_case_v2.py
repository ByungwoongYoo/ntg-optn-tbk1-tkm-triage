#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, os, time
from pathlib import Path
from ortools.sat.python import cp_model

N=100; K=14

def verify(A):
    if len(A)!=K or len(set(A))!=K:return False,{}
    c={t:0 for t in range(1,N)}
    for a in A:
        for b in A:
            if a==b: continue
            c[(a-b)%N]+=1
    return max(c.values(),default=0)<=2,c

def circular_dist(a,b):
    d=(a-b)%N
    return min(d,N-d)

def solve(g:int,limit:int):
    assert 1<=g<=7
    m=cp_model.CpModel(); x=[m.NewBoolVar(f'x_{i}') for i in range(N)]
    m.Add(sum(x)==K)
    # Exact symmetry-normalized case split: in every 14-subset of a 100-cycle,
    # the minimum cyclic gap is <= floor(100/14)=7. Rotate one minimum-gap pair to (0,g).
    m.Add(x[0]==1); m.Add(x[g]==1)
    # g is the minimum cyclic distance among selected points.
    if g>1:
        for a in range(N):
            for b in range(a+1,N):
                if circular_dist(a,b)<g:
                    m.Add(x[a]+x[b] <= 1)
    # Linearize unordered co-selection. Ordered difference multiplicities use the same product.
    y={}
    for a in range(N):
        for b in range(a+1,N):
            v=m.NewBoolVar(f'y_{a}_{b}'); y[a,b]=v
            m.Add(v<=x[a]); m.Add(v<=x[b]); m.Add(v>=x[a]+x[b]-1)
    for t in range(1,N):
        terms=[]
        for b in range(N):
            a=(b+t)%N
            u,v=sorted((a,b))
            terms.append(y[u,v])
        m.Add(sum(terms)<=2)
    # Safe reflection breaker because the normalized pair {0,g} is fixed as a set.
    # Compare presence at g+1 versus -1 only when this does not force either member of the normalized pair.
    if g<49:
        m.Add(x[(g+1)%N] >= x[99])
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=limit
    s.parameters.num_search_workers=min(16,os.cpu_count() or 4); s.parameters.random_seed=20260817+g
    s.parameters.log_search_progress=True
    st=s.Solve(m); name=s.StatusName(st)
    A=[i for i in range(N) if st in (cp_model.OPTIMAL,cp_model.FEASIBLE) and s.Value(x[i])]
    ok,c=verify(A) if A else (False,{})
    return name,A,ok,c,s.ResponseStats()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--gap',type=int,required=True);ap.add_argument('--limit',type=int,default=1500);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    t=time.time(); status,A,ok,c,stats=solve(a.gap,a.limit)
    res={'problem':'14-element B_2[2] subset of Z_100','case_definition':'minimum cyclic gap g after rotation of a minimum-gap pair to {0,g}',
         'gap':a.gap,'coverage_proof':'Every 14-subset of a 100-cycle has minimum cyclic gap <= floor(100/14)=7, hence g=1..7 exhaust all possibilities up to rotation.',
         'solver':'OR-Tools CP-SAT','status':status,'candidate':A,'verified_witness':ok,'difference_counts':c if ok else {},'wall_seconds':time.time()-t,'solver_stats':stats,
         'claim_boundary':'FEASIBLE + independent verifier is a decisive positive witness. INFEASIBLE closes only this normalized gap case computationally; publication-grade global nonexistence additionally requires independently checkable proof/certificate for all seven cases.'}
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
    if ok:(out/'WITNESS.txt').write_text(' '.join(map(str,A))+'\n',encoding='utf-8')
    print(json.dumps(res,indent=2))
if __name__=='__main__': main()
