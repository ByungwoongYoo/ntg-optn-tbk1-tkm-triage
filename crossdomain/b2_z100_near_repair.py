#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,time
from collections import Counter
from pathlib import Path
from ortools.sat.python import cp_model

N=100; K=14
NEAR=[0,2,8,9,15,20,34,44,48,61,77,78,97,99]

def cd(i,j):
    d=abs(i-j)%N
    return min(d,N-d)

def verify(S):
    c=Counter((a-b)%N for a in S for b in S if a!=b)
    bad={d:n for d,n in c.items() if n>2}
    return len(S)==K and len(set(S))==K and not bad, bad, max(c.values(),default=0)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--max-replacements',type=int,required=True); ap.add_argument('--seconds',type=float,default=300); ap.add_argument('--out-dir',required=True); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    model=cp_model.CpModel(); x=[model.NewBoolVar(f'x{i}') for i in range(N)]
    model.Add(sum(x)==K)
    # Both sets have size 14. At most r replacements means overlap with the heuristic near-set >= 14-r.
    model.Add(sum(x[i] for i in NEAR)>=K-a.max_replacements)
    yd={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            y=model.NewBoolVar(f'y{i}_{j}')
            model.Add(y<=x[i]); model.Add(y<=x[j]); model.Add(y>=x[i]+x[j]-1)
            yd[cd(i,j)].append(y)
    for d in range(1,50): model.Add(sum(yd[d])<=2)
    model.Add(sum(yd[50])<=1)
    # Preserve the prior candidate as a search hint only; the exact constraints decide validity.
    ns=set(NEAR)
    for i in range(N): model.AddHint(x[i],1 if i in ns else 0)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=a.seconds; s.parameters.num_search_workers=8; s.parameters.random_seed=20260817+a.max_replacements
    t=time.time(); st=s.Solve(model); elapsed=time.time()-t; name=s.StatusName(st)
    S=[]; check=None
    if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        S=[i for i in range(N) if s.Value(x[i])]; ok,bad,mx=verify(S)
        check={'verified':ok,'violations':bad,'max_ordered_difference_multiplicity':mx,'overlap_with_near':len(set(S)&set(NEAR)),'replacements':K-len(set(S)&set(NEAR))}
    res={'problem':'B2[2] size-14 subset of Z100','near_set':NEAR,'max_replacements':a.max_replacements,'status':name,'elapsed_seconds':elapsed,'branches':s.NumBranches(),'conflicts':s.NumConflicts(),'solution':S,'verification':check,'claim_boundary':'A verified feasible solution is a global existence witness. INFEASIBLE only rules out the specified Hamming neighborhood of the heuristic near-set and says nothing about the remaining search space.'}
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8'); print(json.dumps(res,indent=2),flush=True)
if __name__=='__main__': main()
