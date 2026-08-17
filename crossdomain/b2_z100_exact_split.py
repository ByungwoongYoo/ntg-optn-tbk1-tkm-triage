#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, time
from collections import Counter
from pathlib import Path
from ortools.sat.python import cp_model

N=100
K=14
NEAR=[0,2,8,9,15,20,34,44,48,61,77,78,97,99]


def circ_dist(i,j):
    d=abs(j-i)%N
    return min(d,N-d)


def ordered_difference_counts(S):
    c=Counter()
    for a in S:
        for b in S:
            if a!=b:
                c[(a-b)%N]+=1
    return c


def verify(S):
    if len(S)!=K or len(set(S))!=K:
        return False, {'reason':'wrong_cardinality'}
    c=ordered_difference_counts(S)
    bad={d:n for d,n in c.items() if d!=0 and n>2}
    return not bad, {'max_nonzero_ordered_difference_multiplicity':max(c.values(),default=0),'violations':bad}


def solve(m:int, seconds:float, workers:int, out:Path):
    model=cp_model.CpModel()
    x=[model.NewBoolVar(f'x_{i}') for i in range(N)]
    model.Add(sum(x)==K)

    # Every 14-point subset of a 100-cycle has an adjacent cyclic gap <= floor(100/14)=7.
    # Translate a minimum-gap pair so its endpoints are 0 and m. Cases m=1..7 cover every possible set.
    model.Add(x[0]==1)
    model.Add(x[m]==1)
    for i in range(1,m):
        model.Add(x[i]==0)

    # Exact minimum circular distance m: no selected pair may be closer than m.
    if m>1:
        for i in range(N):
            for j in range(i+1,N):
                if circ_dist(i,j)<m:
                    model.Add(x[i]+x[j] <= 1)

    bydist={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            d=circ_dist(i,j)
            y=model.NewBoolVar(f'y_{i}_{j}')
            model.Add(y<=x[i]); model.Add(y<=x[j]); model.Add(y>=x[i]+x[j]-1)
            bydist[d].append(y)

    # For d=1..49, each selected unordered pair at circular distance d contributes once
    # to ordered difference d and once to 100-d. Hence <=2 pairs per distance class.
    for d in range(1,50):
        model.Add(sum(bydist[d])<=2)
    # At distance 50, a pair contributes twice to the single ordered difference 50.
    model.Add(sum(bydist[50])<=1)

    # A strong near-witness from the preceding 6.34e9-move heuristic is a safe hint only in m=1.
    if m==1:
        ns=set(NEAR)
        for i in range(N):
            model.AddHint(x[i],1 if i in ns else 0)

    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=workers
    solver.parameters.random_seed=20260817+m
    solver.parameters.cp_model_presolve=True
    solver.parameters.log_search_progress=False

    t=time.time(); status=solver.Solve(model); elapsed=time.time()-t
    name=solver.StatusName(status)
    S=[]
    check=None
    if status in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        S=[i for i in range(N) if solver.Value(x[i])]
        ok,detail=verify(S)
        check={'verified':ok,**detail}
    result={
        'problem':'Find a size-14 B2[2] subset of Z_100: every nonzero ordered difference has multiplicity <=2.',
        'symmetry_case_min_gap':m,
        'status':name,
        'elapsed_seconds':elapsed,
        'wall_time_seconds_reported':solver.WallTime(),
        'branches':solver.NumBranches(),
        'conflicts':solver.NumConflicts(),
        'solution':S,
        'verification':check,
        'case_interpretation':(
            'A verified FEASIBLE/OPTIMAL solution is a valid witness for the open instance.' if S else
            'INFEASIBLE eliminates this exact min-gap case under the encoding. UNKNOWN proves nothing.'
        ),
        'global_claim_boundary':'The seven cases m=1..7 cover all 14-point subsets by translation of a minimum cyclic gap. A positive witness must pass independent difference verification. A global impossibility claim would require all seven cases INFEASIBLE and an independent SAT/checker reproduction before being treated as a resolution.'
    }
    out.mkdir(parents=True,exist_ok=True)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
    return 0


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--min-gap',type=int,required=True,choices=range(1,8))
    ap.add_argument('--seconds',type=float,default=900)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--out-dir',required=True)
    a=ap.parse_args()
    raise SystemExit(solve(a.min_gap,a.seconds,a.workers,Path(a.out_dir)))

if __name__=='__main__':
    main()
