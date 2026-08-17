#!/usr/bin/env python3
from __future__ import annotations
import json, os, random, time
from pathlib import Path
from ortools.sat.python import cp_model

N=100; K=14

def verify(A):
    if len(A)!=K or len(set(A))!=K:return False,{}
    counts={t:0 for t in range(1,N)}
    for a in A:
        for b in A:
            if a==b:continue
            t=(a-b)%N
            if t: counts[t]+=1
    return max(counts.values(),default=0)<=2, counts

def cp_search(limit=900):
    model=cp_model.CpModel()
    x=[model.NewBoolVar(f'x_{i}') for i in range(N)]
    model.Add(sum(x)==K); model.Add(x[0]==1)
    # y[a,b] for unordered? r(t) uses ordered pairs. Product x_a*x_b shared across both orientations.
    y={}
    for a in range(N):
        for b in range(a+1,N):
            v=model.NewBoolVar(f'y_{a}_{b}'); y[a,b]=v
            model.Add(v<=x[a]); model.Add(v<=x[b]); model.Add(v>=x[a]+x[b]-1)
    for t in range(1,N):
        terms=[]
        for b in range(N):
            a=(b+t)%N
            if a==b:continue
            u,v=sorted((a,b))
            terms.append(y[u,v])
        # for t=50 each unordered pair appears twice in ordered enumeration; this is correct for r_A(50).
        model.Add(sum(terms)<=2)
    # Symmetry breaker under reflection A -> -A: require weighted bit representation <= reflected one.
    # Use a simple safe condition: x[1] >= x[99]; every reflection orbit has a representative satisfying it.
    model.Add(x[1]>=x[99])
    solver=cp_model.CpSolver(); solver.parameters.max_time_in_seconds=limit
    solver.parameters.num_search_workers=min(16,os.cpu_count() or 4); solver.parameters.random_seed=20260817
    solver.parameters.log_search_progress=True
    status=solver.Solve(model)
    A=[i for i in range(N) if status in (cp_model.OPTIMAL,cp_model.FEASIBLE) and solver.Value(x[i])]
    return solver.StatusName(status),A,solver.ResponseStats()

def local_search(seconds=180):
    rng=random.Random(20260817); best=None; bestcost=10**9; start=time.time()
    def cost(A):
        ok,c=verify(A); return sum(max(0,v-2)**2 for v in c.values()), max(c.values()) if c else 0
    while time.time()-start<seconds:
        A=set([0]+rng.sample(range(1,N),K-1)); c,m=cost(A)
        T=5.0
        for step in range(20000):
            if c==0:return sorted(A),0,m
            old=rng.choice(list(A-{0})); new=rng.choice([i for i in range(1,N) if i not in A])
            B=set(A);B.remove(old);B.add(new);c2,m2=cost(B)
            if c2<c or rng.random()<pow(2.718281828,-max(0,c2-c)/max(T,0.05)):
                A,c,m=B,c2,m2
            T*=0.9997
            if c<bestcost:bestcost=c;best=(sorted(A),c,m)
        if bestcost<=1:break
    return best if best else ([],None,None)

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/b2_z100'));out.mkdir(parents=True,exist_ok=True)
    lsA,lsc,lsm=local_search(120)
    cpstatus,cpA,stats=cp_search(1200)
    candidate=cpA if cpA else lsA
    valid,counts=verify(candidate) if candidate else (False,{})
    result={'problem':'14-element B_2[2] subset of Z_100','local_search':{'candidate':lsA,'cost':lsc,'max_r':lsm},
            'cp_sat_status':cpstatus,'cp_candidate':cpA,'cp_stats':stats,'final_candidate':candidate,'verified':valid,
            'difference_counts':counts if valid else {},
            'claim_boundary':'A verified 14-set is a decisive positive witness. CP-SAT UNKNOWN is not evidence of nonexistence; INFEASIBLE would still need an independently checkable certificate for a publication-grade negative result.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    if valid:(out/'WITNESS.txt').write_text(' '.join(map(str,candidate))+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
