#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,time
from collections import Counter
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

N=100; K=14

def cd(i,j):
    d=abs(i-j)%N
    return min(d,N-d)

def verify(S):
    c=Counter((a-b)%N for a in S for b in S if a!=b)
    bad={d:n for d,n in c.items() if n>2}
    return (len(S)==K and len(set(S))==K and not bad),bad,max(c.values(),default=0)

def build(m):
    v=IDPool(); x=[v.id(f'x{i}') for i in range(N)]; cnf=CNF()
    cnf.extend(CardEnc.equals(x,K,vpool=v,encoding=EncType.seqcounter).clauses)
    cnf.append([x[0]]); cnf.append([x[m]])
    for i in range(1,m): cnf.append([-x[i]])
    if m>1:
        for i in range(N):
            for j in range(i+1,N):
                if cd(i,j)<m: cnf.append([-x[i],-x[j]])
    yd={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            y=v.id(f'y{i}_{j}')
            cnf.extend([[-y,x[i]],[-y,x[j]],[-x[i],-x[j],y]])
            yd[cd(i,j)].append(y)
    for d in range(1,50): cnf.extend(CardEnc.atmost(yd[d],2,vpool=v,encoding=EncType.seqcounter).clauses)
    cnf.extend(CardEnc.atmost(yd[50],1,vpool=v,encoding=EncType.seqcounter).clauses)
    return cnf,x,v

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--min-gap',type=int,required=True,choices=range(1,8)); ap.add_argument('--out-dir',required=True); ap.add_argument('--solver',default='cadical195'); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    cnf,x,v=build(a.min_gap)
    (out/'instance.cnf').write_text(cnf.to_dimacs(),encoding='ascii')
    t=time.time()
    status=None; model=None; solver_used=a.solver; err=None
    try:
        with Solver(name=a.solver,bootstrap_with=cnf.clauses) as s:
            sat=s.solve(); model=s.get_model() if sat else None; status='SAT' if sat else 'UNSAT'
    except Exception as e:
        err=repr(e); solver_used='glucose42'
        with Solver(name=solver_used,bootstrap_with=cnf.clauses) as s:
            sat=s.solve(); model=s.get_model() if sat else None; status='SAT' if sat else 'UNSAT'
    elapsed=time.time()-t
    S=[]; check=None
    if model:
        pos=set(z for z in model if z>0)
        S=[i for i,z in enumerate(x) if z in pos]
        ok,bad,mx=verify(S); check={'verified':ok,'violations':bad,'max_ordered_difference_multiplicity':mx}
    res={'encoding':'Independent CNF/PySAT sequential-cardinality encoding','min_gap_case':a.min_gap,'solver':solver_used,'fallback_error':err,'status':status,'elapsed_seconds':elapsed,'nvars':cnf.nv,'nclauses':len(cnf.clauses),'solution':S,'verification':check,'claim_boundary':'SAT plus independent difference verification gives a witness. UNSAT independently reproduces infeasibility for this symmetry case. A global resolution requires all seven min-gap cases and external review/certificate checking.'}
    (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8'); print(json.dumps(res,indent=2),flush=True)
if __name__=='__main__': main()
