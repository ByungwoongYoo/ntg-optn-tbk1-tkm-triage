#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, time
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Cadical195

N=100; K=14

def verify(A):
    counts={t:0 for t in range(1,N)}
    for a in A:
        for b in A:
            if a!=b: counts[(a-b)%N]+=1
    return len(A)==K and len(set(A))==K and max(counts.values(), default=0)<=2, counts

def build(g):
    pool=IDPool(start_from=1)
    x=[pool.id(f'x_{i}') for i in range(N)]
    cnf=CNF()
    cnf.extend(CardEnc.equals(lits=x, bound=K, vpool=pool, encoding=EncType.cardnetwrk).clauses)
    # Every 14-subset has a cyclic minimum gap <= floor(100/14)=7. Rotate one minimum gap to (0,g).
    cnf.append([x[0]]); cnf.append([x[g]])
    # Enforce minimum circular gap >= g; x0,xg then makes the minimum exactly g.
    for i in range(N):
        for d in range(1,g):
            j=(i+d)%N
            if i<j: cnf.append([-x[i],-x[j]])
            elif j<i: cnf.append([-x[i],-x[j]])
    buckets={d:[] for d in range(1,51)}
    for i in range(N):
        for j in range(i+1,N):
            y=pool.id(f'y_{i}_{j}')
            cnf.append([-y,x[i]]); cnf.append([-y,x[j]]); cnf.append([y,-x[i],-x[j]])
            d=min(j-i, N-(j-i)); buckets[d].append(y)
    for d in range(1,50):
        cnf.extend(CardEnc.atmost(lits=buckets[d], bound=2, vpool=pool, encoding=EncType.seqcounter).clauses)
    cnf.extend(CardEnc.atmost(lits=buckets[50], bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
    return cnf,x,pool.top

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--gap',type=int,required=True); ap.add_argument('--out-dir',required=True)
    a=ap.parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    t=time.time(); cnf,x,nvars=build(a.gap); build_s=time.time()-t
    dimacs=out/f'b2_gap_{a.gap}.cnf'; cnf.to_file(dimacs)
    solve_t=time.time();
    with Cadical195(bootstrap_with=cnf.clauses) as s:
        sat=s.solve(); model=s.get_model() if sat else None; stats=s.accum_stats()
    solve_s=time.time()-solve_t
    A=[]
    if sat:
        pos=set(v for v in model if v>0)
        A=[i for i,v in enumerate(x) if v in pos]
    ok,counts=verify(A) if A else (False,{})
    result={'problem':'14-element B_2[2] subset of Z_100','encoding':'PySAT CNF with exact-14, Tseitin pair products, distance cardinalities, exact minimum-gap normalization','gap':a.gap,'sat':bool(sat),'candidate':A,'verified_witness':ok,'difference_counts':counts if ok else {},'cnf_vars':nvars,'cnf_clauses':len(cnf.clauses),'build_seconds':build_s,'solve_seconds':solve_s,'solver':'CaDiCaL 1.9.5 via PySAT','solver_stats':stats,
    'coverage':'Cases g=1..7 are exhaustive because 14 positive cyclic gaps sum to 100, hence the minimum gap is at most 7; any set can be rotated so one minimum gap is (0,g).',
    'claim_boundary':'SAT plus verified candidate is a decisive positive witness. UNSAT closes this normalized case computationally; publication-grade global nonexistence still requires all seven cases and independently checkable proof logging/certificates.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    if ok:(out/'WITNESS.txt').write_text(' '.join(map(str,A))+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__': main()
