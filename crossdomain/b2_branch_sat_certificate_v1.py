#!/usr/bin/env python3
"""Certificate CNF generator for normalized B2[2] branches in Z/100Z.

The formula uses x_i for membership and y_{i,j} <-> (x_i and x_j).
For circular distances 1..49 at most two selected unordered pairs are
allowed. For distance 50 at most one pair is allowed, because each such
pair contributes twice to the single ordered difference 50.

Normalization:
  * g is the minimum cyclic gap, rotated to the selected pair (0,g);
  * t is the third selected element in increasing order;
  * reflection is fixed by h=t-g <= wrap gap, hence max(A)<=100-h;
  * optional u is the fourth selected element.

The quotient-pair cuts are redundant consequences of the original ordered
B2[2] constraints. They are included only to strengthen SAT propagation.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool

N=100
K=14
MODS=(2,4,5,10,20,25,50)

def build(g:int,t:int|None=None,u:int|None=None):
    if not 1 <= g <= 7:
        raise ValueError('g must be 1..7')
    if t is not None:
        if not (2*g <= t <= 50-5*g):
            raise ValueError(f't must be in [{2*g},{50-5*g}] for g={g}')
    if u is not None:
        if t is None:
            raise ValueError('u requires t')
        umin=t+g
        umax=100-(t-g)-10*g
        if not (umin <= u <= umax):
            raise ValueError(f'u must be in [{umin},{umax}]')

    pool=IDPool(start_from=1)
    x=[pool.id(f'x_{i}') for i in range(N)]
    cnf=CNF()
    cnf.extend(CardEnc.equals(lits=x,bound=K,vpool=pool,encoding=EncType.cardnetwrk).clauses)

    # Rotate a minimum cyclic gap to (0,g), and prohibit every shorter cyclic gap.
    cnf.append([x[0]])
    cnf.append([x[g]])
    for i in range(N):
        for d in range(1,g):
            j=(i+d)%N
            cnf.append([-x[i],-x[j]])

    # Fix the third selected point and reflection orientation.
    if t is not None:
        cnf.append([x[t]])
        for i in range(g+1,t):
            cnf.append([-x[i]])
        h=t-g
        # wrap=100-max(A) >= h, so no selected point may exceed 100-h.
        for i in range(101-h,N):
            cnf.append([-x[i]])

    # Optional fourth-point split.
    if u is not None:
        cnf.append([x[u]])
        for i in range(t+1,u):
            cnf.append([-x[i]])

    # A valid 14-set must contain 5..9 even elements: otherwise the ordered
    # even-difference capacity (49 residues x 2) is exceeded.
    evens=[x[i] for i in range(0,N,2)]
    cnf.extend(CardEnc.atleast(lits=evens,bound=5,vpool=pool,encoding=EncType.seqcounter).clauses)
    cnf.extend(CardEnc.atmost(lits=evens,bound=9,vpool=pool,encoding=EncType.seqcounter).clauses)

    buckets={d:[] for d in range(1,51)}
    pair_var={}
    for i in range(N):
        for j in range(i+1,N):
            y=pool.id(f'y_{i}_{j}')
            pair_var[(i,j)]=y
            # y iff x_i and x_j.
            cnf.append([-y,x[i]])
            cnf.append([-y,x[j]])
            cnf.append([y,-x[i],-x[j]])
            d=min(j-i,N-(j-i))
            buckets[d].append(y)

    for d in range(1,50):
        cnf.extend(CardEnc.atmost(lits=buckets[d],bound=2,vpool=pool,encoding=EncType.seqcounter).clauses)
    cnf.extend(CardEnc.atmost(lits=buckets[50],bound=1,vpool=pool,encoding=EncType.seqcounter).clauses)

    # Redundant quotient-capacity cuts. For modulus m|100, selected pairs
    # in the same residue class contribute twice to quotient difference 0,
    # whose original nonzero preimages have total capacity 2*(100/m-1).
    # For a non-self-inverse quotient distance q, each selected unordered
    # pair contributes once to q and once to -q; capacity is 2*(100/m).
    # At q=m/2, each pair contributes twice, so at most 100/m pairs.
    for m in MODS:
        qb={q:[] for q in range(0,m//2+1)}
        for (i,j),y in pair_var.items():
            r=(j-i)%m
            q=min(r,(m-r)%m)
            qb[q].append(y)
        cnf.extend(CardEnc.atmost(lits=qb[0],bound=100//m-1,vpool=pool,encoding=EncType.seqcounter).clauses)
        for q in range(1,m//2+1):
            if m%2==0 and q==m//2:
                bound=100//m
            else:
                bound=2*(100//m)
            if len(qb[q])>bound:
                cnf.extend(CardEnc.atmost(lits=qb[q],bound=bound,vpool=pool,encoding=EncType.seqcounter).clauses)

    meta={
      'N':N,'K':K,'gap':g,'third':t,'fourth':u,
      'normalization':{
        'minimum_gap_pair':[0,g],
        'third_range':[2*g,50-5*g],
        'reflection_condition':None if t is None else f'max(A) <= {100-(t-g)}',
        'fourth_range':None if t is None else [t+g,100-(t-g)-10*g],
      },
      'variables':pool.top,'clauses':len(cnf.clauses),
      'encoding':'PySAT cardinality network for exact-14; sequential counters for at-most/at-least constraints; Tseitin pair products',
      'claim_boundary':'UNSAT plus a checked DRAT/LRAT proof excludes only this normalized branch; global nonexistence requires complete branch coverage.'
    }
    return cnf,meta

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gap',type=int,required=True)
    ap.add_argument('--third',type=int)
    ap.add_argument('--fourth',type=int)
    ap.add_argument('--out-dir',required=True)
    a=ap.parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    st=time.time(); cnf,meta=build(a.gap,a.third,a.fourth); meta['build_seconds']=time.time()-st
    stem=f'g{a.gap}'
    if a.third is not None: stem+=f'_t{a.third}'
    if a.fourth is not None: stem+=f'_u{a.fourth}'
    cnf_path=out/f'{stem}.cnf'; cnf.to_file(cnf_path)
    meta['cnf_file']=cnf_path.name
    (out/f'{stem}.metadata.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(meta,indent=2))
if __name__=='__main__': main()
