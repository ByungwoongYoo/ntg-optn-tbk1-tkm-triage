#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, itertools, json, time
from pathlib import Path
import z3

def D(v, k):
    return v < k

def verify(table, k):
    n=len(table)
    des=lambda v:v<k
    for x in range(n):
        for y in range(n):
            if des(x) and des(table[x][y]) and not des(y):
                return False, {'failure':'modus_ponens','x':x,'y':y}
    for x,y,z,u in itertools.product(range(n), repeat=4):
        left=table[table[x][y]][z]
        right=table[table[y][table[z][u]]][table[y][u]]
        value=table[left][right]
        if not des(value):
            return False, {'failure':'u4','valuation':[x,y,z,u],'value':value}
    bad=[x for x in range(n) if not des(table[x][x])]
    if not bad:
        return False, {'failure':'reflexivity_not_refuted'}
    return True, {'reflexivity_failures':bad}

def solve(n, k, timeout_ms):
    imp=z3.Function(f'imp_{n}_{k}', z3.IntSort(), z3.IntSort(), z3.IntSort())
    s=z3.Solver(); s.set(timeout=timeout_ms)
    for a in range(n):
        for b in range(n):
            s.add(imp(a,b)>=0, imp(a,b)<n)
    for x in range(n):
        for y in range(n):
            s.add(z3.Implies(z3.And(D(x,k), D(imp(x,y),k)), D(y,k)))
    for x,y,z,u in itertools.product(range(n), repeat=4):
        left=imp(imp(x,y),z)
        right=imp(imp(y,imp(z,u)),imp(y,u))
        s.add(D(imp(left,right),k))
    s.add(z3.Or([z3.Not(D(imp(x,x),k)) for x in range(n)]))
    t=time.time(); status=s.check(); elapsed=time.time()-t
    rec={'n':n,'designated_count':k,'status':str(status).upper(),'elapsed_seconds':elapsed,
         'z3_version':z3.get_version_string(),'timeout_ms':timeout_ms}
    if status==z3.sat:
        m=s.model(); table=[[m.eval(imp(a,b),model_completion=True).as_long() for b in range(n)] for a in range(n)]
        ok,detail=verify(table,k); rec.update({'implication_table':table,'independent_verifier_ok':ok,'verification_detail':detail})
    elif status==z3.unknown:
        rec['reason_unknown']=s.reason_unknown()
    return rec

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,required=True); ap.add_argument('--seconds',type=int,default=600); ap.add_argument('--out-dir',required=True)
    a=ap.parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    result={'problem':'Ulrich u4 finite algebra countermodel search','n':a.n,'runs':[]}
    for k in range(1,a.n):
        r=solve(a.n,k,a.seconds*1000); result['runs'].append(r); print(json.dumps(r),flush=True)
        if r.get('independent_verifier_ok'):
            break
    result['conclusion']='FINITE_COUNTERMODEL_FOUND_AND_VERIFIED' if any(r.get('independent_verifier_ok') for r in result['runs']) else 'NO_VERIFIED_COUNTERMODEL_IN_COMPLETED_BOUNDED_SEARCH'
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (out/'SCRIPT_SHA256.txt').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()+'  '+Path(__file__).name+'\n')
if __name__=='__main__': main()
