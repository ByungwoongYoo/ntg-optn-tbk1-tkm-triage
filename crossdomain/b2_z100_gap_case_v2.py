#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,os
from pathlib import Path
from ortools.sat.python import cp_model
N=100;K=14

def verify(A):
 c={t:0 for t in range(1,N)}
 for a in A:
  for b in A:
   if a!=b:c[(a-b)%N]+=1
 return len(A)==K and len(set(A))==K and max(c.values(),default=0)<=2,c

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--gap',type=int,required=True);ap.add_argument('--seconds',type=int,default=900);ap.add_argument('--out-dir',required=True);a=ap.parse_args();g=a.gap
 if not 1<=g<=7:raise ValueError(g)
 out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
 m=cp_model.CpModel();x=[m.NewBoolVar(f'x{i}') for i in range(N)];m.Add(sum(x)==K);m.Add(x[0]==1);m.Add(x[g]==1)
 # By translation, choose one pair realizing the minimum cyclic adjacent gap as 0,g.
 # Enforce that no two selected points lie at circular distance d<g.
 for d in range(1,g):
  for i in range(N):m.Add(x[i]+x[(i+d)%N]<=1)
 # Product variables shared by difference constraints.
 y={}
 for i in range(N):
  for j in range(i+1,N):
   v=m.NewBoolVar(f'y{i}_{j}');y[i,j]=v;m.Add(v<=x[i]);m.Add(v<=x[j]);m.Add(v>=x[i]+x[j]-1)
 for t in range(1,N):
  terms=[]
  for b in range(N):
   aa=(b+t)%N;u,v=sorted((aa,b));terms.append(y[u,v])
  m.Add(sum(terms)<=2)
 # Remaining reflection symmetry while 0,g fixed: reflect about g/2, mapping i -> g-i mod N.
 # Safe lex-breaking is omitted rather than risk excluding an orbit.
 s=cp_model.CpSolver();s.parameters.max_time_in_seconds=a.seconds;s.parameters.num_search_workers=min(16,os.cpu_count() or 4);s.parameters.random_seed=20260817+g;s.parameters.log_search_progress=True
 st=s.Solve(m);name=s.StatusName(st);A=[i for i in range(N) if st in (cp_model.OPTIMAL,cp_model.FEASIBLE) and s.Value(x[i])];ok,c=verify(A) if A else (False,{})
 res={'gap_case':g,'partition_argument':'Every 14-set on a 100-cycle has minimum cyclic adjacent gap g<=floor(100/14)=7; after translation/orientation a minimum-gap pair may be fixed at {0,g}. Thus cases g=1..7 cover all possible 14-sets.','status':name,'candidate':A,'verified':ok,'max_difference_multiplicity':max(c.values()) if c else None,'solver_stats':s.ResponseStats(),'claim_boundary':'FEASIBLE+verified gives a decisive 14-set witness. INFEASIBLE closes only this gap case and is not a publication-grade UNSAT certificate by itself. UNKNOWN closes nothing.'}
 (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
 if ok:(out/'WITNESS.txt').write_text(' '.join(map(str,A))+'\n')
 print(json.dumps(res,indent=2))
if __name__=='__main__':main()
