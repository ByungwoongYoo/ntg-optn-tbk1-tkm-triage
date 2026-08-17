#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,os,time
from pathlib import Path
from ortools.sat.python import cp_model
N=100;K=14

def verify(A):
 c={t:0 for t in range(1,N)}
 for a in A:
  for b in A:
   if a!=b:c[(a-b)%N]+=1
 return len(A)==K and len(set(A))==K and max(c.values(),default=0)<=2,c

def solve(g,limit):
 m=cp_model.CpModel()
 # Cyclic gaps d0...d13, with d0=g after rotating one minimum gap to the front.
 d=[m.NewIntVar(g,N,f'd{i}') for i in range(K)];m.Add(d[0]==g);m.Add(sum(d)==N)
 # reflection orbit breaker: every cycle with d0 fixed has an orientation with d1 <= d13.
 m.Add(d[1] <= d[13])
 p=[m.NewIntVar(0,N-1,f'p{i}') for i in range(K)];m.Add(p[0]==0)
 for i in range(1,K):m.Add(p[i]==sum(d[:i]))
 # For every unordered pair, circular distance class 1..50. Each class t<50 corresponds to r_A(t)=r_A(100-t), and may occur <=2 times; t=50 may occur <=1 because both orientations contribute to r_A(50).
 buckets={t:[] for t in range(1,51)}
 for i in range(K):
  for j in range(i+1,K):
   delta=m.NewIntVar(1,99,f'delta_{i}_{j}');m.Add(delta==p[j]-p[i])
   comp=m.NewIntVar(1,99,f'comp_{i}_{j}');m.Add(comp==N-delta)
   dist=m.NewIntVar(1,50,f'dist_{i}_{j}');m.AddMinEquality(dist,[delta,comp])
   bs=[]
   for t in range(1,51):
    b=m.NewBoolVar(f'is_{i}_{j}_{t}');m.Add(dist==t).OnlyEnforceIf(b);m.Add(dist!=t).OnlyEnforceIf(b.Not());buckets[t].append(b);bs.append(b)
   m.AddExactlyOne(bs)
 for t in range(1,50):m.Add(sum(buckets[t])<=2)
 m.Add(sum(buckets[50])<=1)
 s=cp_model.CpSolver();s.parameters.max_time_in_seconds=limit;s.parameters.num_search_workers=min(16,os.cpu_count() or 4);s.parameters.random_seed=20260831+g;s.parameters.log_search_progress=True
 st=s.Solve(m);name=s.StatusName(st);A=[];G=[]
 if st in (cp_model.OPTIMAL,cp_model.FEASIBLE):
  G=[s.Value(x) for x in d];A=[s.Value(x) for x in p]
 ok,c=verify(A) if A else (False,{})
 return name,A,G,ok,c,s.ResponseStats()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--gap',type=int,required=True);ap.add_argument('--limit',type=int,default=1500);ap.add_argument('--out-dir',required=True);a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
 t=time.time();st,A,G,ok,c,stats=solve(a.gap,a.limit)
 res={'problem':'14-element B_2[2] subset of Z_100','encoding':'cyclic gap vector + circular pair-distance classes','gap':a.gap,'status':st,'candidate':A,'gaps':G,'verified_witness':ok,'difference_counts':c if ok else {},'wall_seconds':time.time()-t,'solver_stats':stats,
 'coverage':'Combined with g=1..7, minimum-gap normalization is exhaustive because 14 positive cyclic gaps sum to 100, so min gap <=7.',
 'claim_boundary':'A verified FEASIBLE candidate is a decisive construction. INFEASIBLE closes this normalized gap case computationally; global nonexistence would require all g=1..7 plus independently checkable proof artifacts.'}
 (out/'RESULT.json').write_text(json.dumps(res,indent=2),encoding='utf-8');
 if ok:(out/'WITNESS.txt').write_text(' '.join(map(str,A))+'\n',encoding='utf-8')
 print(json.dumps(res,indent=2))
if __name__=='__main__':main()
