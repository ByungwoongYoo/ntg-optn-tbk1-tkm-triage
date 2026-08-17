#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re, subprocess, time
from pathlib import Path

U4='fof(u4,axiom, ![X,Y,Z,U]: p(i(i(i(X, Y), Z), i(i(Y, i(Z, U)), i(Y, U))))).'

def run(cmd,timeout,log):
    t=time.time()
    try:
        p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=timeout)
        text=p.stdout;rc=p.returncode;state='completed'
    except subprocess.TimeoutExpired as e:
        text=(e.stdout or b'').decode(errors='replace') if isinstance(e.stdout,(bytes,bytearray)) else str(e.stdout or '')
        rc=None;state='timeout'
    Path(log).write_text(text,encoding='utf-8')
    statuses=re.findall(r'%\s*SZS status\s+(\w+)',text)
    return {'state':state,'returncode':rc,'seconds':time.time()-t,'command':cmd,'szs_statuses':statuses,'tail':text[-5000:]}

def portfolio(vampire,problem,time_limit):
    # Exact command form written inside the authors' Zenodo u1/u2/u3 input files.
    return [vampire,'--mode','portfolio','--schedule','casc_sat','--saturation_algorithm','lrs','--time_limit',str(time_limit),problem]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--vampire',required=True);ap.add_argument('--u1',required=True);ap.add_argument('--u2',required=True);ap.add_argument('--u3',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    originals={n:Path(getattr(a,n)).read_text(encoding='utf-8') for n in ('u1','u2','u3')}
    templ=originals['u1'];lines=templ.splitlines();replaced=False;new=[]
    for line in lines:
        if re.match(r'\s*fof\(u1\s*,\s*axiom\s*,',line):new.append(U4);replaced=True
        else:new.append(line.replace('candidate u1','candidate u4').replace('u1.p','u4.p'))
    if not replaced:raise RuntimeError('Could not locate u1 axiom in official artifact')
    u4p=out/'u4.p';u4p.write_text('\n'.join(new)+'\n',encoding='utf-8')
    result={'vampire_version':run([a.vampire,'--version'],20,out/'version.log'),'authors_documented_command':'vampire --mode portfolio --schedule casc_sat --saturation_algorithm lrs <problem.p>','controls':{},'u4':{}}
    # Positive controls: must reproduce the authors' decisive saturation/countermodel statuses before u4 is interpretable.
    for n in ('u1','u2','u3'):
        result['controls'][n]=run(portfolio(a.vampire,getattr(a,n),300),330,out/f'{n}.log')
    controls_ok=all(any(s.lower() in ('satisfiable','countersatisfiable') for s in r['szs_statuses']) for r in result['controls'].values())
    result['controls_reproduced']=controls_ok
    # Frozen u4 attempts. First is exactly the documented mode with larger budget. Second is generic portfolio sensitivity.
    result['u4']['documented_casc_sat_1200']=run(portfolio(a.vampire,str(u4p),1200),1230,out/'u4_documented.log')
    result['u4']['generic_casc_600']=run([a.vampire,'--mode','casc','--time_limit','600',str(u4p)],630,out/'u4_casc.log')
    decisive=[]
    for name,r in result['u4'].items():
        for s in r['szs_statuses']:
            if s.lower() in ('satisfiable','countersatisfiable','unsatisfiable','theorem'):decisive.append({'run':name,'status':s})
    result['decisive_candidate']=decisive
    result['interpretable_u4_result']=bool(controls_ok and decisive)
    result['claim_boundary']='Only if u1/u2/u3 controls reproduce and u4 yields a decisive status is there a candidate result. Any u4 countermodel/saturated-set claim then requires an independent checker (E/iProver or explicit model). Timeout/Unknown is no result.'
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
