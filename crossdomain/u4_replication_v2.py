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
    return {'state':state,'returncode':rc,'seconds':time.time()-t,'szs_statuses':statuses,'tail':text[-4000:]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--vampire',required=True);ap.add_argument('--u1',required=True);ap.add_argument('--u2',required=True);ap.add_argument('--u3',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    originals={n:Path(getattr(a,n)).read_text(encoding='utf-8') for n in ('u1','u2','u3')}
    # Use u1 artifact as template: replace its u1 axiom line only; keep exact mp, refl, and any comments/directives.
    templ=originals['u1']
    lines=templ.splitlines(); replaced=False;new=[]
    for line in lines:
        if re.match(r'\s*fof\(u1\s*,\s*axiom\s*,',line):new.append(U4);replaced=True
        else:new.append(line)
    if not replaced: raise RuntimeError('Could not locate u1 axiom in official artifact')
    u4p=out/'u4.p';u4p.write_text('\n'.join(new)+'\n',encoding='utf-8')
    result={'vampire_version':run([a.vampire,'--version'],20,out/'version.log'),'controls':{},'u4':{}}
    # Exact controls using same mode as paper; generous but bounded.
    for n in ('u1','u2','u3'):
        result['controls'][n]=run([a.vampire,'--mode','casc_sat','--time_limit','240',getattr(a,n)],270,out/f'{n}.log')
    # Frozen attempts on u4: exact paper mode and a regular CASC portfolio; no retuning from output.
    result['u4']['casc_sat_900']=run([a.vampire,'--mode','casc_sat','--time_limit','900',str(u4p)],930,out/'u4_casc_sat.log')
    result['u4']['casc_600']=run([a.vampire,'--mode','casc','--time_limit','600',str(u4p)],630,out/'u4_casc.log')
    decisive=[]
    for name,r in result['u4'].items():
        for s in r['szs_statuses']:
            if s.lower() in ('satisfiable','countersatisfiable','unsatisfiable','theorem'):decisive.append({'run':name,'status':s})
    result['decisive_candidate']=decisive
    result['claim_boundary']='A Satisfiable/CounterSatisfiable u4 result is only a candidate resolution until its saturated clause set/model is independently checked (e.g. E/iProver/self-contained model). Timeout/Unknown is no result.'
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
