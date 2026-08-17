#!/usr/bin/env python3
from __future__ import annotations
import json, os, shutil, subprocess, time
from pathlib import Path

FORMULAS={
'u1':'fof(u1,axiom, ![X,Y,Z,U]: p(i(i(i(X,Y),Z),i(i(Z,i(Z,U)),i(Y,U))))).',
'u2':'fof(u2,axiom, ![X,Y,Z,U]: p(i(i(X,i(Y,Z)),i(i(i(U,X),Y),i(X,Z))))).',
'u3':'fof(u3,axiom, ![X,Y,Z,U]: p(i(i(X,Y),i(i(i(Z,X),i(Y,U)),i(X,U))))).',
'u4':'fof(u4,axiom, ![X,Y,Z,U]: p(i(i(i(X,Y),Z),i(i(Y,i(Z,U)),i(Y,U))))).',
}
MP='fof(mp,axiom, ![X,Y]: ((p(i(X,Y)) & p(X)) => p(Y))).'
REFL='fof(refl,conjecture, ![X]: p(i(X,X))).'

def run(cmd,timeout,log):
    t=time.time()
    try:
        p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=timeout)
        text=p.stdout; rc=p.returncode; status='completed'
    except subprocess.TimeoutExpired as e:
        text=(e.stdout or '') if isinstance(e.stdout,str) else ((e.stdout or b'').decode(errors='replace'))
        rc=None; status='timeout'
    Path(log).write_text(text,encoding='utf-8')
    return {'status':status,'returncode':rc,'seconds':time.time()-t,'tail':text[-6000:]}

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/u4_probe')); out.mkdir(parents=True,exist_ok=True)
    for name,formula in FORMULAS.items():
        (out/f'{name}.p').write_text('\n'.join([MP,formula,REFL])+'\n',encoding='utf-8')
    result={'vampire':shutil.which('vampire'),'eprover':shutil.which('eprover'),'tests':{}}
    if result['vampire']:
        for name in ['u1','u2','u3']:
            result['tests'][f'vampire_{name}_casc_sat']=run(['vampire','--mode','casc_sat','--time_limit','90',str(out/f'{name}.p')],110,out/f'vampire_{name}.log')
        # Prespecified u4 attempts: same published mode and generic portfolio. No parameter selection after result.
        result['tests']['vampire_u4_casc_sat']=run(['vampire','--mode','casc_sat','--time_limit','600',str(out/'u4.p')],630,out/'vampire_u4_casc_sat.log')
        result['tests']['vampire_u4_casc']=run(['vampire','--mode','casc','--time_limit','300',str(out/'u4.p')],330,out/'vampire_u4_casc.log')
    if result['eprover']:
        for name in ['u1','u2','u3','u4']:
            result['tests'][f'eprover_{name}']=run(['eprover','--auto','--satauto','--cpu-limit=180',str(out/f'{name}.p')],210,out/f'eprover_{name}.log')
    # Detect any decisive standard TPTP status.
    decisive={}
    for k,v in result['tests'].items():
        tail=v.get('tail','')
        statuses=[]
        for token in ['Satisfiable','CounterSatisfiable','Unsatisfiable','Theorem','ContradictoryAxioms']:
            if token in tail: statuses.append(token)
        decisive[k]=statuses
    result['detected_statuses']=decisive
    u4_status=sum([decisive.get(k,[]) for k in decisive if 'u4' in k],[])
    result['u4_decisive_candidate']=bool(u4_status)
    result['interpretation']='Any u4 decisive status requires independent proof/model verification before claiming resolution.'
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
