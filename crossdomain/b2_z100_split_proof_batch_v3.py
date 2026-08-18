#!/usr/bin/env python3
"""Adaptively split and certify a canonical B2[2]/Z100 prefix branch.

For a parent prefix [0,1,t,...], every admissible next selected residue is
run as a disjoint child CNF. UNSAT is accepted only after drat-trim reports
VERIFIED; SAT is accepted only after the independent model checker passes.
Timeouts are never evidence.
"""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, pathlib, subprocess, time

K=14

def sha(path:pathlib.Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

def validate_parent(prefix:list[int]):
    if len(prefix)<3 or len(prefix)>=K:
        raise ValueError('parent prefix length must be 3..13')
    if prefix[:2]!=[0,1]:
        raise ValueError('prefix must start [0,1]')
    if any(type(x) is not int for x in prefix):
        raise ValueError('prefix entries must be integers')
    if any(prefix[i]>=prefix[i+1] for i in range(len(prefix)-1)):
        raise ValueError('prefix must be strictly increasing')
    t=prefix[2]
    if not 2<=t<=45:
        raise ValueError('third must be 2..45')
    max_selected=101-t
    if prefix[-1]>max_selected-(K-len(prefix)):
        raise ValueError('parent prefix cannot be completed under reflection bound')
    lo=prefix[-1]+1
    hi=max_selected-(K-(len(prefix)+1))
    if lo>hi:
        raise ValueError('no admissible child range')
    return t,max_selected,list(range(lo,hi+1))

def run_cmd(cmd,timeout,stdout_path,stderr_path):
    with stdout_path.open('wb') as out, stderr_path.open('wb') as err:
        try:
            p=subprocess.run(cmd,stdout=out,stderr=err,timeout=timeout,check=False)
            return p.returncode,False
        except subprocess.TimeoutExpired:
            return 124,True

def worker(value,args,prefix,outroot):
    child=prefix+[value]
    cid='p'+'-'.join(map(str,child))
    d=outroot/cid
    d.mkdir(parents=True,exist_ok=True)
    cnf=d/f'{cid}.cnf'
    meta=d/f'{cid}.metadata.json'
    proof=d/f'{cid}.drat'
    model=d/f'{cid}.model'
    rec={'case_id':cid,'parent_prefix':prefix,'child_prefix':child,'split_value':value,'status':'INITIALIZING'}
    st=time.time()
    gcmd=[args.python,args.generator,'--prefix-json',json.dumps(child,separators=(',',':')),'--output',str(cnf),'--metadata',str(meta)]
    grc,_=run_cmd(gcmd,120,d/'generator_stdout.log',d/'generator_stderr.log')
    if grc!=0:
        rec['status']=f'GENERATOR_ERROR_{grc}'
        rec['wall_seconds']=time.time()-st
        return rec
    md=json.loads(meta.read_text())
    if md.get('prefix')!=child or md.get('cnf_sha256')!=sha(cnf):
        rec['status']='GENERATOR_METADATA_MISMATCH'
        rec['wall_seconds']=time.time()-st
        return rec
    rec.update({'cnf_sha256':md['cnf_sha256'],'variables':md['variables'],'clauses':md['clauses']})
    src,did_timeout=run_cmd([args.cadical,'-w',str(model),str(cnf),str(proof)],args.solve_seconds,d/'cadical.log',d/'cadical_stderr.log')
    (d/'cadical_exit_code.txt').write_text(str(src)+'\n')
    if src==20:
        prc,p_to=run_cmd([args.drat_trim,str(cnf),str(proof)],args.check_seconds,d/'drat_check.log',d/'drat_check_stderr.log')
        (d/'drat_check_exit_code.txt').write_text(str(prc)+'\n')
        log=(d/'drat_check.log').read_text(errors='replace') if (d/'drat_check.log').exists() else ''
        if prc==0 and 'VERIFIED' in log:
            rec['status']='VERIFIED_UNSAT'
        elif p_to or prc in (124,143):
            rec['status']='UNSAT_PROOF_CHECK_TIMEOUT_NOT_VERIFIED'
        else:
            rec['status']=f'UNSAT_PROOF_CHECK_FAILED_{prc}'
    elif src==10:
        vrc,_=run_cmd([args.python,args.model_checker,'--cnf',str(cnf),'--model',str(model),'--output',str(d/'model_verification.json')],120,d/'model_verification.log',d/'model_verification_stderr.log')
        rec['status']='VERIFIED_SAT_COUNTEREXAMPLE' if vrc==0 else 'SAT_MODEL_CHECK_FAILED'
    elif did_timeout or src in (124,143):
        rec['status']='TIMEOUT_NOT_EVIDENCE'
    else:
        rec['status']=f'SOLVER_ERROR_{src}'
    if proof.exists() and proof.stat().st_size:
        rec['drat_raw_bytes']=proof.stat().st_size
        rec['drat_raw_sha256']=sha(proof)
        (d/'DRAT_RAW_BYTES.txt').write_text(str(proof.stat().st_size)+'\n')
        (d/'DRAT_RAW_SHA256.txt').write_text(rec['drat_raw_sha256']+'  '+proof.name+'\n')
        if rec['status']=='VERIFIED_UNSAT':
            zrc,_=run_cmd([args.zstd,'-T0','-8','--rm',str(proof)],600,d/'zstd.log',d/'zstd_stderr.log')
            if zrc!=0:
                rec['status']=f'PROOF_COMPRESSION_ERROR_{zrc}'
            else:
                z=pathlib.Path(str(proof)+'.zst')
                rec['drat_zst_bytes']=z.stat().st_size
                rec['drat_zst_sha256']=sha(z)
        else:
            proof.unlink(missing_ok=True)
    rec['wall_seconds']=time.time()-st
    (d/'RESULT.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
    lines=[]
    for p in sorted(d.iterdir()):
        if p.is_file() and p.name!='SHA256SUMS.txt':
            lines.append(f'{sha(p)}  {p.name}')
    (d/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')
    return rec

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--prefix-json',required=True)
    ap.add_argument('--generator',required=True)
    ap.add_argument('--model-checker',required=True)
    ap.add_argument('--cadical',default='cadical')
    ap.add_argument('--drat-trim',default='drat-trim')
    ap.add_argument('--zstd',default='zstd')
    ap.add_argument('--python',default='python3')
    ap.add_argument('--out',type=pathlib.Path,required=True)
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--solve-seconds',type=int,default=300)
    ap.add_argument('--check-seconds',type=int,default=120)
    ap.add_argument('--plan-only',action='store_true')
    args=ap.parse_args()
    prefix=json.loads(args.prefix_json)
    t,max_selected,values=validate_parent(prefix)
    args.out.mkdir(parents=True,exist_ok=True)
    plan={
      'schema':'b2-z100-complete-next-selected-split-v3',
      'parent_prefix':prefix,
      'third':t,
      'reflection_max_selected':max_selected,
      'expected_values':values,
      'expected_child_count':len(values),
      'solve_seconds':args.solve_seconds,
      'check_seconds':args.check_seconds
    }
    (args.out/'SPLIT_PLAN.json').write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
    if args.plan_only:
        print(json.dumps(plan,indent=2,sort_keys=True))
        return
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1,args.workers)) as ex:
        futures={ex.submit(worker,v,args,prefix,args.out):v for v in values}
        for f in concurrent.futures.as_completed(futures):
            try:
                r=f.result()
            except Exception as e:
                v=futures[f]
                r={'split_value':v,'child_prefix':prefix+[v],'status':'WORKER_EXCEPTION','error':repr(e)}
            rows.append(r)
            print(json.dumps(r,sort_keys=True),flush=True)
    rows.sort(key=lambda r:r['split_value'])
    observed=[r['split_value'] for r in rows]
    duplicate=len(observed)!=len(set(observed))
    missing=sorted(set(values)-set(observed))
    extra=sorted(set(observed)-set(values))
    counts={}
    for r in rows:
        counts[r['status']]=counts.get(r['status'],0)+1
    unresolved=[r['child_prefix'] for r in rows if r['status']!='VERIFIED_UNSAT']
    counterexamples=[r for r in rows if r['status']=='VERIFIED_SAT_COUNTEREXAMPLE']
    complete=(not duplicate and not missing and not extra and not unresolved)
    summary={
      'schema':'b2-z100-split-proof-summary-v3',
      'parent_prefix':prefix,
      'expected_child_count':len(values),
      'observed_child_count':len(rows),
      'missing_values':missing,
      'extra_values':extra,
      'duplicate_values':duplicate,
      'counts':counts,
      'unresolved_child_prefixes':unresolved,
      'verified_sat_counterexamples':counterexamples,
      'parent_verified_unsat_by_complete_children':complete,
      'claim_boundary':'The parent branch is excluded only if every disjoint next-selected child is VERIFIED_UNSAT and coverage has no missing/extra/duplicate child.'
    }
    (args.out/'CHILDREN.json').write_text(json.dumps(rows,indent=2,sort_keys=True)+'\n')
    (args.out/'SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (args.out/'UNRESOLVED_CHILDREN.json').write_text(json.dumps(unresolved,indent=2)+'\n')
    status='VERIFIED_UNSAT_BY_COMPLETE_CHILDREN' if complete else ('VERIFIED_SAT_COUNTEREXAMPLE' if counterexamples else 'PARTIAL_UNRESOLVED')
    (args.out/'STATUS.txt').write_text(status+'\n')
    print(json.dumps(summary,indent=2,sort_keys=True))
    if counterexamples:
        raise SystemExit(10)
    if not complete:
        raise SystemExit(2)

if __name__=='__main__':
    main()
