#!/usr/bin/env python3
"""Independent verifier for vesuvius_fiber_trace_benchmark.py outputs."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def metrics(pred, gt, tol=2.0):
    p=np.asarray(pred,dtype=int);g=np.asarray(gt,dtype=int)
    if len(p)==0 or len(g)==0:return {'precision':0.0,'recall':0.0,'f1':0.0}
    dpg=cKDTree(g).query(p,k=1)[0];dgp=cKDTree(p).query(g,k=1)[0]
    pr=float(np.mean(dpg<=tol));rc=float(np.mean(dgp<=tol));f=0.0 if pr+rc==0 else 2*pr*rc/(pr+rc)
    return {'precision':pr,'recall':rc,'f1':f,'mean_pred_to_gt':float(np.mean(dpg)),'mean_gt_to_pred':float(np.mean(dgp)),'p95_pred_to_gt':float(np.percentile(dpg,95))}

def close(a,b,eps=1e-9):return abs(float(a)-float(b))<=eps*max(1,abs(float(a)),abs(float(b)))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('result');ap.add_argument('records');ap.add_argument('out');args=ap.parse_args()
    summary=json.loads(Path(args.result).read_text());records=json.loads(Path(args.records).read_text())
    proto=summary['protocol'];train=set(proto['train_cubes']);val=set(proto['validation_cubes']);test=set(proto['test_cubes'])
    assert train and val and test and not (train&val or train&test or val&test)
    assert all(str(r['cube']) in test for r in records)
    failures=[]
    for i,r in enumerate(records):
        for name,key in [('model','model_path'),('intensity','intensity_path'),('straight','straight_path')]:
            m=metrics(r[key],r['gt'])
            for k,v in m.items():
                if k in r[name] and not close(v,r[name][k],1e-8):failures.append({'index':i,'method':name,'metric':k,'stored':r[name][k],'recomputed':v})
        for key in ['model_path','intensity_path','straight_path']:
            p=np.asarray(r[key],dtype=int)
            if len(p):
                if not np.array_equal(p[0],np.asarray(r['start'])) or not np.array_equal(p[-1],np.asarray(r['end'])):
                    failures.append({'index':i,'path':key,'error':'endpoint_mismatch'})
                jumps=np.linalg.norm(np.diff(p,axis=0),axis=1)
                if np.any(jumps>math.sqrt(3)+1e-9):failures.append({'index':i,'path':key,'error':'nonlocal_jump','max':float(jumps.max())})
    out={'verified':not failures,'record_count':len(records),'split_disjoint':True,'all_records_in_test':all(str(r['cube']) in test for r in records),'failures':failures,
         'result_sha256':hashlib.sha256(Path(args.result).read_bytes()).hexdigest(),'records_sha256':hashlib.sha256(Path(args.records).read_bytes()).hexdigest()}
    Path(args.out).write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
