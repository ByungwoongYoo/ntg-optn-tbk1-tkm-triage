#!/usr/bin/env python3
"""Recursively close a canonical B2[2]/Z100 prefix branch with proof objects.

A node is accepted as UNSAT only through one of two independently auditable
routes:
  (1) CaDiCaL returns UNSAT and drat-trim reports VERIFIED for the exact CNF;
  (2) the node is partitioned by the next selected residue into every
      admissible disjoint child, and every child is closed.
At a fully fixed 14-set, a direct ordered-difference checker either rejects
that one set or reports a verified counterexample. Timeouts are never proof.
"""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, pathlib, subprocess, time

N=100
K=14

class WallBudget(Exception):
    pass

def sha(path:pathlib.Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()

def node_id(prefix):
    return 'p'+'-'.join(map(str,prefix))

def validate_prefix(prefix):
    if not 3<=len(prefix)<=14:
        raise ValueError('prefix length must be 3..14')
    if prefix[:2]!=[0,1] or any(type(x) is not int for x in prefix):
        raise ValueError('prefix must be integer list beginning [0,1]')
    if any(prefix[i]>=prefix[i+1] for i in range(len(prefix)-1)):
        raise ValueError('prefix must be strictly increasing')
    t=prefix[2]
    if not 2<=t<=45:
        raise ValueError('third must be 2..45')
    max_selected=101-t
    remaining=K-len(prefix)
    if prefix[-1]>max_selected-remaining:
        raise ValueError('prefix leaves too few available residues')
    return t,max_selected,remaining

def child_values(prefix):
    _,max_selected,remaining=validate_prefix(prefix)
    if remaining==0:
        return []
    lo=prefix[-1]+1
    hi=max_selected-(remaining-1)
    return list(range(lo,hi+1))

def direct_full_check(prefix):
    validate_prefix(prefix)
    assert len(prefix)==K
    counts=[0]*N
    for a in prefix:
        for b in prefix:
            if a!=b:
                counts[(a-b)%N]+=1
    bad={str(d):counts[d] for d in range(1,N) if counts[d]>2}
    return {
      'ordered_nonzero_difference_total':sum(counts[1:]),
      'max_multiplicity':max(counts[1:]),
      'violations':bad,
      'valid_b2_2':not bad
    }

def run_cmd(cmd,timeout,stdout_path,stderr_path):
    with stdout_path.open('wb') as out,stderr_path.open('wb') as err:
        try:
            p=subprocess.run(cmd,stdout=out,stderr=err,timeout=timeout,check=False)
            return p.returncode,False
        except subprocess.TimeoutExpired:
            return 124,True

def write_hashes(d):
    lines=[]
    for p in sorted(d.iterdir()):
        if p.is_file() and p.name!='SHA256SUMS.txt':
            lines.append(f'{sha(p)}  {p.name}')
    (d/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')

class Closer:
    def __init__(self,args):
        self.a=args
        self.start=time.monotonic()
        self.deadline=self.start+args.wall_seconds
        self.nodes={}
        self.counterexamples=[]
        self.unresolved=[]
        self.node_count=0
    def remaining_wall(self):
        return self.deadline-time.monotonic()
    def ensure_wall(self,need=5):
        if self.remaining_wall()<need:
            raise WallBudget
    def solve_seconds(self,depth):
        return max(self.a.min_solve_seconds,int(self.a.root_solve_seconds*(self.a.decay**depth)))
    def check_seconds(self,depth):
        return max(self.a.min_check_seconds,int(self.a.root_check_seconds*(self.a.decay**depth)))
    def attempt_node(self,prefix,depth):
        self.ensure_wall(15)
        nid=node_id(prefix)
        d=self.a.out/'nodes'/nid
        d.mkdir(parents=True,exist_ok=True)
        rec={
          'node_id':nid,'prefix':prefix,'depth_from_root':depth,
          'prefix_length':len(prefix),'attempted_direct_cnf':True,
          'closure':'UNRESOLVED'
        }
        cnf=d/f'{nid}.cnf'; meta=d/f'{nid}.metadata.json'
        proof=d/f'{nid}.drat'; model=d/f'{nid}.model'
        gcmd=[self.a.python,self.a.generator,'--prefix-json',json.dumps(prefix,separators=(',',':')),'--output',str(cnf),'--metadata',str(meta)]
        grc,_=run_cmd(gcmd,min(120,max(1,int(self.remaining_wall()-5))),d/'generator_stdout.log',d/'generator_stderr.log')
        if grc!=0:
            rec['attempt_status']=f'GENERATOR_ERROR_{grc}'
            return rec
        md=json.load(open(meta))
        if md.get('prefix')!=prefix or md.get('cnf_sha256')!=sha(cnf):
            rec['attempt_status']='GENERATOR_METADATA_MISMATCH'
            return rec
        rec.update({'cnf_sha256':md['cnf_sha256'],'variables':md['variables'],'clauses':md['clauses']})
        stime=min(self.solve_seconds(depth),max(1,int(self.remaining_wall()-10)))
        rec['solver_timeout_seconds']=stime
        rc,to=run_cmd([self.a.cadical,'-w',str(model),str(cnf),str(proof)],stime,d/'cadical.log',d/'cadical_stderr.log')
        (d/'cadical_exit_code.txt').write_text(str(rc)+'\n')
        if rc==20:
            ctime=min(self.check_seconds(depth),max(1,int(self.remaining_wall()-5)))
            rec['checker_timeout_seconds']=ctime
            prc,pto=run_cmd([self.a.drat_trim,str(cnf),str(proof)],ctime,d/'drat_check.log',d/'drat_check_stderr.log')
            (d/'drat_check_exit_code.txt').write_text(str(prc)+'\n')
            log=(d/'drat_check.log').read_text(errors='replace') if (d/'drat_check.log').exists() else ''
            if prc==0 and 'VERIFIED' in log:
                rec['attempt_status']='VERIFIED_UNSAT'
                rec['closure']='VERIFIED_UNSAT_LEAF'
            elif pto or prc in (124,143):
                rec['attempt_status']='UNSAT_PROOF_CHECK_TIMEOUT_NOT_VERIFIED'
            else:
                rec['attempt_status']=f'UNSAT_PROOF_CHECK_FAILED_{prc}'
        elif rc==10:
            vrc,_=run_cmd([self.a.python,self.a.model_checker,'--cnf',str(cnf),'--model',str(model),'--output',str(d/'model_verification.json')],min(120,max(1,int(self.remaining_wall()-5))),d/'model_verification.log',d/'model_verification_stderr.log')
            if vrc==0:
                rec['attempt_status']='VERIFIED_SAT_COUNTEREXAMPLE'
                rec['closure']='VERIFIED_SAT_COUNTEREXAMPLE'
                try:
                    rec['counterexample']=json.load(open(d/'model_verification.json'))
                except Exception:
                    rec['counterexample']={'model_file':model.name}
            else:
                rec['attempt_status']='SAT_MODEL_CHECK_FAILED'
        elif to or rc in (124,143):
            rec['attempt_status']='TIMEOUT_NOT_EVIDENCE'
        else:
            rec['attempt_status']=f'SOLVER_ERROR_{rc}'
        if proof.exists() and proof.stat().st_size:
            rec['drat_raw_bytes']=proof.stat().st_size
            rec['drat_raw_sha256']=sha(proof)
            (d/'DRAT_RAW_BYTES.txt').write_text(str(proof.stat().st_size)+'\n')
            (d/'DRAT_RAW_SHA256.txt').write_text(rec['drat_raw_sha256']+'  '+proof.name+'\n')
            if rec['closure']=='VERIFIED_UNSAT_LEAF':
                zrc,_=run_cmd([self.a.zstd,'-T0','-8','--rm',str(proof)],min(600,max(1,int(self.remaining_wall()-2))),d/'zstd.log',d/'zstd_stderr.log')
                if zrc==0:
                    z=pathlib.Path(str(proof)+'.zst')
                    rec['drat_zst_bytes']=z.stat().st_size
                    rec['drat_zst_sha256']=sha(z)
                else:
                    rec['closure']='UNRESOLVED'
                    rec['attempt_status']=f'PROOF_COMPRESSION_ERROR_{zrc}'
            else:
                proof.unlink(missing_ok=True)
        (d/'ATTEMPT.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
        write_hashes(d)
        return rec
    def close(self,prefix,depth=0,parallel_here=False):
        self.ensure_wall(10)
        validate_prefix(prefix)
        nid=node_id(prefix)
        self.node_count+=1
        if self.node_count>self.a.max_nodes:
            rec={'node_id':nid,'prefix':prefix,'closure':'UNRESOLVED','attempt_status':'MAX_NODE_LIMIT'}
            self.nodes[nid]=rec;self.unresolved.append(prefix);return rec
        if len(prefix)==K:
            chk=direct_full_check(prefix)
            rec={'node_id':nid,'prefix':prefix,'prefix_length':K,'attempted_direct_cnf':False,'direct_check':chk}
            if chk['valid_b2_2']:
                rec['closure']='VERIFIED_SAT_COUNTEREXAMPLE'
                rec['attempt_status']='DIRECT_VALID_14SET'
                self.counterexamples.append({'prefix':prefix,'direct_check':chk})
            else:
                rec['closure']='DIRECTLY_REJECTED_FIXED_14SET'
                rec['attempt_status']='DIRECT_B2_VIOLATION'
            d=self.a.out/'nodes'/nid;d.mkdir(parents=True,exist_ok=True)
            (d/'DIRECT_CHECK.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
            write_hashes(d);self.nodes[nid]=rec;return rec
        rec=self.attempt_node(prefix,depth)
        if rec['closure'] in ('VERIFIED_UNSAT_LEAF','VERIFIED_SAT_COUNTEREXAMPLE'):
            if rec['closure']=='VERIFIED_SAT_COUNTEREXAMPLE':self.counterexamples.append(rec)
            self.nodes[nid]=rec;return rec
        vals=child_values(prefix)
        rec['split_dimension']='next_selected_residue'
        rec['expected_child_values']=vals
        rec['expected_child_count']=len(vals)
        results=[]
        def one(v):
            try:
                return self.close(prefix+[v],depth+1,False)
            except WallBudget:
                return {'node_id':node_id(prefix+[v]),'prefix':prefix+[v],'closure':'UNRESOLVED','attempt_status':'WALL_BUDGET_EXHAUSTED'}
        if parallel_here and self.a.workers>1 and len(vals)>1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.a.workers) as ex:
                futures={ex.submit(one,v):v for v in vals}
                for f in concurrent.futures.as_completed(futures):
                    results.append(f.result())
                    if any(r.get('closure')=='VERIFIED_SAT_COUNTEREXAMPLE' for r in results):
                        break
        else:
            for v in vals:
                results.append(one(v))
                if results[-1].get('closure')=='VERIFIED_SAT_COUNTEREXAMPLE':
                    break
        results.sort(key=lambda r:r['prefix'][-1])
        observed=[r['prefix'][-1] for r in results]
        missing=sorted(set(vals)-set(observed))
        duplicate=len(observed)!=len(set(observed))
        closed_kinds={'VERIFIED_UNSAT_LEAF','VERIFIED_UNSAT_BY_COMPLETE_CHILDREN','DIRECTLY_REJECTED_FIXED_14SET'}
        counter=[r for r in results if r.get('closure')=='VERIFIED_SAT_COUNTEREXAMPLE']
        unresolved=[r for r in results if r.get('closure') not in closed_kinds and r.get('closure')!='VERIFIED_SAT_COUNTEREXAMPLE']
        rec['children']=[r['node_id'] for r in results]
        rec['observed_child_values']=observed
        rec['missing_child_values']=missing
        rec['duplicate_child_values']=duplicate
        if counter:
            rec['closure']='VERIFIED_SAT_COUNTEREXAMPLE'
        elif not missing and not duplicate and len(results)==len(vals) and not unresolved:
            rec['closure']='VERIFIED_UNSAT_BY_COMPLETE_CHILDREN'
        else:
            rec['closure']='UNRESOLVED'
            self.unresolved.extend([r['prefix'] for r in unresolved])
            self.unresolved.extend([prefix+[v] for v in missing])
        self.nodes[nid]=rec
        d=self.a.out/'nodes'/nid;d.mkdir(parents=True,exist_ok=True)
        (d/'NODE_RESULT.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
        write_hashes(d)
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
    ap.add_argument('--root-solve-seconds',type=int,default=600)
    ap.add_argument('--min-solve-seconds',type=int,default=45)
    ap.add_argument('--root-check-seconds',type=int,default=240)
    ap.add_argument('--min-check-seconds',type=int,default=60)
    ap.add_argument('--decay',type=float,default=0.55)
    ap.add_argument('--wall-seconds',type=int,default=19000)
    ap.add_argument('--max-nodes',type=int,default=20000)
    a=ap.parse_args()
    prefix=json.loads(a.prefix_json)
    validate_prefix(prefix)
    a.out.mkdir(parents=True,exist_ok=True)
    c=Closer(a)
    try:
        root=c.close(prefix,0,True)
    except WallBudget:
        root={'node_id':node_id(prefix),'prefix':prefix,'closure':'UNRESOLVED','attempt_status':'WALL_BUDGET_EXHAUSTED_AT_ROOT'}
        c.unresolved.append(prefix)
    unresolved=[];seen=set()
    for p in c.unresolved:
        q=tuple(p)
        if q not in seen:
            unresolved.append(p);seen.add(q)
    payload={
      'schema':'b2-z100-recursive-proof-v4',
      'root_prefix':prefix,
      'root_closure':root.get('closure'),
      'root_node_id':root.get('node_id'),
      'node_count':c.node_count,
      'elapsed_seconds':time.monotonic()-c.start,
      'wall_budget_seconds':a.wall_seconds,
      'unresolved_prefixes':unresolved,
      'verified_sat_counterexamples':c.counterexamples,
      'root_verified_unsat':root.get('closure') in ('VERIFIED_UNSAT_LEAF','VERIFIED_UNSAT_BY_COMPLETE_CHILDREN','DIRECTLY_REJECTED_FIXED_14SET'),
      'claim_boundary':'The root is excluded only when root_verified_unsat is true and the complete proof tree, exact CNFs, checked proof logs, hashes, and direct fixed-set checks are preserved.'
    }
    (a.out/'PROOF_TREE.json').write_text(json.dumps(c.nodes,indent=2,sort_keys=True)+'\n')
    (a.out/'SUMMARY.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    (a.out/'UNRESOLVED_PREFIXES.json').write_text(json.dumps(unresolved,indent=2)+'\n')
    (a.out/'COUNTEREXAMPLES.json').write_text(json.dumps(c.counterexamples,indent=2)+'\n')
    (a.out/'STATUS.txt').write_text(('VERIFIED_UNSAT' if payload['root_verified_unsat'] else ('VERIFIED_SAT_COUNTEREXAMPLE' if c.counterexamples else 'PARTIAL_UNRESOLVED'))+'\n')
    write_hashes(a.out)
    print(json.dumps(payload,indent=2,sort_keys=True))
    if c.counterexamples: raise SystemExit(10)
    if not payload['root_verified_unsat']: raise SystemExit(2)

if __name__=='__main__':
    main()
