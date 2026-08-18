#!/usr/bin/env python3
"""Fail-closed finalizer for the staged proof-producing B2/Z100 workflows."""
from __future__ import annotations
import argparse,csv,hashlib,json,pathlib

WITNESS=[0,5,7,31,58,61,62,63,72,80,84,91,97]
ALL_T=set(range(2,46))

def load(p): return json.loads(pathlib.Path(p).read_text())
def canonical(xs): return sorted(tuple(int(y) for y in x) for x in xs)
def require(cond,label,failures):
    if not cond: failures.append(label)

def witness_audit():
    counts=[0]*100
    for a in WITNESS:
        for b in WITNESS:
            if a!=b: counts[(a-b)%100]+=1
    return {
      'set':WITNESS,'size':len(WITNESS),'distinct':len(set(WITNESS)),
      'ordered_nonzero_difference_total':sum(counts[1:]),
      'max_multiplicity':max(counts[1:]),
      'valid':len(WITNESS)==13 and len(set(WITNESS))==13 and max(counts[1:])<=2,
      'counts':counts
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage1',required=True)
    ap.add_argument('--u-stage',required=True)
    ap.add_argument('--v-stage',required=True)
    ap.add_argument('--recursive-stage',required=True)
    ap.add_argument('--output-dir',type=pathlib.Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    s1,u,v,r=map(load,[a.stage1,a.u_stage,a.v_stage,a.recursive_stage])
    f=[]
    require(s1.get('schema')=='b2-z100-third-branch-drat-v2-summary','stage1 schema',f)
    require(s1.get('expected_third_branches')==44,'stage1 expected count',f)
    require(s1.get('observed_rows')==44,'stage1 observed count',f)
    require(not s1.get('verified_sat_counterexample_thirds'),'stage1 counterexample',f)
    unresolved_t=[int(x) for x in s1.get('unresolved_thirds',[])]
    require(len(unresolved_t)==len(set(unresolved_t)) and set(unresolved_t)<=ALL_T,'stage1 unresolved t domain',f)
    direct_t=sorted(ALL_T-set(unresolved_t))
    require(s1.get('counts',{}).get('VERIFIED_UNSAT',0)==len(direct_t),'stage1 verified count/complement mismatch',f)

    require(u.get('schema')=='b2-z100-adaptive-u-split-v3-summary','u schema',f)
    require(sorted(map(int,u.get('stage1_unresolved_thirds',[])))==sorted(unresolved_t),'u input linkage',f)
    require(not u.get('verified_sat_counterexamples'),'u counterexample',f)
    u_closed=sorted(map(int,u.get('proved_thirds_by_u_split',[])))
    u_unresolved=canonical(u.get('unresolved_prefixes',[]))
    require(all(len(p)==4 and p[:2]==(0,1) and p[2] in unresolved_t for p in u_unresolved),'u unresolved descriptor shape',f)
    represented=set(u_closed)|{p[2] for p in u_unresolved}
    require(represented==set(unresolved_t),'u closed/unresolved parent coverage',f)
    require(not (set(u_closed)&{p[2] for p in u_unresolved}),'u closed/unresolved overlap',f)

    require(v.get('schema')=='b2-z100-adaptive-v-split-v3-summary','v schema',f)
    v_input=canonical(v.get('input_unresolved_prefixes',[]))
    require(v_input==u_unresolved,'v input linkage',f)
    require(not v.get('verified_sat_counterexamples'),'v counterexample',f)
    require(not v.get('missing_parent_prefixes'),'v missing parent',f)
    v_closed=set(canonical(v.get('closed_parent_prefixes',[])))
    v_unresolved=canonical(v.get('unresolved_prefixes',[]))
    require(all(len(p)==5 and p[:4] in set(v_input) for p in v_unresolved),'v unresolved descriptor shape',f)
    represented_v=v_closed|{p[:4] for p in v_unresolved}
    require(represented_v==set(v_input),'v closed/unresolved parent coverage',f)
    require(not (v_closed&{p[:4] for p in v_unresolved}),'v closed/unresolved overlap',f)

    require(r.get('schema')=='b2-z100-recursive-residual-v4-summary','recursive schema',f)
    r_input=canonical(r.get('input_unresolved_prefixes',[]))
    require(r_input==v_unresolved,'recursive input linkage',f)
    require(not r.get('verified_sat_counterexamples'),'recursive counterexample',f)
    require(not r.get('missing_root_prefixes'),'recursive missing root',f)
    require(canonical(r.get('closed_root_prefixes',[]))==r_input,'recursive closed roots mismatch',f)
    require(not r.get('unresolved_prefixes'),'recursive unresolved leaves',f)
    require(r.get('all_input_roots_closed') is True,'recursive gate not closed',f)

    wa=witness_audit()
    require(wa['valid'],'13-witness invalid',f)
    require(wa['ordered_nonzero_difference_total']==13*12,'witness difference total',f)

    all_t_closed=(len(f)==0)
    status=('COMPUTATIONAL_PROOF_CANDIDATE_FOR_M_13_PENDING_INDEPENDENT_EXTERNAL_REVIEW' if all_t_closed else 'PARTIAL_COMPUTATIONAL_RESOLUTION_PENDING_EXTERNAL_REVIEW')
    payload={
      'schema':'b2-z100-final-proof-chain-v5',
      'status':status,
      'witness_valid':wa['valid'],
      'lower_bound':'M>=13' if wa['valid'] else 'not established',
      'upper_bound':'M<=13' if all_t_closed else 'not established by this chain',
      'exact_value':'M=13' if all_t_closed and wa['valid'] else None,
      'normalized_third_domain':[2,45],
      'normalized_third_count':44,
      'stage1_direct_verified_thirds':direct_t,
      'stage1_unresolved_thirds':unresolved_t,
      'u_stage_fully_closed_thirds':u_closed,
      'u_stage_unresolved_fourth_prefix_count':len(u_unresolved),
      'v_stage_fully_closed_fourth_prefix_count':len(v_closed),
      'v_stage_unresolved_fifth_prefix_count':len(v_unresolved),
      'recursive_closed_fifth_prefix_count':len(r_input),
      'global_branch_coverage_and_terminal_evidence_complete':all_t_closed,
      'no_timeout_used_as_unsat':all_t_closed,
      'no_verified_14_set_counterexample':not (s1.get('verified_sat_counterexample_thirds') or u.get('verified_sat_counterexamples') or v.get('verified_sat_counterexamples') or r.get('verified_sat_counterexamples')),
      'source_run_ids':{
        'third_stage':s1.get('run_id'),'u_stage':u.get('run_id'),
        'v_stage':v.get('run_id'),'recursive_stage':r.get('run_id')
      },
      'source_artifacts':{
        'third_stage':{'id':s1.get('artifact_id'),'digest':s1.get('artifact_digest'),'url':s1.get('artifact_url')},
        'u_stage':{'id':u.get('artifact_id'),'digest':u.get('artifact_digest'),'url':u.get('artifact_url')},
        'v_stage':{'id':v.get('artifact_id'),'digest':v.get('artifact_digest'),'url':v.get('artifact_url')},
        'recursive_stage':{'id':r.get('artifact_id'),'digest':r.get('artifact_digest'),'url':r.get('artifact_url')}
      },
      'failures':f,
      'claim_boundary':'A complete chain supports a computational proof claim for M=13. It is not a Lean proof, peer review, or independent external reproduction; external review remains pending.'
    }
    (a.output_dir/'FINAL_PROOF_CHAIN.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    (a.output_dir/'WITNESS_AUDIT.json').write_text(json.dumps({k:v for k,v in wa.items() if k!='counts'},indent=2,sort_keys=True)+'\n')
    with (a.output_dir/'WITNESS_DIFFERENCE_COUNTS.csv').open('w',newline='') as fh:
        w=csv.writer(fh);w.writerow(['difference','ordered_multiplicity'])
        for d in range(1,100):w.writerow([d,wa['counts'][d]])
    (a.output_dir/'STATUS.txt').write_text(status+'\n')
    (a.output_dir/'FAILURES.txt').write_text('\n'.join(f)+('\n' if f else ''))
    print(json.dumps(payload,indent=2,sort_keys=True))
    return 0 if all_t_closed else 1

if __name__=='__main__': raise SystemExit(main())
