#!/usr/bin/env python3
"""Fail-closed finalizer for the independent 44-root B2/Z100 proof chain.

This script does not trust a green GitHub status by itself. It checks exact
stage linkage, canonical coverage, absence of verified counterexamples, and
the 13-element witness. Terminal DRAT validation remains carried by each
preserved job artifact and is indexed by the final workflow.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

WITNESS=[0,5,7,31,58,61,62,63,72,80,84,91,97]

def canon(xs):
    return sorted({tuple(map(int,x)) for x in xs})

def sha(path:Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def verify_witness():
    a=WITNESS
    assert len(a)==13 and len(set(a))==13 and all(0<=x<100 for x in a)
    r=[0]*100
    for x in a:
        for y in a:
            if x!=y:r[(x-y)%100]+=1
    maximum=max(r[1:]);bad=[i for i in range(1,100) if r[i]>2]
    assert maximum<=2 and not bad
    return {'set':a,'size':len(a),'ordered_nonzero_differences':sum(r[1:]),'max_nonzero_multiplicity':maximum,'valid':True,'difference_counts':r[1:]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',required=True);ap.add_argument('--continuation',required=True);ap.add_argument('--normalization-proof',required=True);ap.add_argument('--redundant-cuts-proof',required=True);ap.add_argument('--out',required=True);args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    base=json.load(open(args.base));cont=json.load(open(args.continuation))
    assert base['schema']=='b2-z100-recursive-third-v6-summary',base
    assert cont['schema']=='b2-z100-recursive-third-continuation-v6b-summary',cont
    assert int(cont['base_run_id'])==int(base['run_id']),(cont.get('base_run_id'),base.get('run_id'))
    assert int(base['expected_roots'])==44 and int(base['observed_rows'])==44,base
    closed=set(map(int,base.get('closed_thirds',[])));unresolved_t=set(map(int,base.get('unresolved_thirds',[])))
    assert closed|unresolved_t==set(range(2,46)),(closed,unresolved_t)
    assert not closed&unresolved_t,(closed,unresolved_t)
    assert not base.get('verified_sat_counterexamples'),base
    assert not cont.get('verified_sat_counterexamples'),cont
    assert canon(base.get('unresolved_prefixes',[]))==canon(cont.get('input_unresolved_prefixes',[])),{'base_unresolved':base.get('unresolved_prefixes'),'continuation_input':cont.get('input_unresolved_prefixes')}
    complete=bool(base.get('all_44_roots_verified_unsat')) or bool(cont.get('all_44_roots_proof_complete'))
    assert complete,{'base_complete':base.get('all_44_roots_verified_unsat'),'continuation_complete':cont.get('all_44_roots_proof_complete'),'continuation_unresolved':cont.get('unresolved_prefixes'),'missing':cont.get('missing_input_prefixes')}
    if not base.get('all_44_roots_verified_unsat'):
        assert cont.get('all_input_prefixes_closed'),cont
        assert not cont.get('missing_input_prefixes'),cont
        assert not cont.get('unresolved_prefixes'),cont
        assert canon(cont.get('closed_input_prefixes',[]))==canon(cont.get('input_unresolved_prefixes',[])),cont
    norm=Path(args.normalization_proof);cuts=Path(args.redundant_cuts_proof)
    assert norm.is_file() and norm.stat().st_size>0
    assert cuts.is_file() and cuts.stat().st_size>0
    witness=verify_witness()
    status='COMPUTATIONAL_PROOF_CANDIDATE_FOR_M_13_PENDING_INDEPENDENT_EXTERNAL_REVIEW'
    payload={
      'schema':'b2-z100-final-recursive-proof-v6',
      'status':status,
      'problem':'maximum size of a B2[2] subset of Z/100Z',
      'conclusion':'M=13 computationally, pending independent external review',
      'lower_bound':{'M_at_least':13,'witness':witness},
      'upper_bound':{
        'M_at_most':13,
        'canonical_third_roots':list(range(2,46)),
        'base_run_id':int(base['run_id']),
        'continuation_run_id':int(cont['run_id']),
        'base_closed_thirds':sorted(closed),
        'base_unresolved_thirds':sorted(unresolved_t),
        'continuation_closed_prefix_count':len(cont.get('closed_input_prefixes',[])),
        'all_44_roots_proof_complete':True,
        'normalization_proof_sha256':sha(norm),
        'redundant_cuts_proof_sha256':sha(cuts)
      },
      'no_timeout_used_as_unsat':True,
      'verified_sat_counterexamples':[],
      'claim_boundary':'This is a proof-carrying computational resolution candidate. It is not a Lean proof, peer review, or independent external reproduction. Every terminal UNSAT claim must remain linked to its checked DRAT artifact; timeouts are explicitly non-evidence.'
    }
    (out/'FINAL_PROOF_CHAIN.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    (out/'STATUS.txt').write_text(status+'\n')
    (out/'WITNESS_VERIFICATION.json').write_text(json.dumps(witness,indent=2,sort_keys=True)+'\n')
    print(json.dumps(payload,indent=2,sort_keys=True))

if __name__=='__main__':main()
