#!/usr/bin/env python3
"""Verify the complete fresh frozen-source replay for B2[2]/Z100.

This checker distrusts the aggregate summary. It independently verifies file
hashes, exact normalized branch descriptors, every recorded child range,
TSV/JSON agreement, source hashes, exit codes, node totals, the 13-set witness,
and the arithmetic normalization/pruning lemmas.

A pass supports a computational-proof candidate, pending independent external
review. It is not an external reproduction or a formal proof assistant result.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,re,sys
from collections import Counter
from pathlib import Path

V10_SOURCE='3000f45f9a699445f736815a01bde93dd174353e69c0892c2e3fd9578f3988a1'
V8_SOURCE='93cf177dcfe3f92e24c5f5750af97a73e90cd17673a154188ee7333c466c1d13'
V10_NODES=5880602065
V8_NODES=7787871291
TOTAL_NODES=13668473356
WITNESS=[0,5,7,31,58,61,62,63,72,80,84,91,97]

def sha(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def req(c,msg,errors):
    if not c:errors.append(msg)

def one_hash_line(p:Path):
    lines=p.read_text().strip().splitlines();assert len(lines)==1,(p,len(lines))
    a,b=lines[0].split(None,1);return a,b.strip()

def bval(x):
    if isinstance(x,bool):return x
    if isinstance(x,int):return bool(x)
    return str(x).lower() in ('1','true','yes')

def mathematical_self_checks():
    # Ordered differences and circular unordered-distance capacities are equivalent.
    def ordered_ok(A,n):
        c=Counter((a-b)%n for a in A for b in A if a!=b)
        return max(c.values(),default=0)<=2
    def circular_ok(A,n):
        c=Counter()
        A=sorted(A)
        for i,a in enumerate(A):
            for b in A[i+1:]:
                d=(b-a)%n;d=min(d,n-d);c[d]+=1
        return all(c[d]<=2 for d in range(1,n//2)) and c[n//2]<=1
    for n in (8,10,12):
        for k in range(1,min(6,n+1)):
            for A in itertools.combinations(range(n),k):
                assert ordered_ok(A,n)==circular_ok(A,n),(n,A)
    # If no modulo-10 unit difference occurs, support is one parity or two opposite classes.
    units={1,3,7,9}
    for mask in range(1,1<<10):
        S=[r for r in range(10) if mask>>r&1]
        if all((a-b)%10 not in units for a in S for b in S):
            par={x%2 for x in S}
            assert len(par)==1 or (len(S)==2 and (S[0]-S[1])%10==5),S
    # Parity capacity forces 5..9 even elements.
    feasible=[e for e in range(15) if e*(e-1)+(14-e)*(13-e)<=98]
    assert feasible==[5,6,7,8,9],feasible
    # Reflection range and descriptor counts.
    assert [t for t in range(2,100) if t+11<=101-t]==list(range(2,46))
    assert sum(len(range(t+1,92-t)) for t in (2,3,4))==255
    return {'ordered_circular_equivalence_small_instances':True,'unit_difference_mod10_support_check':True,'parity_feasible_even_counts':feasible,'normalized_thirds':[2,45],'v10_descriptor_count':255,'v8_descriptor_count':41}

def witness_check():
    counts=[0]*100
    for a in WITNESS:
        for b in WITNESS:
            if a!=b:counts[(a-b)%100]+=1
    return {'set':WITNESS,'size':len(WITNESS),'distinct':len(set(WITNESS)),'ordered_pairs':sum(counts[1:]),'max_multiplicity':max(counts[1:]),'valid':len(WITNESS)==13 and len(set(WITNESS))==13 and sum(counts[1:])==156 and max(counts[1:])<=2,'counts':counts}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('--v10-source',type=Path,required=True);ap.add_argument('--v8-source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--counts-csv',type=Path,required=True)
    a=ap.parse_args();root=a.root;errors=[]
    req(root.is_dir(),'replay root missing',errors)
    req(sha(a.v10_source)==V10_SOURCE,'v10 source SHA mismatch',errors)
    req(sha(a.v8_source)==V8_SOURCE,'v8 source SHA mismatch',errors)
    math=mathematical_self_checks();wit=witness_check();req(wit['valid'],'13-set witness invalid',errors)

    # Top aggregate manifest. Its self-line was generated while the file was open,
    # so only all non-self entries are treated as integrity claims.
    top_manifest=root/'ALL_FILES_SHA256SUMS.txt';top_checked=0
    req(top_manifest.is_file(),'aggregate manifest missing',errors)
    if top_manifest.is_file():
        for line in top_manifest.read_text().splitlines():
            if not line.strip():continue
            dig,rel=line.split(None,1);rel=rel.strip()
            if rel=='aggregate/ALL_FILES_SHA256SUMS.txt':continue
            req(rel.startswith('aggregate/'),'bad aggregate manifest path '+rel,errors)
            p=root/rel.removeprefix('aggregate/')
            req(p.is_file(),'aggregate manifest missing '+rel,errors)
            if p.is_file():req(sha(p)==dig,'aggregate manifest hash mismatch '+rel,errors)
            top_checked+=1

    exp10={(t,u) for t in (2,3,4) for u in range(t+1,92-t)};exp8=set(range(5,46))
    got10=set();got8=set();n10=n8=0;results=rows_total=inner_checked=0
    binaries={'v10':set(),'v8':set()};gxx={'v10':set(),'v8':set()};uname={'v10':set(),'v8':set()}
    for top in sorted(p for p in root.iterdir() if p.is_dir()):
        m10=re.fullmatch(r'fresh-v10-t(\d+)-u(\d+)',top.name);m8=re.fullmatch(r'fresh-v8-t(\d+)',top.name)
        req(bool(m10 or m8),'unexpected top directory '+top.name,errors)
        if not (m10 or m8):continue
        eng='v10' if m10 else 'v8';t=int((m10 or m8).group(1));u=int(m10.group(2)) if m10 else None
        child=top/(f't{t}_u{u}' if m10 else f't{t}')
        for x in ('BINARY_SHA256.txt','SOURCE_SHA256.txt','GXX_VERSION.txt','UNAME.txt'):req((top/x).is_file(),f'{top.name} missing {x}',errors)
        for x in ('RESULT.json','BRANCHES.tsv','EXIT_CODE.txt','STDOUT.log','TIME.log','SHA256SUMS.txt'):req((child/x).is_file(),f'{top.name} missing {x}',errors)
        if not (child/'RESULT.json').is_file():continue
        src,_=one_hash_line(top/'SOURCE_SHA256.txt');req(src==(V10_SOURCE if eng=='v10' else V8_SOURCE),f'{top.name} source SHA',errors)
        bh,_=one_hash_line(top/'BINARY_SHA256.txt');binaries[eng].add(bh);gxx[eng].add(sha(top/'GXX_VERSION.txt'));uname[eng].add(sha(top/'UNAME.txt'))
        req((child/'EXIT_CODE.txt').read_text().strip()=='0',f'{top.name} nonzero exit',errors)
        req('Exit status: 0' in (child/'TIME.log').read_text(errors='replace'),f'{top.name} time-log exit',errors)
        for line in (child/'SHA256SUMS.txt').read_text().splitlines():
            if not line.strip():continue
            dig,rel=line.split(None,1);rel=rel.strip();req(rel.startswith('artifact/'),f'{top.name} inner path',errors)
            p=top.joinpath(*Path(rel.removeprefix('artifact/')).parts);req(p.is_file(),f'{top.name} inner missing {rel}',errors)
            if p.is_file():req(sha(p)==dig,f'{top.name} inner hash {rel}',errors)
            inner_checked+=1
        d=json.loads((child/'RESULT.json').read_text());results+=1
        req(d.get('third')==t,f'{top.name} third',errors);req(d.get('completed_exhaustively') is True,f'{top.name} incomplete',errors);req(d.get('timed_out') is False,f'{top.name} timeout',errors);req(d.get('witness_found') is False,f'{top.name} witness',errors);req(d.get('witness')==[],f'{top.name} witness payload',errors)
        with (child/'BRANCHES.tsv').open(newline='') as f:trs=list(csv.DictReader(f,delimiter='\t'))
        bs=d.get('branches',[]);rows_total+=len(trs);req(len(trs)==len(bs),f'{top.name} TSV/JSON count',errors)
        if eng=='v10':
            got10.add((t,u));vals=list(range(u+1,93-t));req(t in (2,3,4) and u in range(t+1,92-t),f'{top.name} descriptor',errors);req(d.get('fourth')==u,f'{top.name} fourth',errors);req(d.get('fifth_range')==[u+1,92-t],f'{top.name} fifth range',errors);req(d.get('expected_fifth_branches')==len(vals)==d.get('observed_fifth_branches'),f'{top.name} fifth count',errors);req([b.get('fifth') for b in bs]==vals,f'{top.name} fifth coverage',errors);req([int(r['fifth']) for r in trs]==vals,f'{top.name} TSV fifth coverage',errors);n=int(d.get('aggregate_nodes_including_fourth_root',-1));n10+=n;req(n==1+sum(int(b['nodes']) for b in bs),f'{top.name} node sum',errors)
        else:
            got8.add(t);vals=list(range(t+1,92-t));req(t in exp8,f'{top.name} descriptor',errors);req(d.get('fourth_range')==[t+1,91-t],f'{top.name} fourth range',errors);req(d.get('expected_fourth_branches')==len(vals)==d.get('observed_fourth_branches'),f'{top.name} fourth count',errors);req([b.get('fourth') for b in bs]==vals,f'{top.name} fourth coverage',errors);req([int(r['fourth']) for r in trs]==vals,f'{top.name} TSV fourth coverage',errors);n=int(d.get('aggregate_nodes_including_third_root',-1));n8+=n;req(n==1+sum(int(b['nodes']) for b in bs),f'{top.name} node sum',errors)
        for b,r in zip(bs,trs):
            req(b.get('completed') is True and b.get('timed') is False and b.get('found') is False and b.get('verified') is False,f'{top.name} child status',errors)
            for jk,tk in (('nodes','nodes'),('prunes','prunes'),('qprunes','qprunes'),('lbprunes','lbprunes')):req(int(r[tk])==int(b[jk]),f'{top.name} TSV {tk}',errors)
            for jk,tk in (('completed','completed'),('timed','timed'),('found','found'),('verified','verified')):req(bval(int(r[tk]))==bool(b[jk]),f'{top.name} TSV {tk}',errors)
    req(got10==exp10,'v10 descriptor coverage mismatch',errors);req(got8==exp8,'v8 descriptor coverage mismatch',errors);req(results==296,f'result count {results}',errors);req(n10==V10_NODES,f'v10 nodes {n10}',errors);req(n8==V8_NODES,f'v8 nodes {n8}',errors);req(n10+n8==TOTAL_NODES,f'total nodes {n10+n8}',errors);req(len(binaries['v10'])==1 and len(binaries['v8'])==1,'binary reproducibility',errors);req(len(gxx['v10'])==len(gxx['v8'])==1,'compiler environment consistency',errors);req(len(uname['v10'])==len(uname['v8'])==1,'runner environment consistency',errors)
    summary=json.loads((root/'FRESH_REPLAY_SUMMARY.json').read_text());req(summary.get('status')=='FRESH_FULL_EXACT_REPLAY_PASS','aggregate summary status',errors);req(summary.get('total_nodes')==TOTAL_NODES,'aggregate summary nodes',errors);req(summary.get('timeouts')==0 and summary.get('witnesses')==0,'aggregate summary outcome',errors)

    a.counts_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.counts_csv.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['difference','ordered_multiplicity']);[w.writerow([i,wit['counts'][i]]) for i in range(1,100)]
    status='COMPUTATIONAL_PROOF_CANDIDATE_FOR_M_13_PENDING_INDEPENDENT_EXTERNAL_REVIEW' if not errors else 'AUDIT_FAILED'
    out={'schema':'b2-z100-fresh-replay-proof-audit-v6','status':status,'witness_valid':wit['valid'],'lower_bound':'M>=13' if wit['valid'] else None,'upper_bound':'M<=13' if not errors else None,'exact_value':'M=13' if not errors and wit['valid'] else None,'result_files':results,'v10_descriptors':len(got10),'v8_descriptors':len(got8),'branch_rows':rows_total,'v10_nodes':n10,'v8_nodes':n8,'total_nodes':n10+n8,'timeouts':0 if not errors else None,'verified_14_set_counterexamples':0 if not errors else None,'aggregate_manifest_entries_checked':top_checked,'inner_manifest_entries_checked':inner_checked,'source_sha256':{'v10':V10_SOURCE,'v8':V8_SOURCE},'binary_hash_counts':{k:len(v) for k,v in binaries.items()},'mathematical_self_checks':math,'witness':{k:v for k,v in wit.items() if k!='counts'},'errors':errors,'claim_boundary':'This establishes an internally audited frozen-source computational proof candidate for M=13. The complete search was freshly replayed by the same research project; independent external reproduction, peer review, and formal proof-assistant verification remain pending.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0 if not errors else 1
if __name__=='__main__':raise SystemExit(main())
