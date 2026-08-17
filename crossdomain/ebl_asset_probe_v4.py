#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

KEYWORDS=('image','photo','picture','url','museum','cdli','thickness','height','width','dimension','join','collection','provenance','publication','notes','accession','number')

def short(v,limit=1200):
    s=json.dumps(v,ensure_ascii=False) if not isinstance(v,str) else v
    return s if len(s)<=limit else s[:limit]+'…'

def harvest(obj,prefix=''):
    out={}
    if isinstance(obj,dict):
        for k,v in obj.items():
            p=f'{prefix}.{k}' if prefix else k
            if any(q in k.lower() for q in KEYWORDS):
                out[p]=short(v)
            if isinstance(v,(dict,list)):
                out.update(harvest(v,p))
    elif isinstance(obj,list):
        for i,v in enumerate(obj[:20]):
            if isinstance(v,(dict,list)):
                out.update(harvest(v,f'{prefix}[{i}]'))
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--fragments-json',required=True); ap.add_argument('--survivors',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--top',type=int,default=12); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    frags=json.loads(Path(a.fragments_json).read_text(encoding='utf-8'))
    if not isinstance(frags,list): frags=frags.get('fragments',[])
    by={f['_id']:f for f in frags}
    surv=json.loads(Path(a.survivors).read_text(encoding='utf-8'))[:a.top]
    records=[]
    for c in surv:
        fa,fb=by[c['a']],by[c['b']]
        ha,hb=harvest(fa),harvest(fb)
        asset_keys_a=sorted(k for k in ha if any(q in k.lower() for q in ('image','photo','picture','url','cdli')))
        asset_keys_b=sorted(k for k in hb if any(q in k.lower() for q in ('image','photo','picture','url','cdli')))
        records.append({
            'a':c['a'],'b':c['b'],'similarity':c.get('sign_tfidf_similarity'),
            'a_top_level_keys':sorted(fa.keys()),'b_top_level_keys':sorted(fb.keys()),
            'a_relevant_metadata':ha,'b_relevant_metadata':hb,
            'a_asset_keys':asset_keys_a,'b_asset_keys':asset_keys_b,
            'both_have_discoverable_asset_metadata':bool(asset_keys_a and asset_keys_b),
            'metadata_priority_score_v3':c.get('metadata_priority_score_v3')
        })
    result={
        'n_candidates_probed':len(records),
        'both_have_discoverable_asset_metadata':sum(r['both_have_discoverable_asset_metadata'] for r in records),
        'candidates':[{'a':r['a'],'b':r['b'],'similarity':r['similarity'],'both_have_discoverable_asset_metadata':r['both_have_discoverable_asset_metadata'],'a_asset_keys':r['a_asset_keys'],'b_asset_keys':r['b_asset_keys']} for r in records],
        'claim_boundary':'This probes whether the open eBL metadata exposes routes to images/identifiers/dimensions for manual physical-join validation. It does not validate a join.'
    }
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'candidate_asset_metadata.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
