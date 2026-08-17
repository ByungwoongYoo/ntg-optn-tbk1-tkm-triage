#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path


def pubnorm(x): return re.sub(r'\s+',' ',str(x or '')).strip()
def dimval(d):
    try:return float((d or {}).get('value'))
    except:return None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fragments-json',required=True);ap.add_argument('--candidates',required=True);ap.add_argument('--out-dir',required=True);a=ap.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    frags=json.loads(Path(a.fragments_json).read_text(encoding='utf-8'))
    if not isinstance(frags,list):frags=frags.get('fragments',[])
    by={f['_id']:f for f in frags}; cand=json.loads(Path(a.candidates).read_text(encoding='utf-8'))
    tiers=[]
    for c in cand:
        fa,fb=by[c['a']],by[c['b']];pa,pb=pubnorm(fa.get('publication')),pubnorm(fb.get('publication'))
        ta,tb=dimval(fa.get('thickness')),dimval(fb.get('thickness'))
        # Strongly penalize two separately published tablets: top example K.2357.A/K.2421 is BAM 554 vs 555, not a new join.
        if pa and pb and pa!=pb: tier='C_BOTH_SEPARATELY_PUBLISHED'
        elif not pa and not pb: tier='A_BOTH_UNPUBLISHED'
        else: tier='B_ONE_UNPUBLISHED'
        thickness_delta=abs(ta-tb) if ta is not None and tb is not None else None
        # Physical joins normally need compatible tablet thickness; a loose <=1 cm filter only ranks, never proves.
        geometry_ok=(thickness_delta is None or thickness_delta<=1.0)
        rec={**c,'a_publication':pa,'b_publication':pb,'a_notes':str(fa.get('notes') or '')[:1500],'b_notes':str(fb.get('notes') or '')[:1500],
             'tier':tier,'thickness_delta_cm':thickness_delta,'loose_thickness_compatible':geometry_ok,
             'discovery_score':float(c['sign_tfidf_similarity'])*(1.0 if tier=='A_BOTH_UNPUBLISHED' else .85 if tier=='B_ONE_UNPUBLISHED' else .15)*(1.0 if geometry_ok else .2)}
        tiers.append(rec)
    tiers.sort(key=lambda x:-x['discovery_score'])
    result={'n_input':len(cand),'tier_counts':{},'top_candidates':tiers[:50],
            'known_false_positive_control':'K.2357.A vs K.2421 is downgraded because CDLI identifies them as separate published composites BAM 6, 554 and BAM 6, 555; textual similarity alone is insufficient.',
            'claim_boundary':'This is triage only. A candidate is not a join until image/edge/curvature/physical or expert confirmation and a literature/CDLI/BM audit find no prior join.'}
    from collections import Counter
    result['tier_counts']=dict(Counter(x['tier'] for x in tiers))
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'ranked_candidates_v2.json').write_text(json.dumps(tiers,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'tier_counts':result['tier_counts'],'top5':[{k:x[k] for k in ('a','b','sign_tfidf_similarity','tier','discovery_score','a_publication','b_publication','thickness_delta_cm')} for x in tiers[:5]]},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
