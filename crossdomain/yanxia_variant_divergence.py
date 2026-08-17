#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re
from pathlib import Path

PUNC_RE=re.compile(r'[\s，,。:：；;、「」『』（）()〔〕【】\[\]<>《》]+')

def norm(s,name=''):
    s=str(s or '')
    if name: s=s.replace(name,'')
    return PUNC_RE.sub('',s)

def ngrams(s,n=2):
    return {s[i:i+n] for i in range(max(0,len(s)-n+1))}

def jac(a,b):
    A,B=ngrams(a),ngrams(b)
    if not A and not B: return 1.0
    if not A or not B: return 0.0
    return len(A&B)/len(A|B)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--census',required=True); ap.add_argument('--out-dir',required=True); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(Path(a.census).read_text(encoding='utf-8'))
    ranked=[]
    for r in rows:
        name=r['name']; yh=[h for h in r.get('hits',[]) if h.get('source')=='煙霞聖效方']; oh=[h for h in r.get('hits',[]) if h.get('source') not in (None,'UNKNOWN','煙霞聖效方')]
        if not yh or not oh: continue
        pairs=[]
        for y in yh:
            ny=norm(y.get('text',''),name)
            if len(ny)<4: continue
            for o in oh:
                no=norm(o.get('text',''),name)
                if len(no)<4: continue
                sim=jac(ny,no)
                pairs.append({'similarity':sim,'other_source':o.get('source'),'yanxia_text':y.get('text',''),'other_text':o.get('text',''),'yanxia_url':y.get('url'),'other_url':o.get('url'),'yanxia_order':y.get('order'),'other_order':o.get('order'),'yanxia_norm_len':len(ny),'other_norm_len':len(no)})
        if not pairs: continue
        best=max(pairs,key=lambda x:x['similarity']); worst=min(pairs,key=lambda x:x['similarity'])
        ranked.append({'name':name,'yanxia_occurrences':len(yh),'other_occurrences':len(oh),'other_sources':sorted({x.get('source') for x in oh if x.get('source')}),'max_similarity_to_any_other':best['similarity'],'min_similarity_to_other':worst['similarity'],'best_match':best,'most_divergent_pair':worst,'candidate_divergence_score':1.0-best['similarity']})
    ranked.sort(key=lambda x:(-x['candidate_divergence_score'],-max(x['best_match']['yanxia_norm_len'],x['best_match']['other_norm_len']),x['name']))
    res={'names_compared':len(ranked),'high_divergence_no_close_match_lt_0_25':sum(x['max_similarity_to_any_other']<0.25 for x in ranked),'moderate_divergence_no_close_match_lt_0_50':sum(x['max_similarity_to_any_other']<0.50 for x in ranked),'top_candidates':[{'name':x['name'],'candidate_divergence_score':x['candidate_divergence_score'],'max_similarity_to_any_other':x['max_similarity_to_any_other'],'other_sources':x['other_sources'],'yanxia_text':x['best_match']['yanxia_text'],'closest_other_text':x['best_match']['other_text'],'closest_other_source':x['best_match']['other_source']} for x in ranked[:15]],'claim_boundary':'This is a paragraph-level character-bigram triage of same-name occurrences, not a critical edition. Low similarity can reflect different formulae sharing a name, OCR/orthographic variation, or paragraph-boundary differences. Candidates require manual source collation and ingredient-level segmentation before any claim of a Yanxia-specific prescription variant.'}
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'ranked_variant_candidates.json').write_text(json.dumps(ranked,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
