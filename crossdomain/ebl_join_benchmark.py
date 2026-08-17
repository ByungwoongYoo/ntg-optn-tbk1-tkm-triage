#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer

KREF=re.compile(r'\bK\.?\s*0*(\d{2,5})\b',re.I)

def is_medical(f):
    s=str(f.get('genres') or '').lower()
    return any(x in s for x in ('medicine','medical','therapeutic'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fragments-json',required=True);ap.add_argument('--out-dir',required=True)
    a=ap.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    frags=json.loads(Path(a.fragments_json).read_text(encoding='utf-8'))
    if not isinstance(frags,list):frags=frags.get('fragments',[])
    ids={f['_id'] for f in frags}; id2idx={f['_id']:i for i,f in enumerate(frags)}
    strong=set(); anyrel=set()
    for f in frags:
        src=f['_id']; note=str(f.get('notes') or ''); pub=str(f.get('publication') or '')
        refs=[]
        for n in KREF.findall(note+' '+pub):
            rid='K.'+str(int(n))
            if rid in ids and rid!=src: refs.append(rid); anyrel.add(tuple(sorted((src,rid))))
        if re.search(r'\bjoins?\b|\bjoined\b',note,re.I):
            for rid in refs: strong.add(tuple(sorted((src,rid))))
    texts=[' '.join(str(f.get('signs') or '').split()) for f in frags]
    vec=TfidfVectorizer(token_pattern=r'(?u)\b\S+\b',ngram_range=(1,3),min_df=2,max_features=120000,sublinear_tf=True,norm='l2')
    X=vec.fit_transform(texts)
    adj=defaultdict(set)
    for x,y in strong: adj[x].add(y);adj[y].add(x)
    sources=[s for s in adj if texts[id2idx[s]] and any(texts[id2idx[m]] for m in adj[s])]
    metrics=[]
    for start in range(0,len(sources),50):
        ss=sources[start:start+50]; idx=[id2idx[s] for s in ss]; sims=(X[idx]@X.T).toarray()
        for r,s in enumerate(ss):
            sims[r,id2idx[s]]=-1
            mates=[id2idx[m] for m in adj[s] if m in id2idx]
            best=float(np.max(sims[r,mates])); rank=int(np.sum(sims[r]>best)+1)
            metrics.append({'source':s,'best_known_join_rank':rank,'best_known_join_similarity':best,'known_mates':sorted(adj[s]),'medical':is_medical(frags[id2idx[s]])})
    def summarize(rows):
        ranks=np.array([x['best_known_join_rank'] for x in rows],int)
        return {'n_sources':len(rows),'top1':float(np.mean(ranks<=1)) if len(rows) else None,'top10':float(np.mean(ranks<=10)) if len(rows) else None,'top100':float(np.mean(ranks<=100)) if len(rows) else None,'median_rank':float(np.median(ranks)) if len(rows) else None}
    # Discovery candidates: medical source, same museum+collection+period, no existing explicit cross-reference in notes/publication.
    med_idx=[i for i,f in enumerate(frags) if is_medical(f) and texts[i]]; cand={}
    for start in range(0,len(med_idx),40):
        idx=med_idx[start:start+40]; sims=(X[idx]@X.T).toarray()
        for rr,ii in enumerate(idx):
            fa=frags[ii]; src=fa['_id']; sims[rr,ii]=-1
            k=min(150,len(frags)-1); top=np.argpartition(-sims[rr],k)[:k]; top=top[np.argsort(-sims[rr,top])]
            kept=0
            for jj in top:
                fb=frags[jj]; dst=fb['_id']; pair=tuple(sorted((src,dst)))
                if pair in anyrel:continue
                if fa.get('museum')!=fb.get('museum') or fa.get('collection')!=fb.get('collection'):continue
                pa=(fa.get('script') or {}).get('period');pb=(fb.get('script') or {}).get('period')
                if pa and pb and pa!=pb:continue
                score=float(sims[rr,jj]);
                if score<=0:continue
                key=pair; prev=cand.get(key)
                rec={'a':pair[0],'b':pair[1],'sign_tfidf_similarity':score,'a_medical':is_medical(frags[id2idx[pair[0]]]),'b_medical':is_medical(frags[id2idx[pair[1]]]),
                     'a_genres':frags[id2idx[pair[0]]].get('genres'),'b_genres':frags[id2idx[pair[1]]].get('genres'),
                     'a_collection':frags[id2idx[pair[0]]].get('collection'),'b_collection':frags[id2idx[pair[1]]].get('collection'),
                     'a_script':frags[id2idx[pair[0]]].get('script'),'b_script':frags[id2idx[pair[1]]].get('script'),
                     'a_dimensions':[frags[id2idx[pair[0]]].get(x) for x in ('width','length','thickness')], 'b_dimensions':[frags[id2idx[pair[1]]].get(x) for x in ('width','length','thickness')],
                     'a_atf_preview':str(frags[id2idx[pair[0]]].get('atf') or '')[:1200], 'b_atf_preview':str(frags[id2idx[pair[1]]].get('atf') or '')[:1200],
                     'warning':'High textual/sign similarity can indicate a parallel copy rather than a physical join. This is a discovery candidate only.'}
                if prev is None or score>prev['sign_tfidf_similarity']:cand[key]=rec
                kept+=1
                if kept>=3:break
    candidates=sorted(cand.values(),key=lambda x:-x['sign_tfidf_similarity'])[:100]
    result={'fragment_count':len(frags),'strong_explicit_join_pairs':len(strong),'all_explicit_cross_reference_pairs':len(anyrel),'medical_fragment_count':sum(is_medical(f) for f in frags),
            'retrieval_all':summarize(metrics),'retrieval_medical_sources':summarize([x for x in metrics if x['medical']]),
            'method':'TF-IDF of sign-token 1-3 grams only; notes/publication and join fields excluded from retrieval features.',
            'candidate_count_saved':len(candidates),
            'claim_boundary':'Known-join recovery validates candidate generation only. New pairs require literature audit and physical/image/tablet-expert verification before any discovery claim.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'known_join_retrieval.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'medical_join_candidates.json').write_text(json.dumps(candidates,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
