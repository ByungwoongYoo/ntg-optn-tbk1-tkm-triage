#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re,random,math,time
from pathlib import Path
from collections import defaultdict
import numpy as np
import requests
from bs4 import BeautifulSoup
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.model_selection import train_test_split

KREF=re.compile(r'\bK\.?\s*0*(\d{2,5})\b',re.I)

def med(f): return 'medicine' in str(f.get('genres') or '').lower()
def note(f):
    x=f.get('notes') or {}; return x.get('text','') if isinstance(x,dict) else str(x)
def dim(f,k):
    x=f.get(k)
    if isinstance(x,dict): return x.get('value')
    return x if isinstance(x,(int,float)) else None
def norm_lines(f):
    out=[]
    for s in str(f.get('atf') or '').splitlines():
        if re.match(r"\s*\d+['.]?",s):
            out.append(re.sub(r'[^0-9A-Za-zšṣṭḫĝŋ]+',' ',s.lower()).strip())
    return out

def line_jaccard(a,b):
    best=0.0
    for x in norm_lines(a):
        X=set(x.split())
        if not X: continue
        for y in norm_lines(b):
            Y=set(y.split())
            if Y: best=max(best,len(X&Y)/max(1,len(X|Y)))
    return best

def feat(a,b,sign_jac=None):
    sa=set(str(a.get('signs') or '').split()); sb=set(str(b.get('signs') or '').split())
    sj=len(sa&sb)/max(1,len(sa|sb)) if sign_jac is None else sign_jac
    lj=line_jaccard(a,b); ds=[]
    for k in ('width','length','thickness'):
        x,y=dim(a,k),dim(b,k); ds.append(abs(x-y) if x is not None and y is not None else np.nan)
    ca='@colophon' in str(a.get('atf') or '').lower(); cb='@colophon' in str(b.get('atf') or '').lower()
    return [sj,lj,*ds,int(ca),int(cb),int(ca and cb),int(med(a)),int(med(b))]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fragments-json',required=True);ap.add_argument('--out-dir',required=True);a=ap.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    data=json.loads(Path(a.fragments_json).read_text(encoding='utf-8')); fr=data if isinstance(data,list) else data.get('fragments',[])
    id2={f['_id']:f for f in fr}; ids=set(id2)
    direct=set(); uncertain=set(); explicit_neg=set(); globally_mentioned=set()
    for f in fr:
        src=f['_id']; nt=note(f); pub=str(f.get('publication') or '')
        for clause in re.split(r'[\n;]',nt+';'+pub):
            refs=['K.'+str(int(n)) for n in KREF.findall(clause)]; refs=[r for r in refs if r in ids and r!=src]
            for r in refs: globally_mentioned.add(tuple(sorted((src,r))))
            low=clause.lower()
            if re.search(r'\b(?:joins?|joined(?:\s+to)?|join\s+by)\b',clause,re.I):
                if any(x in low for x in ('no join','non-physical','indirect','might','maybe','possible','probably','sandwich','?')):
                    for r in refs: uncertain.add(tuple(sorted((src,r))))
                else:
                    for r in refs: direct.add(tuple(sorted((src,r))))
            if any(x in low for x in ('no join','non-physical','indirect join','parallel','dup.','duplicate')):
                for r in refs: explicit_neg.add(tuple(sorted((src,r))))
    pos=sorted(direct-explicit_neg-uncertain); neg0=sorted(explicit_neg-direct)
    # Build random matched negatives to keep classifier from learning only edited examples.
    available=[x for x,f in id2.items() if str(f.get('signs') or '').strip()]
    rng=random.Random(20260817); neg=set(neg0)
    target=max(len(neg0),3*len(pos))
    attempts=0
    while len(neg)<target and attempts<300000:
        attempts+=1; x,y=rng.sample(available,2); A,B=id2[x],id2[y]
        if A.get('museum')!=B.get('museum') or A.get('collection')!=B.get('collection'): continue
        if (A.get('script') or {}).get('period')!=(B.get('script') or {}).get('period'): continue
        p=tuple(sorted((x,y)))
        if p in direct or p in uncertain or p in globally_mentioned: continue
        neg.add(p)
    pairs=pos+sorted(neg); Y=np.array([1]*len(pos)+[0]*len(neg))
    X=np.array([feat(id2[x],id2[y]) for x,y in pairs],float)
    tr,te=train_test_split(np.arange(len(Y)),test_size=.30,random_state=20260817,stratify=Y)
    model=make_pipeline(SimpleImputer(strategy='median'),HistGradientBoostingClassifier(max_depth=4,learning_rate=.07,max_iter=220,l2_regularization=1,random_state=20260817))
    model.fit(X[tr],Y[tr]); pp=model.predict_proba(X[te])[:,1]
    held={'roc_auc':float(roc_auc_score(Y[te],pp)),'average_precision':float(average_precision_score(Y[te],pp)),'n_test':int(len(te))}
    # Retrieval pool from medical fragments using sign n-gram similarity, but final rank uses physical-join classifier.
    texts=[' '.join(str(f.get('signs') or '').split()) for f in fr]
    vec=TfidfVectorizer(token_pattern=r'(?u)\b\S+\b',ngram_range=(1,3),min_df=2,max_features=120000,sublinear_tf=True,norm='l2'); Z=vec.fit_transform(texts)
    ididx={f['_id']:i for i,f in enumerate(fr)}; mids=[i for i,f in enumerate(fr) if med(f) and texts[i]]
    cand={}
    for st in range(0,len(mids),30):
        ix=mids[st:st+30]; sims=(Z[ix]@Z.T).toarray()
        for rr,ii in enumerate(ix):
            F=fr[ii]; src=F['_id']; sims[rr,ii]=-1; k=min(250,len(fr)-1); top=np.argpartition(-sims[rr],k)[:k];top=top[np.argsort(-sims[rr,top])]
            saved=0
            for jj in top:
                G=fr[jj]; dst=G['_id']; p=tuple(sorted((src,dst)))
                if p in direct or p in uncertain or p in globally_mentioned: continue
                if F.get('museum')!=G.get('museum') or F.get('collection')!=G.get('collection'): continue
                if (F.get('script') or {}).get('period')!=(G.get('script') or {}).get('period'): continue
                s=float(sims[rr,jj]); ff=feat(F,G); prob=float(model.predict_proba(np.array(ff,float).reshape(1,-1))[0,1])
                rec={'a':p[0],'b':p[1],'physical_join_model_score':prob,'sign_tfidf_similarity':s,'features':{'sign_jaccard':ff[0],'max_line_jaccard':ff[1],'width_abs_diff':ff[2],'length_abs_diff':ff[3],'thickness_abs_diff':ff[4],'both_colophon':bool(ff[7])},
                     'a_medical':med(id2[p[0]]),'b_medical':med(id2[p[1]]),'a_cdli':(id2[p[0]].get('externalNumbers') or {}).get('cdliNumber'),'b_cdli':(id2[p[1]].get('externalNumbers') or {}).get('cdliNumber'),
                     'a_publication':id2[p[0]].get('publication'),'b_publication':id2[p[1]].get('publication'),'a_notes':note(id2[p[0]])[:1000],'b_notes':note(id2[p[1]])[:1000],
                     'a_atf':str(id2[p[0]].get('atf') or '')[:1800],'b_atf':str(id2[p[1]].get('atf') or '')[:1800]}
                old=cand.get(p)
                if old is None or prob>old['physical_join_model_score']:cand[p]=rec
                saved+=1
                if saved>=8: break
    ranked=sorted(cand.values(),key=lambda r:(-r['physical_join_model_score'],-r['sign_tfidf_similarity']))[:75]
    # Fetch CDLI public artifact pages for top 20 and record image/asset links without claiming visual contact.
    ses=requests.Session();ses.headers.update({'User-Agent':'Mozilla/5.0 eBL-medical-join-audit/1.0'})
    for r in ranked[:20]:
        r['cdli_pages']={}
        for side in ('a','b'):
            pnum=r.get(side+'_cdli') or ''
            num=re.sub(r'\D','',str(pnum))
            if not num: continue
            url=f'https://cdli.earth/artifacts/{num}'
            try:
                q=ses.get(url,timeout=45); soup=BeautifulSoup(q.text,'html.parser')
                links=[]
                for tag in soup.find_all(['a','img']):
                    u=tag.get('href') or tag.get('src') or ''
                    if re.search(r'asset|image|photo|\.jpe?g|\.png|\.tif',u,re.I):links.append(u)
                r['cdli_pages'][side]={'url':url,'status':q.status_code,'asset_like_links':list(dict.fromkeys(links))[:30]}
            except Exception as e:r['cdli_pages'][side]={'url':url,'error':repr(e)}
            time.sleep(.1)
    result={'status':'STRICT_PHYSICAL_JOIN_TRIAGE_COMPLETE','fragment_count':len(fr),'medical_fragments':sum(med(f) for f in fr),'direct_physical_training_pairs':len(pos),'explicit_negative_or_parallel_pairs':len(neg0),'training_total':len(Y),'heldout_random_split':held,
            'candidate_count':len(ranked),'top_candidates':[{k:v for k,v in x.items() if k not in ('a_atf','b_atf')} for x in ranked[:15]],
            'claim_boundary':'This is candidate triage, not a join discovery. Direct eBL note labels are used as noisy training data. Every candidate still requires targeted bibliography audit plus inspection of tablet photographs/edges and ideally curator/Assyriologist confirmation.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');(out/'candidates_v3.json').write_text(json.dumps(ranked,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
