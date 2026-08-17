#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from collections import Counter
from pathlib import Path
import requests

API='https://zh.wikisource.org/w/api.php'
S=requests.Session(); S.headers['User-Agent']='CrossDomainResearchProbe/1.0 (public research; contact via GitHub)'

def api(params):
    p={'format':'json','formatversion':2,**params}; r=S.get(API,params=p,timeout=90);r.raise_for_status();return r.json()

def allprefix(prefix):
    out=[]; cont=None
    while True:
        params={'action':'query','list':'allpages','apprefix':prefix,'apnamespace':0,'aplimit':'max'}
        if cont: params['apcontinue']=cont
        j=api(params);out += [x['title'] for x in j['query']['allpages']]
        cont=j.get('continue',{}).get('apcontinue')
        if not cont:break
    return out

def content(title):
    j=api({'action':'query','prop':'revisions','titles':title,'rvprop':'content','rvslots':'main'})
    pages=j.get('query',{}).get('pages',[])
    if not pages:return ''
    return pages[0].get('revisions',[{}])[0].get('slots',{}).get('main',{}).get('content','')

def strip_markup(t):
    t=re.sub(r'<ref[^>]*>.*?</ref>','',t,flags=re.S|re.I);t=re.sub(r'<[^>]+>','',t)
    t=re.sub(r'\{\{[^{}]*\}\}','',t);t=re.sub(r'\[\[(?:[^\]|]+\|)?([^\]]+)\]\]',r'\1',t)
    return re.sub(r'[\'\"《》〈〉「」『』\s]+','',t)

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/uibangyuchwi'));out.mkdir(parents=True,exist_ok=True)
    titles=allprefix('醫方類聚/')
    # Ensure base citation-list page if prefix API behaves differently.
    for x in ['醫方類聚/引用諸書']:
        if x not in titles:titles.append(x)
    corpus={};
    for i,t in enumerate(titles):
        try: corpus[t]=content(t)
        except Exception as e: corpus[t]=f'ERROR:{e!r}'
        if i%30==0: print('fetched',i,'/',len(titles),flush=True)
        time.sleep(.05)
    (out/'pages.json').write_text(json.dumps(corpus,ensure_ascii=False),encoding='utf-8')
    citepage=corpus.get('醫方類聚/引用諸書','')
    # Extract linked or list-item candidate source titles from the citation list.
    candidates=[]
    for line in citepage.splitlines():
        if line.startswith(('*','**','#')):
            x=re.sub(r'^\*+|^#+','',line).strip()
            x=re.sub(r'\[\[(?:[^\]|]+\|)?([^\]]+)\]\]',r'\1',x)
            x=re.sub(r'\{\{.*?\}\}','',x);x=re.sub(r'<.*?>','',x)
            x=re.split(r'[：:（(]',x)[0].strip(' \t;；,，。')
            if 1<=len(strip_markup(x))<=30:candidates.append(strip_markup(x))
    candidates=list(dict.fromkeys(x for x in candidates if x))
    body='\n'.join(v for k,v in corpus.items() if k!='醫方類聚/引用諸書' and not v.startswith('ERROR:'))
    plain=strip_markup(body)
    counts=[]
    for c in candidates:
        counts.append({'source':c,'occurrences_raw':body.count(c),'occurrences_normalized':plain.count(strip_markup(c))})
    counts=sorted(counts,key=lambda x:(-x['occurrences_normalized'],-x['occurrences_raw'],x['source']))
    # Known positive control from the published Gangibang reconstruction.
    gang_raw=body.count('簡奇方'); gang_norm=plain.count('簡奇方')
    # Capture context windows for exact source witness inspection.
    contexts=[]
    for m in re.finditer('簡奇方',body):
        contexts.append(body[max(0,m.start()-120):m.start()+600].replace('\n',' '))
        if len(contexts)>=400:break
    result={'status':'PUBLIC_CC_CORPUS_ACCESSED','page_count':len(titles),'citation_source_candidates':len(candidates),
            'gangibang_occurrences_raw':gang_raw,'gangibang_occurrences_normalized':gang_norm,
            'published_control_target':'287 Uibangyuchwi Gangibang excerpts reported in restoration study; raw occurrence count need not equal excerpt count because one excerpt may name title multiple/zero times after editorial markup.',
            'top_source_counts':counts[:100],
            'next_exact_step':'Segment source-labelled passages and compare blinded reconstruction against the published Gangibang 287-excerpt corpus; then rank other lost source titles only after excluding prior reconstructions.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'source_counts.json').write_text(json.dumps(counts,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'gangibang_contexts.txt').write_text('\n\n---\n\n'.join(contexts),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
