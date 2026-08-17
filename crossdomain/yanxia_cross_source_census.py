#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re,time
from collections import defaultdict,Counter
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)
YANXIA_VARIANTS={'煙霞聖效方','烟霞聖效方'}

def pars(raw):
    out=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"'); obj=json.loads(inner)
            text=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            out.append({'id':m.group(1),'order':int(m.group(3)),'text':text})
        except Exception: pass
    return sorted(out,key=lambda x:x['order'])

def source_titles(first):
    titles=[]
    for p in first:
        if 14<=p['order']<=163:
            for x in re.split(r'[、，,。；;]',p['text']):
                x=re.sub(r'[\s：:]+','',x)
                if 1<len(x)<=16 and x not in titles: titles.append(x)
    for x in YANXIA_VARIANTS:
        if x not in titles: titles.append(x)
    return sorted(titles,key=lambda x:(-len(x),x))

def source_start(text,titles):
    for t in titles:
        if text.startswith(t): return t
    return None

def canon_source(s):
    return '煙霞聖效方' if s in YANXIA_VARIANTS else (s or 'UNKNOWN')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--inventory',required=True); ap.add_argument('--out-dir',required=True); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    inv=json.loads(Path(a.inventory).read_text(encoding='utf-8'))
    names=[x['name'] for x in inv if x.get('strong_detail')]
    sess=requests.Session(); sess.headers.update({'User-Agent':'Mozilla/5.0 historical-medical-text-census/1.1'})
    front=[]; raw0=None
    for book,ch in SEEDS:
        r=sess.get(BASE.format(book=book,chapter=ch),timeout=60); r.raise_for_status(); raw0=raw0 or r.text
        front += [(book,x) for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(ch[:4])]
    front=list(dict.fromkeys(front+SEEDS)); titles=source_titles(pars(raw0))
    hits=defaultdict(list); errs=[]
    for n,(book,ch) in enumerate(front,1):
        try:
            url=BASE.format(book=book,chapter=ch); r=sess.get(url,timeout=60)
            if r.status_code!=200 or '醫方類聚' not in r.text: continue
            current=None
            for p in pars(r.text):
                s=source_start(p['text'],titles)
                if s: current=s
                for name in names:
                    if name in p['text']:
                        hits[name].append({'source_raw':current,'source':canon_source(current),'url':url,'paragraph_id':p['id'],'order':p['order'],'text':p['text']})
            if n%30==0: print(n,len(front),'hits',sum(map(len,hits.values())),flush=True)
            time.sleep(.07)
        except Exception as e: errs.append({'book':book,'chapter':ch,'error':repr(e)})
    rows=[]
    for name in names:
        hh=hits.get(name,[]); sc=Counter(h.get('source') or 'UNKNOWN' for h in hh)
        other={k:v for k,v in sc.items() if k!='煙霞聖效方'}
        rows.append({'name':name,'total_occurrences':len(hh),'yanxia_occurrences':sc.get('煙霞聖效方',0),'other_source_occurrences':sum(other.values()),'other_sources':other,'only_yanxia_within_scanned_uibang':sc.get('煙霞聖效方',0)>0 and not other,'hits':hh})
    rows.sort(key=lambda x:(not x['only_yanxia_within_scanned_uibang'],-x['yanxia_occurrences'],x['name']))
    res={'target':'煙霞聖效方','source_variants_canonicalized':sorted(YANXIA_VARIANTS),'strong_names_input':len(names),'candidate_chapter_urls':len(front),'names_found_anywhere':sum(bool(x['total_occurrences']) for x in rows),'names_with_yanxia_occurrence':sum(bool(x['yanxia_occurrences']) for x in rows),'names_only_yanxia_within_scanned_uibang':sum(x['only_yanxia_within_scanned_uibang'] for x in rows),'names_with_other_source_occurrences':sum(bool(x['other_source_occurrences']) for x in rows),'only_yanxia_names':[x['name'] for x in rows if x['only_yanxia_within_scanned_uibang']],'cross_source_names':[{'name':x['name'],'yanxia_occurrences':x['yanxia_occurrences'],'other_sources':x['other_sources']} for x in rows if x['yanxia_occurrences'] and x['other_source_occurrences']],'errors':errs[:100],'claim_boundary':'Only-Yanxia means no occurrence of the same exact formula name was found under another source heading in the scanned Uibangyuchwi pages after canonicalizing 煙/烟 in the Yanxia source title. It is not a world-literature novelty claim and does not exclude orthographic variants, OCR errors, or occurrences in other books.'}
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'formula_source_census.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
