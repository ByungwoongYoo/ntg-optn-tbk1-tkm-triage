#!/usr/bin/env python3
from __future__ import annotations
import html as htmllib, json, os, re, time
from collections import Counter, defaultdict
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)

def parse_paragraphs(raw):
    rows=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"');obj=json.loads(inner)
            text=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            rows.append({'paragraph_id':m.group(1),'order':int(m.group(3)),'text':text})
        except Exception: pass
    rows.sort(key=lambda x:x['order']);return rows

def clean_title(x):
    return re.sub(r'[\s、，,。；;：:]+','',x).strip()

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/uibang_source_inventory_v2'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 research-source-inventory/1.0'})
    frontier=[];seed_raw=None
    for book,chapter in SEEDS:
        r=s.get(BASE.format(book=book,chapter=chapter),timeout=60);r.raise_for_status()
        if seed_raw is None:seed_raw=r.text
        prefix=chapter[:4];ids=[x for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(prefix)]
        frontier += [(book,x) for x in ids]
    frontier=list(dict.fromkeys(frontier+SEEDS))
    # Citation list is explicitly paragraphs 14..163 of volume 3 after '引用諸書' and before the table of contents.
    first=parse_paragraphs(seed_raw)
    source_titles=[]
    for p in first:
        if 14<=p['order']<=163:
            for x in re.split(r'[、，,。；;]',p['text']):
                x=clean_title(x)
                if 1<len(x)<=16 and x not in source_titles:source_titles.append(x)
    # Longest-first prevents short source names from stealing a longer prefix.
    source_titles=sorted(source_titles,key=lambda x:(-len(x),x))
    counts=Counter(); examples=defaultdict(list); pages=0;pars_total=0;errors=[]
    for n,(book,ch) in enumerate(frontier,1):
        try:
            url=BASE.format(book=book,chapter=ch);r=s.get(url,timeout=60)
            if r.status_code!=200:errors.append({'url':url,'status':r.status_code});continue
            if '醫方類聚' not in r.text:continue
            pars=parse_paragraphs(r.text);pages+=1;pars_total+=len(pars)
            for i,p in enumerate(pars):
                text=p['text']
                if not text:continue
                matched=None
                for title in source_titles:
                    if text.startswith(title):matched=title;break
                if matched:
                    counts[matched]+=1
                    if len(examples[matched])<12:
                        tail=''.join(q['text'] for q in pars[i:i+4])[:1800]
                        examples[matched].append({'url':url,'chapter_title_match':re.search(r'<title[^>]*>(.*?)</title>',r.text,re.I|re.S).group(1) if re.search(r'<title[^>]*>(.*?)</title>',r.text,re.I|re.S) else '',
                                                  'paragraph_order':p['order'],'context':tail})
            if n%25==0:print('pages',n,'/',len(frontier),'parsed',pages,'paragraphs',pars_total,flush=True)
            time.sleep(.10)
        except Exception as e:errors.append({'book':book,'chapter':ch,'error':repr(e)})
    ranked=[]
    for title in source_titles:
        ranked.append({'source':title,'explicit_paragraph_starts':counts[title],'examples':examples[title]})
    ranked.sort(key=lambda x:(-x['explicit_paragraph_starts'],x['source']))
    gang=next((x for x in ranked if x['source']=='簡奇方'),None)
    result={'status':'SOURCE_PREFIX_INVENTORY_COMPLETE','pages_parsed':pages,'paragraphs_parsed':pars_total,'citation_source_titles':len(source_titles),
            'sources_with_explicit_starts':sum(x['explicit_paragraph_starts']>0 for x in ranked),'gangibang_control':gang,
            'top_sources':[{'source':x['source'],'explicit_paragraph_starts':x['explicit_paragraph_starts']} for x in ranked[:100]],'errors':errors[:100],
            'interpretation':'Counts measure explicit source-title paragraph starts, not total excerpts. Uibangyuchwi often abbreviates or omits repeated source labels, so this is a ranking/inventory layer only.',
            'next_step':'For high-count sources, build alias dictionaries and source blocks, then perform an external prior-reconstruction/survival audit before calling any passage lost or newly recovered.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'source_inventory.json').write_text(json.dumps(ranked,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('pages_parsed','paragraphs_parsed','citation_source_titles','sources_with_explicit_starts','gangibang_control')},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
