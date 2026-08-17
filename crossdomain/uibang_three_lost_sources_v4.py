#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from collections import defaultdict
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)
FORMULA_RE=re.compile(r'([一-龥]{2,12}(?:丸|散|湯|汤|丹|膏|飲|饮|煎|方|餅|饼|錠|锭))')
TARGET_ALIASES={
 '澹軒方':['澹軒方','澹轩方'],
 '醫林方':['醫林方','医林方'],
 '煙霞聖效方':['煙霞聖效方','烟霞聖效方','煙霞圣效方','烟霞圣效方'],
}

def parse_paragraphs(raw):
    rows=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"'); obj=json.loads(inner)
            txt=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            if txt: rows.append({'id':m.group(1),'order':int(m.group(3)),'text':txt})
        except Exception: pass
    rows.sort(key=lambda x:x['order']); return rows

def clean(x): return re.sub(r'[\s、，,。；;：:]+','',x).strip()

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/uibang_three_lost_v4')); out.mkdir(parents=True,exist_ok=True)
    s=requests.Session(); s.headers.update({'User-Agent':'Mozilla/5.0 lost-medical-source-reconstruction/4.0'})
    seed_raw=None; frontier=[]
    for book,ch in SEEDS:
        r=s.get(BASE.format(book=book,chapter=ch),timeout=60); r.raise_for_status()
        if seed_raw is None: seed_raw=r.text
        prefix=ch[:4]
        ids=[x for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(prefix)]
        frontier += [(book,x) for x in ids]
    frontier=list(dict.fromkeys(frontier+SEEDS))
    first=parse_paragraphs(seed_raw)
    source_titles=[]
    for p in first:
        if 14<=p['order']<=163:
            for x in re.split(r'[、，,。；;]',p['text']):
                x=clean(x)
                if 1<len(x)<=16 and x not in source_titles: source_titles.append(x)
    # aliases so OCR/simplification does not miss target starts
    all_target_aliases=[a for vals in TARGET_ALIASES.values() for a in vals]
    starts=defaultdict(list); errors=[]; pages=0; paragraphs=0
    source_prefixes=sorted(source_titles,key=lambda x:(-len(x),x))
    for idx,(book,ch) in enumerate(frontier,1):
        url=BASE.format(book=book,chapter=ch)
        try:
            r=s.get(url,timeout=60)
            if r.status_code!=200 or '醫方類聚' not in r.text: continue
            pars=parse_paragraphs(r.text); pages+=1; paragraphs+=len(pars)
            for i,p in enumerate(pars):
                c=clean(p['text']); canonical=None
                for name,aliases in TARGET_ALIASES.items():
                    if any(c.startswith(a) for a in aliases): canonical=name; break
                if not canonical: continue
                body=[]
                for j in range(i,min(len(pars),i+140)):
                    q=pars[j]; qc=clean(q['text'])
                    if j>i:
                        # stop at a later explicit source-title prefix, but not target aliases inside body
                        other=None
                        for st in source_prefixes:
                            if qc.startswith(st): other=st; break
                        if other and not any(qc.startswith(a) for a in TARGET_ALIASES[canonical]): break
                    body.append(q)
                text='\n'.join(q['text'] for q in body)
                starts[canonical].append({'url':url,'book':book,'chapter':ch,'start_paragraph_id':p['id'],'start_order':p['order'],'paragraph_count':len(body),'characters':len(text),'formula_name_candidates':sorted(set(FORMULA_RE.findall(text))),'text':text})
            if idx%30==0: print('progress',idx,'/',len(frontier),flush=True)
            time.sleep(.05)
        except Exception as e: errors.append({'url':url,'error':repr(e)})
    # dedupe exact starts
    result_sources={}; all_blocks={}
    for name in TARGET_ALIASES:
        blocks=list({(b['url'],b['start_paragraph_id']):b for b in starts[name]}.values())
        formulas=sorted(set(x for b in blocks for x in b['formula_name_candidates']))
        result_sources[name]={'explicit_blocks':len(blocks),'paragraphs_captured':sum(b['paragraph_count'] for b in blocks),'characters_captured':sum(b['characters'] for b in blocks),'formula_name_candidates':len(formulas),'sample_formula_names':formulas[:100]}
        all_blocks[name]=blocks
        (out/(name+'_blocks.json')).write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf-8')
    result={'status':'THREE_LOST_SOURCE_EXPLICIT_CORPORA_EXTRACTED','pages_parsed':pages,'paragraphs_parsed':paragraphs,'sources':result_sources,'errors':errors[:100],
      'historical_scope':'These three titles are treated as lost-source targets because historical bibliography (e.g. Zhongguo Yiji Kao) records 煙霞聖效方, 澹軒方, and 醫林方 as not extant and notes their citation in Uibangyuchwi.',
      'claim_boundary':'This reconstructs explicit source-attributed blocks only. It is not a claim of complete book reconstruction, novelty, or corrected critical text. Repeated source labels may be omitted in Uibangyuchwi and OCR can be corrupt. Page-image verification and modern Chinese/Korean/Japanese prior-art collation are required.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
if __name__=='__main__': main()
