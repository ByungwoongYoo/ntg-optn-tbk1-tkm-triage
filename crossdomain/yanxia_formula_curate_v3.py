#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from pathlib import Path
from collections import Counter,defaultdict
import requests

BASE='https://www.shidianguji.com/{prefix}book/{book}/chapter/{chapter}'
KNOWN=[
 ('zh/','7589183229018505251','1lql44ullvwnn'),('zh/','7589183229018505251','1lql47ofy836y'),('zh/','7589183229018505251','1lql4nk8whk4q'),('zh/','7589183229018505251','1lql4sn1t44cz'),('zh/','7589183229018505251','1lql4titao8q2'),('zh/','7589183229018505251','1lql4tsip4k6c'),('zh/','7589183229018505251','1lql4vbwmnhub'),('', '7589226716518678568','1lqv89y58adqn'),('zh/','7589226716518678568','1lqv89m6e81lv'),('', '7589226716518678568','1lqv8b3zuul1c'),('', '7589226716518678568','1lqv8dbrc3tis'),('zh/','7589226716518678568','1lqv8ditrmkhg'),('', '7589226716518678568','1lqv8hvfwvpq8'),('', '7589226716518678568','1lqv8i138mgbm'),('', '7589226716518678568','1lqv8ive83y8y'),('', '7589226716518678568','1lqv8lc32unf3'),('', '7589226716518678568','1lqv8m67cod1c'),('', '7589226716518678568','1lqv8p9uollig'),('', '7589226716518678568','1lqv8rblqrehu'),('', '7589226716518678568','1lqv8tro7gjz6')]
TARGETS=('煙霞聖效方','烟霞圣效方','烟霞聖效方','煙霞圣效方')
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)
# High-precision prescription-name endings. 方 is excluded because it is too often a source-title/generic word.
NAME_END='(?:丸|散|湯|汤|丹|膏|飲|饮|餅|饼|錠|锭)'
HEADER_RE=re.compile(r'^\s*([一-龥]{2,10}'+NAME_END+r')(?=\s*(?:治|療|疗|主|出|，|。|：|:|$))')
INLINE_RE=re.compile(r'(?:^|[。；;])\s*([一-龥]{2,10}'+NAME_END+r')(?=\s*(?:治|療|疗|主|出|，|。|：|:|$))')
BAD_PREFIX=('每服','毎服','用','右','以','同','又','或','服','煎','溫','温','冷','熱','热','水','酒','米','生薑','姜汁','空心','臨卧','临卧','分作','煉蜜','炼蜜','熔','熬','和藥','和药','不論','勿服','至','入','取','治')

def pars(raw):
    out=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"');obj=json.loads(inner);txt=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            if txt: out.append({'id':m.group(1),'order':int(m.group(3)),'text':txt})
        except Exception: pass
    return sorted(out,key=lambda x:x['order'])

def target(text): return any(t in text for t in TARGETS)
def clean_name(n): return re.sub(r'\s+','',n)
def plausible(n):
    n=clean_name(n)
    if not (2<=len(n)<=10):return False
    if any(n.startswith(x) for x in BAD_PREFIX):return False
    if re.search(r'[一二三四五六七八九十百千卜]{2,}',n):return False
    return True

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/yanxia_formula_v3'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 historical-medical-text-curation/3.0'})
    pages=[]
    for pref,book,ch in KNOWN:
        url=BASE.format(prefix=pref,book=book,chapter=ch);r=s.get(url,timeout=60);r.raise_for_status();pages.append({'url':url,'paragraphs':pars(r.text)});time.sleep(.05)
    # Target spans: from explicit target heading to next paragraph that looks like another short source heading.
    all_rows=[]; spans=[]
    for pg in pages:
        pp=pg['paragraphs'];
        for p in pp: all_rows.append({'url':pg['url'],**p})
        starts=[i for i,p in enumerate(pp) if target(p['text'])]
        for i in starts:
            body=[]
            for j in range(i,min(len(pp),i+100)):
                x=pp[j]['text'].strip()
                if j>i and len(re.sub(r'\W','',x))<=16 and x.endswith('方') and not target(x): break
                body.append(pp[j])
            spans.append({'url':pg['url'],'start_id':pp[i]['id'],'rows':body})
    found=[]
    for sp in spans:
        rows=sp['rows']
        for ix,p in enumerate(rows):
            candidates=[]
            m=HEADER_RE.search(p['text'])
            if m:candidates.append((m.group(1),2 if re.search(r'(?:治|療|疗)',p['text'][m.end():m.end()+3]) else 1))
            for q in INLINE_RE.finditer(p['text']):candidates.append((q.group(1),2 if re.search(r'(?:治|療|疗)',p['text'][q.end():q.end()+3]) else 1))
            for n,conf in candidates:
                n=clean_name(n)
                if not plausible(n):continue
                # Capture up to four following paragraphs, stopping at next high-confidence formula header.
                body=[p['text']]
                for k in range(ix+1,min(len(rows),ix+5)):
                    if HEADER_RE.search(rows[k]['text']):break
                    body.append(rows[k]['text'])
                found.append({'name':n,'confidence':conf,'url':sp['url'],'paragraph_id':p['id'],'order':p['order'],'text':'\n'.join(body)})
    # Deduplicate witnesses.
    uniq={ (r['name'],r['url'],r['paragraph_id']):r for r in found}; found=list(uniq.values())
    names=sorted(set(r['name'] for r in found)); alltext='\n'.join(r['text'] for r in all_rows)
    occ={n:alltext.count(n) for n in names}; counts=Counter(r['name'] for r in found)
    by=defaultdict(list)
    for r in found:by[r['name']].append(r)
    summary=[]
    for n in names:
        ws=by[n];summary.append({'name':n,'source_witnesses':len(ws),'max_confidence':max(x['confidence'] for x in ws),'occurrences_in_20_located_pages':occ[n],
                                 'appears_only_in_extracted_target_witnesses_within_local_page_set':occ[n]==sum(x['text'].count(n) for x in ws),'example':ws[0]})
    summary.sort(key=lambda x:(-x['max_confidence'],-int(x['appears_only_in_extracted_target_witnesses_within_local_page_set']),x['occurrences_in_20_located_pages'],x['name']))
    result={'status':'YANXIA_FORMULA_HIGH_PRECISION_CURATION_COMPLETE','pages':len(pages),'target_spans':len(spans),'curated_formula_names':len(summary),'high_confidence_names':sum(x['max_confidence']>=2 for x in summary),
            'locally_target_only_names':sum(x['appears_only_in_extracted_target_witnesses_within_local_page_set'] for x in summary),'top_candidates':summary[:40],
            'claim_boundary':'High-precision formula segmentation removes many obvious OCR false positives, but local uniqueness across 20 located Uibangyuchwi pages is not historical novelty. Before any discovery claim, each top formula must be collated against the full Uibangyuchwi corpus, other surviving medical books, modern reconstructions, and page images.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');(out/'curated_formula_witnesses.json').write_text(json.dumps(found,ensure_ascii=False,indent=2),encoding='utf-8');(out/'formula_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
