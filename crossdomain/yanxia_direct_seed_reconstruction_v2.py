#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from pathlib import Path
import requests

BASE='https://www.shidianguji.com/{prefix}book/{book}/chapter/{chapter}'
# Prespecified after a location-only web search for the requested source title; these are discovery locations, not outcome/novelty labels.
KNOWN=[
 ('zh/','7589183229018505251','1lql44ullvwnn'),
 ('zh/','7589183229018505251','1lql47ofy836y'),
 ('zh/','7589183229018505251','1lql4nk8whk4q'),
 ('zh/','7589183229018505251','1lql4sn1t44cz'),
 ('zh/','7589183229018505251','1lql4titao8q2'),
 ('zh/','7589183229018505251','1lql4tsip4k6c'),
 ('zh/','7589183229018505251','1lql4vbwmnhub'),
 ('','7589226716518678568','1lqv89y58adqn'),
 ('zh/','7589226716518678568','1lqv89m6e81lv'),
 ('','7589226716518678568','1lqv8b3zuul1c'),
 ('','7589226716518678568','1lqv8dbrc3tis'),
 ('zh/','7589226716518678568','1lqv8ditrmkhg'),
 ('','7589226716518678568','1lqv8hvfwvpq8'),
 ('','7589226716518678568','1lqv8i138mgbm'),
 ('','7589226716518678568','1lqv8ive83y8y'),
 ('','7589226716518678568','1lqv8lc32unf3'),
 ('','7589226716518678568','1lqv8m67cod1c'),
 ('','7589226716518678568','1lqv8p9uollig'),
 ('','7589226716518678568','1lqv8rblqrehu'),
 ('','7589226716518678568','1lqv8tro7gjz6'),
]
TARGETS=('煙霞聖效方','烟霞圣效方','烟霞聖效方','煙霞圣效方')
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)
FORMULA_RE=re.compile(r'([一-龥]{2,12}(?:丸|散|湯|汤|丹|膏|飲|饮|煎|方|餅|饼|錠|锭))')
SOURCE_END=re.compile(r'(?:方|書|书|論|论|錄|录|集驗方|集验方|神方|要方|秘方)$')

def pars(raw):
    out=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"');obj=json.loads(inner)
            text=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            if text:out.append({'id':m.group(1),'order':int(m.group(3)),'text':text})
        except Exception:pass
    return sorted(out,key=lambda x:x['order'])

def target_in(text):return any(t in text for t in TARGETS)
def clean_source(x):return re.sub(r'[\s：:。；;、，,]+','',x)

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/yanxia_v2'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 historical-medical-text-reconstruction/2.0'})
    pages=[];errors=[]
    for prefix,book,ch in KNOWN:
        url=BASE.format(prefix=prefix,book=book,chapter=ch)
        try:
            r=s.get(url,timeout=60);r.raise_for_status();pp=pars(r.text)
            pages.append({'url':url,'book':book,'chapter':ch,'paragraphs':pp,'raw_has_target':target_in(r.text)})
        except Exception as e:errors.append({'url':url,'error':repr(e)})
        time.sleep(.08)
    # Learn source-heading vocabulary from short paragraphs across all located pages. This is only used to stop a block.
    source_titles=set(TARGETS)
    for pg in pages:
        for p in pg['paragraphs']:
            x=clean_source(p['text'])
            if 2<=len(x)<=18 and SOURCE_END.search(x):source_titles.add(x)
    source_titles=sorted(source_titles,key=lambda x:(-len(x),x))
    blocks=[]
    for pg in pages:
        pp=pg['paragraphs'];starts=[]
        for i,p in enumerate(pp):
            if target_in(clean_source(p['text'])):starts.append(i)
        for i in starts:
            body=[];start=pp[i]
            # include the target heading, then continue until a later short source heading.
            for j in range(i,min(len(pp),i+100)):
                p=pp[j];x=clean_source(p['text'])
                if j>i and x in source_titles and not target_in(x):break
                body.append(p)
            text='\n'.join(x['text'] for x in body)
            names=sorted(set(FORMULA_RE.findall(text)))
            blocks.append({'url':pg['url'],'book':pg['book'],'chapter':pg['chapter'],'start_paragraph_id':start['id'],'start_order':start['order'],'paragraph_count':len(body),'characters':len(text),'formula_name_candidates':names,'text':text})
    # de-duplicate identical starts
    blocks=list({(b['url'],b['start_paragraph_id']):b for b in blocks}.values())
    formulas=sorted(set(n for b in blocks for n in b['formula_name_candidates']))
    result={'status':'YANXIA_DIRECT_LOCATIONS_EXTRACTED' if blocks else 'NO_BLOCKS_PARSED',
            'located_pages_requested':len(KNOWN),'located_pages_downloaded':len(pages),'pages_with_raw_target_string':sum(p['raw_has_target'] for p in pages),
            'explicit_target_starts':len(blocks),'paragraphs_captured':sum(b['paragraph_count'] for b in blocks),'characters_captured':sum(b['characters'] for b in blocks),'formula_name_candidate_count':len(formulas),'formula_name_candidates':formulas[:300],'errors':errors,
            'claim_boundary':'Web search was used only to locate pages containing the requested source title. Extracted blocks/formula names are witnesses to source occurrence, not evidence that a formula is unique, unpublished, or that the lost book is completely reconstructed. Prior-art and parallel-text collation are mandatory.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'blocks.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'formulas.json').write_text(json.dumps(formulas,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
