#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from collections import Counter
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
TARGETS=('煙霞聖效方','烟霞圣效方')
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\\\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\\d+,"content":"((?:\\\\.|[^"])*)".*?"inChapterOrder":(\\d+)',re.S)
FORMULA_RE=re.compile(r'([一-龥]{2,10}(?:丸|散|湯|汤|丹|膏|飲|饮|煎|方|餅|饼|錠|锭))')

def pars(raw):
    out=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"');obj=json.loads(inner)
            text=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            out.append({'id':m.group(1),'order':int(m.group(3)),'text':text})
        except Exception:pass
    return sorted(out,key=lambda x:x['order'])

def source_titles(first):
    titles=[]
    for p in first:
        if 14<=p['order']<=163:
            for x in re.split(r'[、，,。；;]',p['text']):
                x=re.sub(r'[\\s：:]+','',x)
                if 1<len(x)<=18 and x not in titles:titles.append(x)
    return sorted(titles,key=lambda x:(-len(x),x))

def match_source(text,titles):
    for t in titles:
        if text.startswith(t):return t
    return None

def is_target(s):
    return any(t in (s or '') for t in TARGETS)

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/yanxia_blocks'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 historical-medical-text-research/1.0'})
    front=[];raw0=None
    for book,ch in SEEDS:
        r=s.get(BASE.format(book=book,chapter=ch),timeout=60);r.raise_for_status();raw0=raw0 or r.text
        front += [(book,x) for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(ch[:4])]
    front=list(dict.fromkeys(front+SEEDS));titles=source_titles(pars(raw0))
    blocks=[];all_pars=[];errs=[]
    for n,(book,ch) in enumerate(front,1):
        try:
            url=BASE.format(book=book,chapter=ch);r=s.get(url,timeout=60)
            if r.status_code!=200 or '醫方類聚' not in r.text:continue
            pp=pars(r.text);all_pars.extend([{'url':url,'chapter_id':ch,**p} for p in pp])
            current=None;block=None
            for p in pp:
                src=match_source(p['text'],titles)
                if src:
                    if block and is_target(block['source']):blocks.append(block)
                    current=src;block={'source':src,'url':url,'chapter_id':ch,'start_paragraph_id':p['id'],'start_order':p['order'],'paragraphs':[]}
                if current and is_target(current) and block:
                    block['paragraphs'].append({'id':p['id'],'order':p['order'],'text':p['text']})
            if block and is_target(block['source']):blocks.append(block)
            if n%30==0:print(n,len(front),'yanxia_blocks',len(blocks),flush=True)
            time.sleep(.08)
        except Exception as e:errs.append({'book':book,'chapter':ch,'error':repr(e)})
    uniq={(b['url'],b['start_paragraph_id']):b for b in blocks};blocks=list(uniq.values())
    target_names=[]
    for b in blocks:
        text='\\n'.join(x['text'] for x in b['paragraphs']);b['text']=text;b['paragraph_count']=len(b['paragraphs']);b['characters']=len(text)
        names=[]
        for p in b['paragraphs']:
            names += FORMULA_RE.findall(p['text'])
        b['formula_name_candidates']=sorted(set(names));target_names+=names;b.pop('paragraphs',None)
    target_names=sorted(set(target_names))
    outside=Counter()
    target_para_ids={b['start_paragraph_id'] for b in blocks}
    for p in all_pars:
        # Conservative: count name occurrence anywhere in scanned corpus; this is not a source attribution by itself.
        for name in target_names:
            if name and name in p['text']:outside[name]+=1
    name_rows=[{'name':n,'scanned_corpus_occurrences':outside[n]} for n in target_names]
    result={'status':'YANXIA_EXPLICIT_BLOCKS_EXTRACTED' if blocks else 'YANXIA_NOT_FOUND_IN_SCANNED_FRONTIER',
            'target_variants':TARGETS,'source_title_inventory_contains_target':any(is_target(t) for t in titles),
            'explicit_source_blocks':len(blocks),'paragraphs_in_blocks':sum(b['paragraph_count'] for b in blocks),'characters_in_blocks':sum(b['characters'] for b in blocks),
            'formula_name_candidate_count':len(target_names),'formula_name_candidates':target_names[:200],
            'errors':errs[:100],
            'claim_boundary':'These are explicit source-attributed witnesses and heuristic formula-name candidates only. Formula boundaries, text variants, parallel-source collation, page-image verification, and prior-art audit are required before any lost-text or novel-formula claim.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'yanxia_blocks.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'formula_name_candidates.json').write_text(json.dumps(name_rows,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
