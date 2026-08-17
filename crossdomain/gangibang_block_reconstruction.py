#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)

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
                x=re.sub(r'[\s：:]+','',x)
                if 1<len(x)<=16 and x not in titles:titles.append(x)
    return sorted(titles,key=lambda x:(-len(x),x))

def match_source(text,titles):
    for t in titles:
        if text.startswith(t):return t
    return None

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/gangibang_blocks'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 research-lost-medical-text/1.0'})
    front=[];raw0=None
    for book,ch in SEEDS:
        r=s.get(BASE.format(book=book,chapter=ch),timeout=60);r.raise_for_status();raw0=raw0 or r.text
        front += [(book,x) for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(ch[:4])]
    front=list(dict.fromkeys(front+SEEDS)); titles=source_titles(pars(raw0))
    blocks=[];errs=[]
    for n,(book,ch) in enumerate(front,1):
        try:
            url=BASE.format(book=book,chapter=ch);r=s.get(url,timeout=60)
            if r.status_code!=200 or '醫方類聚' not in r.text:continue
            pp=pars(r.text);current=None;block=None
            for p in pp:
                src=match_source(p['text'],titles)
                if src:
                    if block and block['source']=='簡奇方':blocks.append(block)
                    current=src;block={'source':src,'url':url,'chapter_id':ch,'start_paragraph_id':p['id'],'start_order':p['order'],'paragraphs':[]}
                if current=='簡奇方' and block:
                    block['paragraphs'].append({'id':p['id'],'order':p['order'],'text':p['text']})
            if block and block['source']=='簡奇方':blocks.append(block)
            if n%30==0:print(n,len(front),'blocks',len(blocks),flush=True)
            time.sleep(.1)
        except Exception as e:errs.append({'book':book,'chapter':ch,'error':repr(e)})
    # De-duplicate by URL/start paragraph and quantify content. The published control contains 290 formulae across Uibangyuchwi + two other books, so no exact 290 equality is expected here.
    uniq={ (b['url'],b['start_paragraph_id']):b for b in blocks};blocks=list(uniq.values())
    for b in blocks:
        text='\n'.join(x['text'] for x in b['paragraphs']);b['text']=text;b['paragraph_count']=len(b['paragraphs']);b['characters']=len(text)
        b['treatment_marker_count']=len(re.findall(r'治',text));b['prescription_right_marker_count']=len(re.findall(r'右(?:件|藥|爲|各|用|將)?',text))
        b.pop('paragraphs',None)
    total_text='\n\n'.join(b['text'] for b in blocks)
    result={'status':'GANGIBANG_BLOCKS_EXTRACTED','explicit_source_blocks':len(blocks),'paragraphs_in_blocks':sum(b['paragraph_count'] for b in blocks),'characters_in_blocks':sum(b['characters'] for b in blocks),
            'treatment_marker_count':sum(b['treatment_marker_count'] for b in blocks),'prescription_right_marker_count':sum(b['prescription_right_marker_count'] for b in blocks),
            'published_control':'Ahn 2008 reports 290 medical formulae recovered from Uibangyuchwi, Hyangyakjipseongbang, and Changjinjip combined. This run extracts Uibangyuchwi explicit source blocks only, so it is a component/control, not a 290-formula reproduction.',
            'errors':errs[:100],'claim_boundary':'Blocks are source-attributed text witnesses. Formula segmentation, OCR correction, other cited books, and comparison with Ahn 2008 are required before claiming a complete reconstruction.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'gangibang_blocks.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'gangibang_text.txt').write_text(total_text,encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
