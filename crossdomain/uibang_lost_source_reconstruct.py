#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from pathlib import Path
import requests

SEEDS=[('7589183229018505251','1lql41s1l5hc6'),('7589226716518678568','1lqv8opt0f95v')]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
PAR_RE=re.compile(r'"paragraphId":"([^"]+)","paragraphType":\d+,"content":"((?:\\.|[^"])*)".*?"inChapterOrder":(\d+)',re.S)
FORMULA_NAME_RE=re.compile(r'^.{0,18}(?:方|丸|散|湯|丹|膏|飲|酒|煎|粉|餅|貼|元|圓)(?:$|[：:治])')


def pars(raw):
    out=[]
    for m in PAR_RE.finditer(raw):
        try:
            inner=json.loads('"'+m.group(2)+'"')
            obj=json.loads(inner)
            text=''.join(str(x.get('content','')) for x in obj.get('lines',[])).strip()
            out.append({'id':m.group(1),'order':int(m.group(3)),'text':text})
        except Exception:
            pass
    return sorted(out,key=lambda x:x['order'])


def source_titles(first):
    titles=[]
    for p in first:
        if 14<=p['order']<=163:
            for x in re.split(r'[、，,。；;]',p['text']):
                x=re.sub(r'[\s：:]+','',x)
                if 1<len(x)<=16 and x not in titles:
                    titles.append(x)
    return sorted(titles,key=lambda x:(-len(x),x))


def match_source(text,titles,target):
    if text.startswith(target):
        return target
    for t in titles:
        if text.startswith(t):
            return t
    return None


def marker_windows(paragraphs):
    windows=[]
    for i,p in enumerate(paragraphs):
        text=p['text'].strip()
        if not text:
            continue
        formulaish=bool(FORMULA_NAME_RE.search(text))
        treatment=('治' in text[:80])
        right=('右' in text[:30])
        if formulaish or treatment or right:
            lo=max(0,i-1); hi=min(len(paragraphs),i+3)
            windows.append({
                'anchor_id':p['id'],'anchor_order':p['order'],
                'formula_name_like':formulaish,'treatment_marker':treatment,'right_marker':right,
                'window':'\n'.join(x['text'] for x in paragraphs[lo:hi])
            })
    return windows


def main():
    target=os.environ.get('TARGET_SOURCE','煙霞聖效方').strip()
    out=Path(os.environ.get('OUT_DIR','artifact/uibang_lost_source'))
    out.mkdir(parents=True,exist_ok=True)
    sess=requests.Session(); sess.headers.update({'User-Agent':'Mozilla/5.0 research-lost-medical-text/1.1'})
    front=[]; raw0=None
    for book,ch in SEEDS:
        r=sess.get(BASE.format(book=book,chapter=ch),timeout=60); r.raise_for_status(); raw0=raw0 or r.text
        front += [(book,x) for x in sorted(set(CHAPTER_RE.findall(r.text))) if x.startswith(ch[:4])]
    front=list(dict.fromkeys(front+SEEDS))
    titles=source_titles(pars(raw0))
    if target not in titles:
        titles=[target]+titles
    blocks=[]; errs=[]
    for n,(book,ch) in enumerate(front,1):
        try:
            url=BASE.format(book=book,chapter=ch)
            r=sess.get(url,timeout=60)
            if r.status_code!=200 or '醫方類聚' not in r.text:
                continue
            pp=pars(r.text); current=None; block=None
            for p in pp:
                src=match_source(p['text'],titles,target)
                if src:
                    if block and block['source']==target:
                        blocks.append(block)
                    current=src
                    block={'source':src,'url':url,'chapter_id':ch,'start_paragraph_id':p['id'],'start_order':p['order'],'paragraphs':[]}
                if current==target and block:
                    block['paragraphs'].append({'id':p['id'],'order':p['order'],'text':p['text']})
            if block and block['source']==target:
                blocks.append(block)
            if n%30==0:
                print(n,len(front),'target',target,'blocks',len(blocks),flush=True)
            time.sleep(.08)
        except Exception as e:
            errs.append({'book':book,'chapter':ch,'error':repr(e)})

    uniq={(b['url'],b['start_paragraph_id']):b for b in blocks}
    blocks=list(uniq.values())
    all_windows=[]
    for b in blocks:
        paras=b['paragraphs']
        text='\n'.join(x['text'] for x in paras)
        wins=marker_windows(paras)
        all_windows += [{**w,'url':b['url'],'chapter_id':b['chapter_id']} for w in wins]
        b['text']=text
        b['paragraph_count']=len(paras)
        b['characters']=len(text)
        b['treatment_marker_count']=len(re.findall(r'治',text))
        b['prescription_right_marker_count']=len(re.findall(r'右(?:件|藥|爲|各|用|將)?',text))
        b['marker_window_count']=len(wins)
        b.pop('paragraphs',None)

    control=None
    if target=='簡奇方':
        control={
            'reference_run':{'blocks':27,'paragraphs':360,'characters':15470},
            'observed':{
                'blocks':len(blocks),
                'paragraphs':sum(b['paragraph_count'] for b in blocks),
                'characters':sum(b['characters'] for b in blocks)
            }
        }
        control['exact_match']=(control['reference_run']==control['observed'])

    result={
        'status':'SOURCE_BLOCKS_EXTRACTED',
        'target_source':target,
        'candidate_chapter_urls':len(front),
        'explicit_source_blocks':len(blocks),
        'paragraphs_in_blocks':sum(b['paragraph_count'] for b in blocks),
        'characters_in_blocks':sum(b['characters'] for b in blocks),
        'treatment_marker_count':sum(b['treatment_marker_count'] for b in blocks),
        'prescription_right_marker_count':sum(b['prescription_right_marker_count'] for b in blocks),
        'marker_windows':len(all_windows),
        'formula_name_like_windows':sum(1 for w in all_windows if w['formula_name_like']),
        'control':control,
        'errors':errs[:100],
        'claim_boundary':'These are source-attributed quotation witnesses in Uibangyuchwi. Marker windows are heuristic segmentation aids, not authenticated formula counts. A lost-book reconstruction or novelty claim requires textual collation, OCR correction, other surviving witnesses, and a dedicated prior-literature audit.'
    }
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'source_blocks.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'marker_windows.json').write_text(json.dumps(all_windows,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'source_text.txt').write_text('\n\n'.join(b['text'] for b in blocks),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
