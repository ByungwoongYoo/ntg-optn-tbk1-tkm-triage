#!/usr/bin/env python3
from __future__ import annotations
import html as htmllib, json, os, re, time
from collections import Counter
from pathlib import Path
import requests
from bs4 import BeautifulSoup

SEEDS=[
('7589183229018505251','1lql41s1l5hc6'),
('7589226716518678568','1lqv8opt0f95v'),
]
BASE='https://www.shidianguji.com/zh/book/{book}/chapter/{chapter}'
CHAPTER_RE=re.compile(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',re.I)
VOL_RE=re.compile(r'卷([一二三四五六七八九十百〇零○兩两0-9]+)')

def clean_raw_context(s):
    s=htmllib.unescape(s)
    s=s.replace('\\u7c21\\u5947\\u65b9','簡奇方')
    s=re.sub(r'\\[nrt]',' ',s);s=re.sub(r'\\"','"',s);s=re.sub(r'<[^>]+>',' ',s)
    return re.sub(r'\s+',' ',s).strip()

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/shidian_bulk'));out.mkdir(parents=True,exist_ok=True)
    sess=requests.Session();sess.headers.update({'User-Agent':'Mozilla/5.0 research-text-mining/1.0'})
    frontier=[]
    for book,chapter in SEEDS:
        r=sess.get(BASE.format(book=book,chapter=chapter),timeout=60);r.raise_for_status()
        ids=sorted(set(CHAPTER_RE.findall(r.text)))
        # Most actual chapter ids share the seed's prefix; keep the dominant 4-char prefix.
        prefix=chapter[:4]
        ids=[x for x in ids if x.startswith(prefix)]
        for x in ids:frontier.append((book,x))
    # Include seeds explicitly and deduplicate.
    frontier=list(dict.fromkeys(frontier+SEEDS))
    records=[];gang_contexts=[];source_heading_contexts=[];errors=[]
    source_pattern=re.compile(r'([一-龥]{1,12}(?:方|書|論|經|集|法|要|錄|訣|編|本草|脈|疏|傳|記|說|圖|鈔|抄|類|目))(?=[：:。\s<\\])')
    for n,(book,ch) in enumerate(frontier,1):
        try:
            url=BASE.format(book=book,chapter=ch);r=sess.get(url,timeout=60)
            if r.status_code!=200:errors.append({'url':url,'status':r.status_code});continue
            soup=BeautifulSoup(r.text,'html.parser');title=soup.title.get_text(' ',strip=True) if soup.title else ''
            if '醫方類聚' not in title:continue
            raw=r.text
            gc=raw.count('簡奇方') + raw.count('\\u7c21\\u5947\\u65b9')
            rec={'book':book,'chapter_id':ch,'title':title,'url':url,'bytes':len(r.content),'gangibang_raw_occurrences':gc}
            records.append(rec)
            # Extract only bounded contexts, not the whole third-party corpus.
            for needle in ['簡奇方','\\u7c21\\u5947\\u65b9']:
                for m in re.finditer(re.escape(needle),raw):
                    gang_contexts.append({'title':title,'url':url,'context':clean_raw_context(raw[max(0,m.start()-500):m.start()+3500])})
            # Capture a bounded set of source-like labels from decoded page text/scripts for source inventory only.
            decoded=clean_raw_context(raw)
            labels=[]
            for m in source_pattern.finditer(decoded):
                lab=m.group(1)
                if lab not in labels:labels.append(lab)
                if len(labels)>=300:break
            rec['source_like_labels']=labels
            if n%20==0:print('fetched',n,'/',len(frontier),'valid',len(records),'gang contexts',len(gang_contexts),flush=True)
            time.sleep(.12)
        except Exception as e:errors.append({'book':book,'chapter':ch,'error':repr(e)})
    label_counts=Counter()
    for r in records:label_counts.update(r.get('source_like_labels',[]))
    result={'status':'PUBLIC_BOUNDED_EXCERPT_EXTRACTION','candidate_chapter_urls':len(frontier),'valid_uibangyuchwi_pages':len(records),
            'total_page_bytes_fetched':sum(r['bytes'] for r in records),'gangibang_raw_occurrence_sum':sum(r['gangibang_raw_occurrences'] for r in records),
            'gangibang_context_count':len(gang_contexts),'published_positive_control':'Gangibang restoration reports 287 excerpts from Uibangyuchwi; this raw name/context count is a feasibility check, not yet an excerpt-level reproduction.',
            'top_source_like_labels':label_counts.most_common(200),'errors':errors[:100],
            'claim_boundary':'Only bounded source-name/context evidence and metadata are saved; no full third-party corpus is redistributed. Novel lost-text claims require excerpt-level segmentation and prior-reconstruction audit.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'chapter_inventory.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'gangibang_contexts.json').write_text(json.dumps(gang_contexts,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
