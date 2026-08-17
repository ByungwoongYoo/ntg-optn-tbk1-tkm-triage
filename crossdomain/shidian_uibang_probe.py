#!/usr/bin/env python3
from __future__ import annotations
import json, os, re
from pathlib import Path
import requests
from bs4 import BeautifulSoup

URLS=[
'https://www.shidianguji.com/zh/book/7589183229018505251/chapter/1lql41s1l5hc6',
'https://www.shidianguji.com/zh/book/7589183229018505251/chapter/1lql4nk8whk4q',
'https://www.shidianguji.com/zh/book/7589226716518678568/chapter/1lqv8opt0f95v',
]

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/shidian_probe'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 research-access-probe/1.0'})
    allrec=[]
    for i,url in enumerate(URLS):
        r=s.get(url,timeout=60); rec={'url':url,'status':r.status_code,'bytes':len(r.content),'final_url':r.url}
        html=r.text; (out/f'page{i}.html').write_text(html,encoding='utf-8')
        soup=BeautifulSoup(html,'html.parser');
        rec['title']=soup.title.get_text(' ',strip=True) if soup.title else ''
        links=[]
        for a in soup.find_all('a',href=True):
            href=a['href'];txt=a.get_text(' ',strip=True)
            if '/chapter/' in href or '/book/' in href:links.append({'href':href,'text':txt[:100]})
        rec['chapter_links']=links[:1000]
        # Next/Nuxt embedded payloads and strings containing chapter ids.
        scripts='\n'.join(x.get_text() for x in soup.find_all('script'))
        rec['script_bytes']=len(scripts)
        rec['chapter_strings']=sorted(set(re.findall(r'(?:/chapter/|chapterId[\\\"\':= ]+)([0-9a-z]{8,30})',html,re.I)))[:2000]
        rec['book_strings']=sorted(set(re.findall(r'/book/(\d{10,30})',html)))
        rec['gangibang_count']=html.count('簡奇方')
        # Visible text sample around Gangibang.
        text=soup.get_text('\n',strip=True); pos=text.find('簡奇方');rec['gangibang_context']=text[max(0,pos-300):pos+1800] if pos>=0 else ''
        allrec.append(rec)
    result={'status':'PUBLIC_PAGE_PROBE','pages':allrec,'claim_boundary':'No bulk crawling was performed in this probe. If pages expose a stable public chapter graph, the next step will use a low request rate and store only source-labelled excerpts/metadata needed for reconstruction.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
