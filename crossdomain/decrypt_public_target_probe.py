#!/usr/bin/env python3
from __future__ import annotations
import json, os, re, time
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

BASE='https://de-crypt.org'
QUERY=BASE+'/decrypt-web/RecordsQuery'
LIST=BASE+'/decrypt-web/RecordsList'

def hidden(soup,name):
    x=soup.find('input',{'name':name});return x.get('value','') if x else ''

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/decrypt_targets'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 public-research-probe/1.0'})
    q=s.get(QUERY,timeout=60);q.raise_for_status();qs=BeautifulSoup(q.text,'html.parser')
    common={'csrf_name':hidden(qs,'csrf_name'),'csrf_value':hidden(qs,'csrf_value'),'t':'records','action':'search','rules':''}
    searches=[]
    specs=[
      ('nondecrypted_cipher',{'x_record_type':'1','x_status':'2'}),
      ('partially_cipher',{'x_record_type':'1','x_status':'3'}),
      ('nondecrypted_transcription',{'x_record_type':'1','x_status':'2','x_document_types[]':'7'}),
      ('partial_transcription',{'x_record_type':'1','x_status':'3','x_document_types[]':'7'}),
    ]
    all_ids=set()
    for label,extra in specs:
        data={**common,**extra};r=s.post(LIST,data=data,timeout=90,allow_redirects=True)
        soup=BeautifulSoup(r.text,'html.parser');ids=sorted(set(int(x) for x in re.findall(r'/decrypt-web/RecordsView/(\d+)',r.text)))
        # Capture pagination links; current page suffices as a feasibility probe.
        searches.append({'label':label,'status':r.status_code,'final_url':r.url,'bytes':len(r.content),'ids_first_page':ids,'id_count_first_page':len(ids),
                         'title':soup.title.get_text(' ',strip=True) if soup.title else ''})
        all_ids.update(ids)
        (out/f'{label}.html').write_text(r.text,encoding='utf-8')
        time.sleep(.5)
    records=[]
    keywords=re.compile(r'alchem|medic|physic|recipe|remed|pharmac|iatro|disease|health|doctor|surgeon',re.I)
    for rid in sorted(all_ids)[:120]:
        r=s.get(f'{BASE}/decrypt-web/RecordsView/{rid}',timeout=60);text=BeautifulSoup(r.text,'html.parser').get_text(' ',strip=True)
        rec={'id':rid,'status':r.status_code,'bytes':len(r.content),'authentication_required':bool(re.search(r'Authentication required',text,re.I)),
             'mentions_transcription':bool(re.search(r'Transcription',text,re.I)),'mentions_plaintext':bool(re.search(r'Plaintext',text,re.I)),
             'medical_alchemical_keyword':bool(keywords.search(text)),'keyword_contexts':[]}
        for m in list(keywords.finditer(text))[:5]:rec['keyword_contexts'].append(text[max(0,m.start()-180):m.start()+500])
        records.append(rec);time.sleep(.12)
    result={'status':'PUBLIC_SEARCH_POST_EXECUTED','searches':searches,'unique_ids_first_pages':len(all_ids),'records_audited':len(records),
            'records_without_authentication_gate':sum(not x['authentication_required'] for x in records),
            'records_with_transcription_mention':sum(x['mentions_transcription'] for x in records),
            'medical_alchemical_keyword_records':[x for x in records if x['medical_alchemical_keyword']],
            'claim_boundary':'Only public search/results and record pages were queried. No authentication-protected document was accessed or bypassed. A target survives only if a ciphertext/transcription is openly retrievable and an exact key/plaintext can later be independently verified.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
