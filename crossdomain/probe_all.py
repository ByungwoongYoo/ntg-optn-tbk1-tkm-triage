#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, re, subprocess, sys, zipfile
from collections import Counter
from pathlib import Path

import requests


def get_json(url: str):
    r=requests.get(url,timeout=90); r.raise_for_status(); return r.json()

def get_text(url: str):
    r=requests.get(url,timeout=90); r.raise_for_status(); return r.text

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out-dir',required=True); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    result={'tracks':{}}

    # T1 eBL/Nineveh: download the open Zenodo JSON and inspect joins + medical subset.
    try:
        meta=get_json('https://zenodo.org/api/records/10018951')
        files=meta.get('files',[])
        fj=next(f for f in files if f.get('key')=='fragments.json')
        url=fj['links']['self']
        p=out/'ebl_fragments.json'
        with requests.get(url,stream=True,timeout=180) as r:
            r.raise_for_status()
            with p.open('wb') as h:
                for chunk in r.iter_content(1<<20):
                    if chunk:h.write(chunk)
        data=json.loads(p.read_text(encoding='utf-8'))
        fragments=data.get('fragments',data if isinstance(data,list) else [])
        join_records=[]; med=[]; nonempty=0
        for f in fragments:
            joins=f.get('joins') or []
            if joins:
                nonempty+=1; join_records.append({'id':f.get('_id'),'joins':joins[:5]})
            genres=f.get('genres') or []
            cats=[]
            for g in genres:
                cats += g.get('category',[]) if isinstance(g,dict) else []
            if any(str(x).lower()=='medicine' for x in cats): med.append(f)
        result['tracks']['T1_nineveh_ebl']={
            'status':'DATA_ACCESSED', 'fragment_count':len(fragments),
            'fragments_with_nonempty_joins':nonempty,'medical_fragment_count':len(med),
            'example_join_records':join_records[:20],
            'medical_with_joins':sum(bool(f.get('joins')) for f in med),
            'medical_publication_counts':Counter(str(f.get('publication','')) for f in med).most_common(20),
            'json_bytes':p.stat().st_size,
        }
        (out/'ebl_medical_join_examples.json').write_text(json.dumps([
            {'id':f.get('_id'),'publication':f.get('publication'),'description':f.get('description'),'joins':f.get('joins'),
             'atf':(f.get('atf') or '')[:1000]} for f in med if f.get('joins')
        ][:200],ensure_ascii=False,indent=2),encoding='utf-8')
    except Exception as e:
        result['tracks']['T1_nineveh_ebl']={'status':'FAILED','error':repr(e)}

    # T3 DECRYPT: public query surface + count known metadata; no authentication bypass.
    try:
        html=get_text('https://de-crypt.org/decrypt-web/RecordsQuery')
        home=get_text('https://de-crypt.org/decrypt-web/')
        result['tracks']['T3_decrypt']={
            'status':'PUBLIC_WEB_ACCESS', 'query_html_bytes':len(html), 'home_html_bytes':len(home),
            'mentions_transcription':bool(re.search('transcription',html,re.I)),
            'mentions_decrypted':bool(re.search('decrypt',html,re.I)),
            'forms':re.findall(r'<form[^>]*>',html,re.I)[:10],
            'candidate_endpoints':sorted(set(re.findall(r'(?:href|action)=[\"\']([^\"\']+)',html,re.I)))[:100],
        }
        (out/'decrypt_records_query.html').write_text(html,encoding='utf-8')
    except Exception as e:
        result['tracks']['T3_decrypt']={'status':'FAILED','error':repr(e)}

    # T5 RePAIR: inspect open record only. Do not download licensed 1.1 GB/53.5 GB data without explicit acceptance.
    try:
        meta=get_json('https://zenodo.org/api/records/13993089')
        files=[{'key':f.get('key'),'size':f.get('size'),'checksum':f.get('checksum')} for f in meta.get('files',[])]
        result['tracks']['T5_repair']={
            'status':'METADATA_ONLY_LEGAL_GATE','files':files,
            'license_requires_user_acceptance':True,
            'reason':'Dataset record contains custom license terms; full data not downloaded in this automated run.'
        }
    except Exception as e:
        result['tracks']['T5_repair']={'status':'FAILED','error':repr(e)}

    # T6 Vesuvius: inspect public repository/catalog and test remote metadata access only.
    try:
        api=get_json('https://api.github.com/repos/ScrollPrize/open-data/contents')
        names=[x.get('name') for x in api]
        readme=get_text('https://raw.githubusercontent.com/ScrollPrize/open-data/main/README.md')
        result['tracks']['T6_vesuvius']={
            'status':'OPEN_DATA_TOOLING_ACCESS','root_files':names,
            'readme_bytes':len(readme),'mentions_zarr':'zarr' in readme.lower(),'mentions_s3':'s3' in readme.lower(),
        }
        (out/'vesuvius_open_data_README.md').write_text(readme,encoding='utf-8')
    except Exception as e:
        result['tracks']['T6_vesuvius']={'status':'FAILED','error':repr(e)}

    # T7 CASP17: active target inventory; submission itself needs registration/account.
    try:
        html=get_text('https://predictioncenter.org/casp17/targetlist.cgi')
        target_ids=sorted(set(re.findall(r'\b(?:T|H|R|M|D|L|E|A)\d{2,4}(?:v\d+)?\b',html)))
        # Keep current relevant deadline snippets for 2026-08-17 onward.
        snippets=[]
        for tid in target_ids:
            m=re.search(r'(.{0,250}\b'+re.escape(tid)+r'\b.{0,500})',html,re.S)
            if m and re.search(r'2026-08-(?:1[7-9]|2\d|3[01])',m.group(1)):
                snippets.append(re.sub('<[^>]+>',' ',m.group(1)))
        seq_index=get_text('https://predictioncenter.org/download_area/CASP17/sequences/')
        result['tracks']['T7_casp17']={
            'status':'PUBLIC_TARGET_LIST_ACCESS_REGISTRATION_GATE_FOR_SUBMISSION',
            'target_id_count':len(target_ids),'late_august_snippets':snippets[:40],
            'public_sequence_index_bytes':len(seq_index),'submission_requires_registration':True,
        }
        (out/'casp17_targetlist.html').write_text(html,encoding='utf-8')
    except Exception as e:
        result['tracks']['T7_casp17']={'status':'FAILED','error':repr(e)}

    # T2 Uibangyuchwi: verify public web-indexed corpus endpoints; no scraping of blocked/private interfaces.
    urls=['https://mediclassics.kr/','https://mediclassics.kr/books/','https://mediclassics.kr/search']
    probes=[]
    for url in urls:
        try:
            r=requests.get(url,timeout=45,allow_redirects=True)
            probes.append({'url':url,'status_code':r.status_code,'final_url':r.url,'bytes':len(r.content),'title':(re.search(r'<title[^>]*>(.*?)</title>',r.text,re.I|re.S).group(1).strip() if re.search(r'<title[^>]*>(.*?)</title>',r.text,re.I|re.S) else '')})
        except Exception as e: probes.append({'url':url,'error':repr(e)})
    result['tracks']['T2_uibangyuchwi']={'status':'WEB_ACCESS_PROBE','probes':probes}

    (out/'PROBE_RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
