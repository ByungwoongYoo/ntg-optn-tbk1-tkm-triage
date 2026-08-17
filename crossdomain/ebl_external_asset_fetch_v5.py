#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re,time
from pathlib import Path
import requests

IMG_RE=re.compile(r'https?://[^\"\'<> ]+?(?:\.jpg|\.jpeg|\.png|\.webp)(?:\?[^\"\'<> ]*)?',re.I)
META_IMG_RE=re.compile(r'<meta[^>]+(?:property|name)=[\"\'](?:og:image|twitter:image)[\"\'][^>]+content=[\"\']([^\"\']+)',re.I)

def fetch(s,url):
    try:
        r=s.get(url,timeout=40,allow_redirects=True)
        txt=r.text[:2_000_000] if 'text' in r.headers.get('content-type','') or 'html' in r.headers.get('content-type','') else ''
        imgs=list(dict.fromkeys(META_IMG_RE.findall(txt)+IMG_RE.findall(txt)))[:20]
        return {'url':url,'status':r.status_code,'final_url':r.url,'content_type':r.headers.get('content-type',''),'bytes':len(r.content),'title':(re.search(r'<title[^>]*>(.*?)</title>',txt,re.I|re.S).group(1).strip() if re.search(r'<title[^>]*>(.*?)</title>',txt,re.I|re.S) else ''),'image_urls':imgs}
    except Exception as e:
        return {'url':url,'error':repr(e)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--asset-metadata',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--top',type=int,default=3); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(Path(a.asset_metadata).read_text(encoding='utf-8'))[:a.top]
    s=requests.Session(); s.headers.update({'User-Agent':'Mozilla/5.0 open-research-join-audit/1.0'})
    results=[]
    for row in rows:
        rr={'a':row['a'],'b':row['b'],'sides':[]}
        for side in ('a','b'):
            md=row[f'{side}_relevant_metadata']
            cdli=md.get('externalNumbers.cdliNumber','')
            bmid=md.get('externalNumbers.bmIdNumber','')
            urls=[]
            if cdli:
                urls += [
                    f'https://cdli.mpiwg-berlin.mpg.de/artifacts/{cdli}',
                    f'https://cdli.ucla.edu/{cdli}',
                    f'https://cdli.earth/artifacts/{cdli}',
                ]
            if bmid:
                urls.append(f'https://www.britishmuseum.org/collection/object/{bmid}')
            probes=[fetch(s,u) for u in urls]
            rr['sides'].append({'side':side,'id':row[side],'cdli':cdli,'bm_id':bmid,'probes':probes})
            time.sleep(.15)
        results.append(rr)
    summary=[]
    for r in results:
        for sd in r['sides']:
            ok=[p for p in sd['probes'] if p.get('status')==200]
            img=sum(len(p.get('image_urls',[])) for p in ok)
            summary.append({'candidate_pair':f"{r['a']} + {r['b']}",'fragment':sd['id'],'http200_endpoints':len(ok),'image_urls_found':img,'successful_endpoints':[p['final_url'] for p in ok]})
    res={'fragments_probed':len(summary),'fragments_with_http200':sum(bool(x['http200_endpoints']) for x in summary),'fragments_with_image_urls':sum(bool(x['image_urls_found']) for x in summary),'summary':summary,'claim_boundary':'Endpoint/image discovery only. No geometric or philological join validation is inferred.'}
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'endpoint_details.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
