#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import requests

TYPES=[
 ('photo','https://cdli.ucla.edu/dl/photo/{p}.jpg'),
 ('photo_detail','https://cdli.ucla.edu/dl/photo/{p}_d.jpg'),
 ('photo_edge','https://cdli.ucla.edu/dl/photo/{p}_e.jpg'),
 ('lineart','https://cdli.ucla.edu/dl/lineart/{p}_l.jpg'),
 ('lineart_detail','https://cdli.ucla.edu/dl/lineart/{p}_ld.jpg'),
 ('lineart_side','https://cdli.ucla.edu/dl/lineart/{p}_ls.jpg'),
 ('thumb','https://cdli.ucla.edu/dl/tn_photo/{p}.jpg'),
]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--asset-metadata',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--top',type=int,default=3); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(Path(a.asset_metadata).read_text(encoding='utf-8'))[:a.top]
    sess=requests.Session(); sess.headers.update({'User-Agent':'Mozilla/5.0 open-research-cuneiform/1.0'})
    seen=set(); rec=[]
    for row in rows:
        for side in ('a','b'):
            frag=row[side]
            if frag in seen: continue
            seen.add(frag)
            md=row[f'{side}_relevant_metadata']; p=md.get('externalNumbers.cdliNumber','')
            rr={'fragment':frag,'p_number':p,'files':[]}
            if p:
                for typ,tmpl in TYPES:
                    url=tmpl.format(p=p)
                    try:
                        r=sess.get(url,timeout=60,allow_redirects=True)
                        ok=r.status_code==200 and r.headers.get('content-type','').lower().startswith('image/') and len(r.content)>1000
                        item={'type':typ,'url':url,'status':r.status_code,'content_type':r.headers.get('content-type',''),'bytes':len(r.content),'ok_image':ok}
                        if ok:
                            ext='.jpg'; fn=f'{frag.replace(".","_")}_{p}_{typ}{ext}'; (out/fn).write_bytes(r.content)
                            item['file']=fn; item['sha256']=hashlib.sha256(r.content).hexdigest()
                        rr['files'].append(item)
                    except Exception as e:
                        rr['files'].append({'type':typ,'url':url,'error':repr(e),'ok_image':False})
            rec.append(rr)
    res={'fragments':len(rec),'fragments_with_any_image':sum(any(f.get('ok_image') for f in r['files']) for r in rec),'images_downloaded':sum(sum(bool(f.get('ok_image')) for f in r['files']) for r in rec),'records':rec,'claim_boundary':'Images are fetched from the CDLI digital-library file patterns documented by CDLI. Image availability enables visual edge/orientation screening but is not itself evidence of a join.'}
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
