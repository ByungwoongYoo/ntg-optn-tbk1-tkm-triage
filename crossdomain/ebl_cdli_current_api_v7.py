#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re,time
from pathlib import Path
import requests

ART_RE=re.compile(r'/artifacts/(\d+)')
ASSET_WORDS=('image','photo','picture','media','resource','external','iiif','thumbnail','lineart','scan','url','uri')

def harvest(o,p=''):
    out={}
    if isinstance(o,dict):
        for k,v in o.items():
            q=f'{p}.{k}' if p else k
            if any(w in k.lower() for w in ASSET_WORDS): out[q]=v
            if isinstance(v,(dict,list)): out.update(harvest(v,q))
    elif isinstance(o,list):
        for i,v in enumerate(o[:100]):
            if isinstance(v,(dict,list)): out.update(harvest(v,f'{p}[{i}]'))
    return out

def req(sess,url,headers=None):
    try:
        r=sess.get(url,headers=headers or {},timeout=45,allow_redirects=True)
        rec={'url':url,'status':r.status_code,'final_url':r.url,'content_type':r.headers.get('content-type',''),'bytes':len(r.content)}
        txt=r.text[:1000000] if ('text' in rec['content_type'] or 'json' in rec['content_type'] or 'html' in rec['content_type']) else ''
        rec['artifact_ids_in_body']=list(dict.fromkeys(ART_RE.findall(txt)))[:20]
        if 'json' in rec['content_type']:
            try:
                j=r.json(); rec['json']=j; rec['asset_fields']=harvest(j)
            except Exception as e: rec['json_error']=repr(e)
        return rec
    except Exception as e:
        return {'url':url,'error':repr(e)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--asset-metadata',required=True); ap.add_argument('--out-dir',required=True); ap.add_argument('--top',type=int,default=3); a=ap.parse_args()
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(Path(a.asset_metadata).read_text(encoding='utf-8'))[:a.top]
    sess=requests.Session(); sess.headers.update({'User-Agent':'Mozilla/5.0 open-research-cuneiform/1.1'})
    frags=[]; seen=set()
    known_internal={'P398727':'628988','P426107':'812220','P401572':'631811'}
    for row in rows:
        for side in ('a','b'):
            frag=row[side]
            if frag in seen: continue
            seen.add(frag)
            md=row[f'{side}_relevant_metadata']; p=md.get('externalNumbers.cdliNumber','')
            probes=[]
            if p:
                probes.append(req(sess,f'https://cdli.earth/{p}'))
                probes.append(req(sess,f'https://cdli.earth/{p}',{'Accept':'application/json'}))
                probes.append(req(sess,f'https://cdli.earth/artifacts/{p}.json'))
                probes.append(req(sess,f'https://cdli.earth/artifacts/{p}',{'Accept':'application/json'}))
                iid=known_internal.get(p)
                if iid:
                    probes.append(req(sess,f'https://cdli.earth/artifacts/{iid}.json'))
                    probes.append(req(sess,f'https://cdli.earth/artifacts/{iid}',{'Accept':'application/json'}))
                    probes.append(req(sess,f'https://cdli.earth/artifactsExternalResources/{iid}.json'))
            ids=[]
            for pr in probes:
                ids += pr.get('artifact_ids_in_body',[])
                m=ART_RE.search(pr.get('final_url',''))
                if m: ids.append(m.group(1))
            for iid in list(dict.fromkeys(ids))[:5]:
                u=f'https://cdli.earth/artifacts/{iid}.json'
                if not any(x.get('url')==u for x in probes): probes.append(req(sess,u))
                u2=f'https://cdli.earth/artifactsExternalResources/{iid}.json'
                if not any(x.get('url')==u2 for x in probes): probes.append(req(sess,u2))
            frags.append({'fragment':frag,'p_number':p,'known_internal_id':known_internal.get(p),'probes':probes})
            time.sleep(.1)
    summary=[]
    for f in frags:
        js=[p for p in f['probes'] if p.get('status')==200 and 'json' in p.get('content_type','')]
        af=sum(len(p.get('asset_fields',{})) for p in js)
        summary.append({'fragment':f['fragment'],'p_number':f['p_number'],'known_internal_id':f['known_internal_id'],'json_200_endpoints':len(js),'asset_field_count':af,'successful_json_urls':[p['final_url'] for p in js]})
    res={'fragments_probed':len(frags),'fragments_with_json_api':sum(bool(x['json_200_endpoints']) for x in summary),'fragments_with_asset_fields':sum(bool(x['asset_field_count']) for x in summary),'summary':summary,'claim_boundary':'Current CDLI API-route discovery only. Asset metadata can enable later image retrieval; no join claim is made.'}
    (out/'RESULT.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'api_probe_details.json').write_text(json.dumps(frags,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
