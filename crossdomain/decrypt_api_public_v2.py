#!/usr/bin/env python3
from __future__ import annotations
import json,os,re,time
from pathlib import Path
import requests
IDS=[990,991,993,994,995,1001,1002,1003,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1057,1072,1073,1074,1076,1325,1408,1469,1470,1874,1889,1890,1893,1894,1953,2078,2294,8342,8345,8755,8756,9322,9323,9403,9407,9408,9410,9411,9412,9413,9414,9415,9416,9417,9424,9427,9652,9772,9787,9952,9957,9959,9960,9962,9963,9966,9967,9970,10170,10171,10172,10173,10174,10175,10176,10177,10178,10179,10180,10181,10198]
BASE='https://de-crypt.org/decrypt-web'
KEYWORDS=re.compile(r'alchem|medic|physic|recipe|remed|pharmac|doctor|surgeon|disease|health|iatro',re.I)
def flatten(x,p=''):
    out=[]
    if isinstance(x,dict):
        for k,v in x.items():out+=flatten(v,p+'.'+str(k) if p else str(k))
    elif isinstance(x,list):
        for i,v in enumerate(x):out+=flatten(v,p+f'[{i}]')
    else: out.append((p,x))
    return out
def main():
    out=Path(os.environ.get('OUT_DIR','artifact/decrypt_api_v2'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 public-DECRYPT-API-research/2.0'})
    records=[]
    for rid in IDS:
        url=f'{BASE}/api/view/records/{rid}'
        try:
            r=s.get(url,timeout=45); rec={'id':rid,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content)}
            if r.ok:
                try:
                    j=r.json(); rec['json']=j; flat=flatten(j); rec['interesting_fields']=[{'path':p,'value':str(v)[:1500]} for p,v in flat if re.search(r'trans|plain|cipher|status|language|title|comment|document|file|author|date',p,re.I)][:100]
                    blob=' '.join(str(v) for _,v in flat);rec['medical_alchemical']=bool(KEYWORDS.search(blob));rec['keyword_context']=blob[max(0,(KEYWORDS.search(blob).start() if KEYWORDS.search(blob) else 0)-200):][:1200] if rec['medical_alchemical'] else ''
                except Exception as e:rec['parse_error']=repr(e);rec['body']=r.text[:3000]
            records.append(rec)
        except Exception as e:records.append({'id':rid,'error':repr(e)})
        time.sleep(.05)
    openjson=[x for x in records if x.get('json')]
    fields=set()
    for x in openjson:
        for f in x.get('interesting_fields',[]):fields.add(f['path'])
    result={'status':'PUBLIC_API_AUDIT_COMPLETE','ids_requested':len(IDS),'json_records_returned':len(openjson),'medical_alchemical_records':[x['id'] for x in openjson if x.get('medical_alchemical')],
            'interesting_field_paths':sorted(fields),'candidate_records':[{'id':x['id'],'interesting_fields':x.get('interesting_fields',[]),'medical_alchemical':x.get('medical_alchemical')} for x in openjson[:40]],
            'claim_boundary':'Only unauthenticated GET endpoints documented by DECRYPT were queried. File/transcription availability is treated as public only if the API actually returns it without credentials; no authentication bypass is attempted.'}
    (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');(out/'records_raw.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
