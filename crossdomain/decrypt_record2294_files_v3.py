#!/usr/bin/env python3
from __future__ import annotations
import json,os,re,time
from pathlib import Path
import requests
FIELDS=['ciphertext','transcription','plaintext','cleartext','image','images','document','documents','file','files','ciphertext_file','transcription_file','plaintext_file']
BASE='https://de-crypt.org/decrypt-web/api/file/records/{field}/2294'
def main():
 out=Path(os.environ.get('OUT_DIR','artifact/decrypt2294'));out.mkdir(parents=True,exist_ok=True)
 s=requests.Session();s.headers.update({'User-Agent':'Mozilla/5.0 DECRYPT-public-API-file-audit/3.0'})
 rows=[]
 for field in FIELDS:
  url=BASE.format(field=field)
  try:
   r=s.get(url,timeout=45);rec={'field':field,'url':url,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content),'body_preview':r.text[:2500] if 'text' in r.headers.get('content-type','') or 'json' in r.headers.get('content-type','') else None}
   try:rec['json']=r.json()
   except Exception:pass
   rows.append(rec)
  except Exception as e:rows.append({'field':field,'url':url,'error':repr(e)})
  time.sleep(.05)
 result={'record_id':2294,'status':'DOCUMENTED_FILE_ENDPOINTS_PROBED','rows':rows,'claim_boundary':'Only documented unauthenticated GET file-info endpoints were tested. No login or access-control bypass is attempted.'}
 (out/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
