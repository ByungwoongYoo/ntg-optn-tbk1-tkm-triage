#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,os,re,time
from pathlib import Path
import requests
from bs4 import BeautifulSoup
WANTED={'E2448':139,'E2449':92}
ESM='https://api.esmatlas.com/foldSequence/v1/pdb/'
def sha(b):return hashlib.sha256(b).hexdigest()
def text(html):return BeautifulSoup(html,'html.parser').get_text('\n',strip=True)
def parse_name(t):
 m=re.search(r'Target:\s*([A-Z]\d+)',t);return m.group(1) if m else None
def parse_seq(html):
 t=text(html);m=re.search(r'>[^\n]*\n?\s*\|?\s*([ACDEFGHIKLMNPQRSTVWYXOUZBJ\s]{8,1000})',t,re.I)
 if m:
  s=re.sub('[^A-Z]','',m.group(1).upper());
  if len(s)>=8:return s
 m=re.search(r'Sequence:\s*\([^\n]*\)\s*>[^|]*\|\s*([A-Z\s]{8,1000})',t,re.I)
 return re.sub('[^A-Z]','',m.group(1).upper()) if m else ''
def main():
 out=Path(os.environ.get('OUT_DIR','artifact/casp17_active_v2'));out.mkdir(parents=True,exist_ok=True)
 s=requests.Session();s.headers.update({'User-Agent':'CASP17-active-blind-freezer/2.0'})
 found={};inventory=[]
 for tid in range(252,281):
  try:
   r=s.get(f'https://predictioncenter.org/casp17/target.cgi?id={tid}&view=all',timeout=45);t=text(r.text);name=parse_name(t);inventory.append({'id':tid,'status':r.status_code,'name':name})
   if name in WANTED:
    seq=parse_seq(r.text);found[name]={'id':tid,'page_url':r.url,'page_sha256':sha(r.content),'expected_length':WANTED[name],'sequence':seq,'sequence_length':len(seq),'sequence_sha256':sha(seq.encode()) if seq else None}
  except Exception as e:inventory.append({'id':tid,'error':repr(e)})
 if set(found)!=set(WANTED):
  (out/'inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
  raise RuntimeError(f'Could not locate all wanted targets: found {found.keys()}')
 for name,rec in found.items():
  seq=rec['sequence'];rec['blind_prediction_frozen']=False
  if len(seq)!=rec['expected_length']:
   rec['skip_reason']='sequence-length mismatch';continue
  last=None
  for k in range(1,4):
   try:
    rr=s.post(ESM,data=seq,timeout=300);last=rr;rec.setdefault('esmfold_attempts',[]).append({'attempt':k,'status':rr.status_code,'bytes':len(rr.content)})
    if rr.ok and 'ATOM' in rr.text:break
   except Exception as e:rec.setdefault('esmfold_attempts',[]).append({'attempt':k,'error':repr(e)})
   time.sleep(5*k)
  if last is not None and last.ok and 'ATOM' in last.text:
   pdb=last.text;(out/f'{name}_ESMFold_blind.pdb').write_text(pdb,encoding='utf-8');rec['pdb_sha256']=sha(pdb.encode());rec['blind_prediction_frozen']=True
 (out/'inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
 result={'freeze_date_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'targets':found,'status':'ACTIVE_BLIND_PREDICTIONS_FROZEN' if any(x['blind_prediction_frozen'] for x in found.values()) else 'NO_PDB_FROZEN','claim_boundary':'Local hashes predate experimental truth but are not official CASP submissions. E2448/E2449 are active ensemble targets according to the official target list; official recognition requires a registered CASP17 group and timely submission through the proper gateway.'}
 (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
