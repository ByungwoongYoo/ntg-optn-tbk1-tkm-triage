#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, os, re, time
from pathlib import Path
import requests
from bs4 import BeautifulSoup

BASE='https://predictioncenter.org/casp17/target.cgi?id={id}&view=all'
TARGETS={
 'E2446':(250,107),'E2447':(251,97),'E2458':(267,131),'E2459':(268,117),'E2460':(269,205),'T2463':(272,16),'T2464':(273,17),
}
ESM='https://api.esmatlas.com/foldSequence/v1/pdb/'

def sha(b):return hashlib.sha256(b).hexdigest()
def extract_sequence(html):
    text=BeautifulSoup(html,'html.parser').get_text('\n',strip=True);soup=BeautifulSoup(html,'html.parser')
    candidates=[tag.get_text('\n',strip=True) for tag in soup.find_all(['pre','textarea','tt','code'])]+[text]
    for chunk in candidates:
        m=re.search(r'>[^\n]*\n([A-Z\n\r ]{8,10000})',chunk)
        if m:
            seq=re.sub('[^A-Z]','',m.group(1));
            if len(seq)>=8:return seq
        m=re.search(r'(?:Target Sequence|Sequence)\s*[:\n]\s*([ACDEFGHIKLMNPQRSTVWYXOUZBJ\s]{8,10000})',chunk,re.I)
        if m:
            seq=re.sub('[^A-Z]','',m.group(1).upper());
            if len(seq)>=8:return seq
    for m in re.finditer(r'([ACDEFGHIKLMNPQRSTVWYX]{8,1000})',html):
        s=m.group(1)
        if len(s)>=8 and not any(w in s for w in ('TARGET','SEQUENCE','BACKGROUND')):return s
    return ''

def pdb_stats(text):
    seen=set()
    for line in text.splitlines():
        if line.startswith(('ATOM  ','HETATM')) and line[12:16].strip()=='CA':seen.add((line[21:22],line[22:26].strip(),line[26:27]))
    return {'ca_atoms':len(seen),'unique_residues':len(seen),'pdb_bytes':len(text.encode())}

def main():
    out=Path(os.environ.get('OUT_DIR','artifact/casp17_blind'));out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers.update({'User-Agent':'CASP17-public-blind-feasibility/1.0'})
    results={}
    for target,(tid,expected_len) in TARGETS.items():
        page=s.get(BASE.format(id=tid),timeout=60);rec={'page_status':page.status_code,'page_url':page.url,'page_sha256':sha(page.content),'expected_length_from_target_list':expected_len}
        (out/f'{target}_target.html').write_bytes(page.content)
        seq=extract_sequence(page.text);rec['sequence']=seq;rec['sequence_length']=len(seq);rec['sequence_length_matches_target_list']=(len(seq)==expected_len);rec['sequence_sha256']=sha(seq.encode()) if seq else None
        if seq and len(seq)==expected_len:
            last=None
            for attempt in range(1,4):
                try:
                    r=s.post(ESM,data=seq,timeout=300);last=r;rec.setdefault('esmfold_attempts',[]).append({'attempt':attempt,'http_status':r.status_code,'bytes':len(r.content)})
                    if r.ok and ('ATOM' in r.text or 'HETATM' in r.text):break
                except Exception as e:rec.setdefault('esmfold_attempts',[]).append({'attempt':attempt,'error':repr(e)})
                time.sleep(5*attempt)
            if last is not None and last.ok and ('ATOM' in last.text or 'HETATM' in last.text):
                pdb=last.text;(out/f'{target}_ESMFold_blind.pdb').write_text(pdb,encoding='utf-8');rec['pdb_sha256']=sha(pdb.encode());rec['pdb_stats']=pdb_stats(pdb);rec['blind_prediction_frozen']=rec['pdb_stats']['ca_atoms']==expected_len
            else:
                if last is not None:(out/f'{target}_ESMFold_error.txt').write_bytes(last.content[:200000])
                rec['blind_prediction_frozen']=False
        else:
            rec['blind_prediction_frozen']=False;rec['skip_reason']='sequence extraction failed or target-list length mismatch'
        results[target]=rec;print(target,rec,flush=True);time.sleep(.5)
    result={'freeze_date_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'source':'official CASP17 target pages + length cross-check from official target list','predictor':'public ESMFold API, no target-specific tuning; HTTP failures retried only','targets':results,
            'claim_boundary':'These are pre-outcome blind candidate structures, not official CASP submissions. Official scoring requires Prediction Center account/group registration and submission before the target deadline. No structural accuracy claim is possible until experimental targets are released.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    manifest={}
    for p in sorted(out.glob('*')):
        if p.is_file() and p.name!='MANIFEST.json':manifest[p.name]={'bytes':p.stat().st_size,'sha256':sha(p.read_bytes())}
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
