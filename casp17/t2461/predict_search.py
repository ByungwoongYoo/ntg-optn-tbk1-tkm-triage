#!/usr/bin/env python3
"""Predict T2461 monomers with the public ESMFold API and search Foldseek."""
from __future__ import annotations
import hashlib, json, os, statistics, time
from pathlib import Path
import requests

TARGET='T2461'
FULL='MHHHHHHGSMAIGIVELSSIAMGLKLADEMLKAADVKLLVSRPILPGKFLIILGGETEAIRKAIAVATEAAGSKLVRSALIEDIHPSVLPAISGINPVEERQAVGIVETESLEAAILAANAAVKGSNVTLVRIRMLSGITGKCYIVVAGDVDDVALAVVVAAEVAASRGKLIYAALIPRPHPAIWPLIVEG'
CORE=FULL[9:]
OUT=Path(os.getenv('OUT_DIR','artifact')); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({'User-Agent':'ByungwoongYoo-CASP17-T2461/1.0'})

def sh(b:bytes): return hashlib.sha256(b).hexdigest()

def esmfold(name,seq):
    url='https://api.esmatlas.com/foldSequence/v1/pdb/'
    last=None
    for attempt in range(1,5):
        r=S.post(url,data=seq,timeout=900); last=r
        if r.status_code==200 and b'ATOM' in r.content: break
        time.sleep(20*attempt)
    rec={'name':name,'sequence_length':len(seq),'status':last.status_code,'bytes':len(last.content),'sha256':sh(last.content)}
    if last.status_code!=200 or b'ATOM' not in last.content:
        rec['error']=last.text[:1000]; return rec,None
    path=OUT/f'{name}.pdb'; path.write_bytes(last.content); rec['file']=path.name
    ca=[]; bf=[]
    for line in last.text.splitlines():
        if line.startswith('ATOM'):
            try: bf.append(float(line[60:66]))
            except Exception: pass
            if line[12:16].strip()=='CA': ca.append(line)
    rec['ca_atoms']=len(ca); rec['mean_atom_bfactor']=statistics.fmean(bf) if bf else None
    return rec,path

def foldseek(query:Path,label:str,mode:str,databases:list[str]):
    url='https://search.foldseek.com/api/ticket'
    data=[('mode',mode)]+[('database[]',x) for x in databases]
    with query.open('rb') as f:
        r=S.post(url,files={'q':(query.name,f,'chemical/x-pdb')},data=data,timeout=300)
    ticket=r.json() if r.ok else {'status':'HTTP_ERROR','body':r.text[:1000]}
    rec={'label':label,'mode':mode,'databases':databases,'submission_status':r.status_code,'ticket':ticket}
    tid=ticket.get('id')
    if not tid: return rec
    status=ticket.get('status')
    for _ in range(120):
        if status in ('COMPLETE','ERROR','MAINTENANCE'): break
        time.sleep(5)
        p=S.get(f'https://search.foldseek.com/api/ticket/{tid}',timeout=120)
        try: ticket=p.json(); status=ticket.get('status')
        except Exception: status='POLL_PARSE_ERROR'; break
    rec['final_ticket']=ticket
    if status!='COMPLETE': return rec
    rr=S.get(f'https://search.foldseek.com/api/result/{tid}/0',timeout=300)
    rec['result_status']=rr.status_code; rec['result_bytes']=len(rr.content); rec['result_sha256']=sh(rr.content)
    (OUT/f'foldseek_{label}_{mode}_result.json').write_bytes(rr.content)
    dl=S.get(f'https://search.foldseek.com/api/result/download/{tid}',timeout=600)
    rec['download_status']=dl.status_code; rec['download_bytes']=len(dl.content); rec['download_sha256']=sh(dl.content)
    if dl.status_code==200: (OUT/f'foldseek_{label}_{mode}_download.bin').write_bytes(dl.content)
    return rec

def main():
    predictions=[]; paths={}
    for name,seq in [('T2461_full_esmfold',FULL),('T2461_core_esmfold',CORE)]:
        rec,p=esmfold(name,seq); predictions.append(rec)
        if p: paths[name]=p
    searches=[]
    for name,p in paths.items():
        searches.append(foldseek(p,name,'3diaa',['pdb100','cath50','afdb-swissprot']))
        searches.append(foldseek(p,name,'tmalign',['pdb100']))
    summary={'target':TARGET,'full_length':len(FULL),'core_start':10,'core_length':len(CORE),
             'predictions':predictions,'searches':searches,
             'claim_boundary':'Public ESMFold/Foldseek outputs are candidate models/templates, not experimental truth or official CASP submission.'}
    (OUT/'PREDICTION_SEARCH_SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
