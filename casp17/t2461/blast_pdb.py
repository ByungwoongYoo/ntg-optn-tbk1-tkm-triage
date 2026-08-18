#!/usr/bin/env python3
"""Public NCBI BLASTP search of CASP17 T2461 against the PDB database."""
import hashlib,json,os,re,time,xml.etree.ElementTree as ET
from pathlib import Path
import requests
FULL='MHHHHHHGSMAIGIVELSSIAMGLKLADEMLKAADVKLLVSRPILPGKFLIILGGETEAIRKAIAVATEAAGSKLVRSALIEDIHPSVLPAISGINPVEERQAVGIVETESLEAAILAANAAVKGSNVTLVRIRMLSGITGKCYIVVAGDVDDVALAVVVAAEVAASRGKLIYAALIPRPHPAIWPLIVEG'
CORE=FULL[9:]
OUT=Path(os.getenv('OUT_DIR','artifact'));OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session();S.headers.update({'User-Agent':'ByungwoongYoo CASP17 T2461; yoonge3@gmail.com'})
URL='https://blast.ncbi.nlm.nih.gov/Blast.cgi'

def run(label,seq):
    p={'CMD':'Put','PROGRAM':'blastp','DATABASE':'pdb','QUERY':seq,'EXPECT':'1000','HITLIST_SIZE':'200','FILTER':'L'}
    r=S.post(URL,data=p,timeout=300); text=r.text
    (OUT/f'{label}_put_response.txt').write_text(text)
    m=re.search(r'RID\s*=\s*([^\s]+)',text); e=re.search(r'RTOE\s*=\s*(\d+)',text)
    if not m:return {'label':label,'error':'RID not found','status':r.status_code,'preview':text[:1000]}
    rid=m.group(1); wait=max(10,int(e.group(1)) if e else 20); time.sleep(min(wait,60))
    status='WAITING'
    for _ in range(90):
        x=S.get(URL,params={'CMD':'Get','RID':rid,'FORMAT_OBJECT':'SearchInfo'},timeout=120); st=x.text
        if 'Status=READY' in st: status='READY';break
        if 'Status=FAILED' in st or 'Status=UNKNOWN' in st: status=st.strip();break
        time.sleep(10)
    rec={'label':label,'rid':rid,'poll_status':status}
    if status!='READY':return rec
    x=S.get(URL,params={'CMD':'Get','RID':rid,'FORMAT_TYPE':'XML'},timeout=600)
    (OUT/f'{label}_blast.xml').write_bytes(x.content); rec['xml_status']=x.status_code;rec['xml_bytes']=len(x.content);rec['xml_sha256']=hashlib.sha256(x.content).hexdigest()
    hits=[]
    try:
        root=ET.fromstring(x.content)
        for h in root.findall('.//Hit')[:200]:
            hs=h.find('Hit_hsps/Hsp')
            hits.append({'id':h.findtext('Hit_id'),'def':h.findtext('Hit_def'),'accession':h.findtext('Hit_accession'),
             'length':int(h.findtext('Hit_len') or 0),'evalue':float(hs.findtext('Hsp_evalue') or 'inf'),
             'bit_score':float(hs.findtext('Hsp_bit-score') or 0),'identity':int(hs.findtext('Hsp_identity') or 0),
             'align_len':int(hs.findtext('Hsp_align-len') or 0),'q_from':int(hs.findtext('Hsp_query-from') or 0),
             'q_to':int(hs.findtext('Hsp_query-to') or 0),'h_from':int(hs.findtext('Hsp_hit-from') or 0),'h_to':int(hs.findtext('Hsp_hit-to') or 0)})
    except Exception as z: rec['parse_error']=repr(z)
    rec['hits']=hits;(OUT/f'{label}_hits.json').write_text(json.dumps(hits,indent=2));return rec

def main():
    out=[run('T2461_core',CORE),run('T2461_full',FULL)]
    (OUT/'BLAST_SUMMARY.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
