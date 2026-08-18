#!/usr/bin/env python3
"""Public-template discovery for CASP17 target T2461."""
from __future__ import annotations
import hashlib, json, os, re, time, urllib.error, urllib.request
from pathlib import Path

TARGET='T2461'
SEQ='MHHHHHHGSMAIGIVELSSIAMGLKLADEMLKAADVKLLVSRPILPGKFLIILGGETEAIRKAIAVATEAAGSKLVRSALIEDIHPSVLPAISGINPVEERQAVGIVETESLEAAILAANAAVKGSNVTLVRIRMLSGITGKCYIVVAGDVDDVALAVVVAAEVAASRGKLIYAALIPRPHPAIWPLIVEG'
OUT=Path(os.getenv('OUT_DIR','artifact')); OUT.mkdir(parents=True,exist_ok=True)
HEAD={'User-Agent':'ByungwoongYoo-CASP17-T2461/1.0'}

def req(url,data=None,timeout=180):
    h=dict(HEAD)
    if data is not None: h['Content-Type']='application/json'
    r=urllib.request.Request(url,data=data,headers=h)
    try:
        with urllib.request.urlopen(r,timeout=timeout) as x:
            return x.status,x.read(),dict(x.headers)
    except urllib.error.HTTPError as e:
        return e.code,e.read(),dict(e.headers)
    except Exception as e:
        return 0,repr(e).encode(),{}

def get_json(url):
    c,b,_=req(url)
    if c!=200: return {'_status':c,'_error':b.decode(errors='replace')[:500]}
    try: return json.loads(b)
    except Exception: return {'_status':c,'_error':'non-json'}

def sha(b): return hashlib.sha256(b).hexdigest()

def main():
    (OUT/'TARGET.fasta').write_text(f'>{TARGET}\n{SEQ}\n')
    q={'query':{'type':'terminal','service':'sequence','parameters':{
       'evalue_cutoff':10,'identity_cutoff':0.1,'sequence_type':'protein','value':SEQ}},
       'return_type':'polymer_entity','request_options':{'paginate':{'start':0,'rows':200},
       'results_verbosity':'verbose','scoring_strategy':'sequence'}}
    code,body,_=req('https://search.rcsb.org/rcsbsearch/v2/query',json.dumps(q).encode())
    raw=json.loads(body) if code==200 else {'error':body.decode(errors='replace'),'status':code}
    (OUT/'rcsb_sequence_search.json').write_text(json.dumps(raw,indent=2))
    hits=[]; candidates=[]
    for rank,h in enumerate(raw.get('result_set',[]),1):
        ident=h.get('identifier','')
        m=re.fullmatch(r'([0-9A-Za-z]{4})_(\d+)',ident)
        if not m: continue
        pdb=m.group(1).lower(); ent=m.group(2)
        entity=get_json(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{ent}')
        entry=get_json(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
        info=entity.get('entity_poly') or {}
        record={'rank':rank,'identifier':ident,'pdb':pdb,'entity_id':ent,'score':h.get('score'),
          'sequence_identity_context':h.get('matching_context'),'length':info.get('rcsb_sample_sequence_length'),
          'description':(entity.get('rcsb_polymer_entity') or {}).get('pdbx_description'),
          'title':(entry.get('struct') or {}).get('title'),'assemblies':[]}
        n=int((entry.get('rcsb_entry_info') or {}).get('assembly_count') or 0)
        for aid in range(1,n+1):
            a=get_json(f'https://data.rcsb.org/rest/v1/core/assembly/{pdb}/{aid}')
            aa=a.get('pdbx_struct_assembly') or {}; sy=a.get('rcsb_struct_symmetry') or []
            ar={'assembly_id':aid,'oligomeric_count':aa.get('oligomeric_count'),
                'details':aa.get('details'),'symmetry':[{'symbol':s.get('symbol'),
                'oligomeric_state':s.get('oligomeric_state'),'stoichiometry':s.get('stoichiometry')} for s in sy]}
            record['assemblies'].append(ar)
            try: olig=int(ar['oligomeric_count'] or 0)
            except Exception: olig=0
            if olig in (12,24,60) or any(str(s.get('oligomeric_state','')).startswith(('12','24','60')) for s in ar['symmetry']):
                candidates.append({**{k:record[k] for k in ('rank','identifier','pdb','entity_id','score','length','description','title')},'assembly':ar})
        hits.append(record); time.sleep(.02)
    (OUT/'template_hits.json').write_text(json.dumps(hits,indent=2))
    (OUT/'oligomeric_candidates.json').write_text(json.dumps(candidates,indent=2))

    probes=[]
    urls=[
      'https://predictioncenter.org/download_area/CASP17/predictions/oligo/T2461.tar.gz',
      'https://predictioncenter.org/download_area/CASP17/predictions/regular/T2461.tar.gz',
      'https://predictioncenter.org/download_area/CASP17/predictions/oligo/',
      'https://predictioncenter.org/download_area/CASP17/predictions/regular/',
      'https://predictioncenter.org/casp17/target.cgi?target=T2461&view=template']
    for i,u in enumerate(urls):
        c,b,h=req(u)
        name=f'official_probe_{i}' + ('.tar.gz' if u.endswith('.tar.gz') else '.html')
        if c==200: (OUT/name).write_bytes(b)
        probes.append({'url':u,'status':c,'bytes':len(b),'sha256':sha(b) if c==200 else None,
                       'content_type':h.get('Content-Type'),'file':name if c==200 else None,
                       'error':None if c==200 else b.decode(errors='replace')[:500]})
    (OUT/'official_archive_probes.json').write_text(json.dumps(probes,indent=2))

    downloads=[]
    for x in candidates[:20]:
        p=x['pdb']; aid=x['assembly']['assembly_id']; u=f'https://files.rcsb.org/download/{p}-assembly{aid}.cif'
        c,b,_=req(u); fn=f'template_{p}_assembly{aid}.cif'
        if c==200: (OUT/fn).write_bytes(b)
        downloads.append({'pdb':p,'assembly':aid,'status':c,'bytes':len(b),'sha256':sha(b) if c==200 else None,'file':fn if c==200 else None})
    (OUT/'template_downloads.json').write_text(json.dumps(downloads,indent=2))
    summary={'target':TARGET,'sequence_length':len(SEQ),'rcsb_hits':len(hits),
             'oligomeric_candidates':len(candidates),'top_candidates':candidates[:20],'probes':probes}
    (OUT/'DISCOVERY_SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
