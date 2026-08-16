#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, re, subprocess, sys, time
from pathlib import Path
from urllib.parse import quote
import requests

OUT=Path(os.environ.get('OUT_DIR','artifact/probe')); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({'User-Agent':'BioOpenProblemProbe/1.0 (independent research; yoonge3@gmail.com)'})

def req(url, method='GET', **kw):
    r=S.request(method,url,timeout=kw.pop('timeout',60),allow_redirects=True,**kw)
    return r

def check_url(url, range_bytes=None):
    try:
        headers={'Range':f'bytes=0-{range_bytes-1}'} if range_bytes else {}
        r=req(url,headers=headers,stream=True)
        n=0
        if r.ok:
            for chunk in r.iter_content(65536):
                n+=len(chunk)
                if range_bytes and n>=range_bytes: break
        return {'url':url,'status':r.status_code,'ok':r.ok,'content_length':r.headers.get('content-length'),'content_range':r.headers.get('content-range'),'bytes_read':n,'final_url':r.url}
    except Exception as e:
        return {'url':url,'ok':False,'error':repr(e)}

res={}
# ProteinGym reference and bulk resources
pg_ref='https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv'
r=req(pg_ref); r.raise_for_status(); (OUT/'DMS_substitutions.csv').write_bytes(r.content)
rows=list(csv.DictReader(r.text.splitlines()))
res['proteingym']={
 'reference_ok':True,'n_assays':len(rows),'n_human_assays':sum(str(x.get('taxon','')).lower() in {'human','homo sapiens','9606'} or 'human' in str(x.get('taxon','')).lower() for x in rows),
 'dms_zip':check_url('https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip',range_bytes=1024*1024),
 'score_zip':check_url('https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_substitutions_scores.zip',range_bytes=1024*1024),
}
# Arc Virtual Cell public bucket API and documented dataset page
arc_api='https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/2025/&maxResults=100'
try:
    ar=req(arc_api); arc={'status':ar.status_code,'ok':ar.ok}
    if ar.ok:
        j=ar.json(); arc['n_listed']=len(j.get('items',[])); arc['objects']=[{'name':x.get('name'),'size':x.get('size')} for x in j.get('items',[])[:30]]
    else: arc['body']=ar.text[:500]
except Exception as e: arc={'ok':False,'error':repr(e)}
res['virtual_cell']={'bucket_probe':arc,'documented_training_size_gb':18.55,'documented_train_perturbations':150,'documented_validation':50,'documented_test':100}
# CAGI access gates
res['cagi_splicing']={'challenge_page':check_url('https://genomeinterpretation.org/cagi7-splicing.html',range_bytes=200000),'data_access':'registered Synapse users only','n_variants':9133,'n_labeled_sample':1257}
res['cagi_missense']={'challenge_page':check_url('https://genomeinterpretation.org/cagi7-annotate-all-missense.html',range_bytes=200000),'data_zip':check_url('https://genomeinterpretation.org/download/dbNSFP5.1_nsSNV.zip',range_bytes=1024*1024),'n_variants':82198516}
res['cagi_lentimpra']={'challenge_page':check_url('https://genomeinterpretation.org/cagi7-lenti-mpra.html',range_bytes=200000),'data_access':'Synapse/registration expected'}
res['rare_genomes']={'challenge_page':check_url('https://www.genomeinterpretation.org/cagi7-rgp-cram.html',range_bytes=200000),'data_access':'encrypted; institutional signature required'}
# CAFA status
res['cafa']={'site':check_url('https://biofunctionprediction.org/',range_bytes=200000),'note':'future-annotation benchmark, but no verified active 2026 submission round in probe'}
# CAMDA AMR agreement gate
res['camda_amr']={'site':check_url('https://bipress.boku.ac.at/camda2026/the-camda-contest-challenges/',range_bytes=200000),'data_access':'download agreement required','train_n':4800,'test_n':1500}
# CAMI III current public pages and toy data page
res['cami_iii']={'site':check_url('https://cami-challenge.org/cami-iii-challenges/',range_bytes=200000),'toy_page':check_url('https://cami-challenge.org/datasets/toy-human-gut/',range_bytes=200000),'status':'current challenge; toy set has gold standard; challenge set may require registration'}

# Transparent feasibility scores (0-3 per axis)
ratings={
 'ProteinGym_missense':{'objective_truth':3,'open_data':3,'free_compute_feasibility':2,'headroom':3,'immediate_benchmark':3},
 'Virtual_cell':{'objective_truth':3,'open_data':2 if arc.get('ok') else 1,'free_compute_feasibility':0,'headroom':3,'immediate_benchmark':2},
 'CAGI_splicing':{'objective_truth':3,'open_data':1,'free_compute_feasibility':2,'headroom':3,'immediate_benchmark':1},
 'CAGI_lentiMPRA':{'objective_truth':3,'open_data':1,'free_compute_feasibility':2,'headroom':3,'immediate_benchmark':1},
 'CAFA':{'objective_truth':3,'open_data':2,'free_compute_feasibility':1,'headroom':3,'immediate_benchmark':1},
 'CAMDA_AMR':{'objective_truth':3,'open_data':1,'free_compute_feasibility':1,'headroom':2,'immediate_benchmark':1},
 'Rare_genomes':{'objective_truth':3,'open_data':0,'free_compute_feasibility':0,'headroom':3,'immediate_benchmark':0},
 'CAMI_III':{'objective_truth':3,'open_data':2,'free_compute_feasibility':1,'headroom':2,'immediate_benchmark':3},
}
for v in ratings.values(): v['total']=sum(v.values())
res['ratings']=ratings
res['ranking']=sorted([{'track':k,**v} for k,v in ratings.items()],key=lambda x:x['total'],reverse=True)
(OUT/'probe_results.json').write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding='utf-8')
with (OUT/'PROBE_REPORT.md').open('w',encoding='utf-8') as f:
    f.write('# Bio/medical open-problem smoke probe\n\n')
    f.write('|Rank|Track|Total /15|Objective truth|Open data|Free compute|Headroom|Immediate benchmark|\n|---:|---|---:|---:|---:|---:|---:|---:|\n')
    for i,x in enumerate(res['ranking'],1):
        f.write(f"|{i}|{x['track']}|{x['total']}|{x['objective_truth']}|{x['open_data']}|{x['free_compute_feasibility']}|{x['headroom']}|{x['immediate_benchmark']}|\n")
    f.write('\n## Decision\n\nProteinGym missense prediction advances to an actual benchmark experiment because it combines public labels, official scoring, moderate compute, and substantial headroom. Access- or compute-gated tracks remain documented but do not advance by assumption.\n')
print(json.dumps(res['ranking'],indent=2))
