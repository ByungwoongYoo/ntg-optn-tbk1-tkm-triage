#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path

BONUS={'T1':2,'T2':2,'T3':1,'T4':0,'T5':0,'T6':0,'T7':2,'T8':0}

def load(path):
    try:return json.loads(Path(path).read_text(encoding='utf-8'))
    except Exception:return None

def find(root,name):
    xs=list(Path(root).rglob(name));return xs[0] if xs else None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--out-dir',required=True);a=ap.parse_args()
    root=Path(a.root);out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    rows=[]
    def add(t,title,stage,control,candidate,note,details=None):
        score=10*stage+3*int(control)+2*int(candidate)+BONUS[t]
        rows.append({'track':t,'title':title,'evidence_stage':stage,'control_pass':bool(control),'candidate_ready':bool(candidate),'medical_bonus':BONUS[t],'score':score,'note':note,'details':details or {}})

    # T1 eBL/Nineveh
    p=find(root,'T1_EBL_RESULT.json') or find(root,'RESULT_EBL.json')
    d=load(p) if p else None
    if d:
        med=(d.get('retrieval_medical_sources') or {}); n=med.get('n_sources') or 0; top100=med.get('top100')
        control=bool(n and top100 is not None and top100>0)
        cand=(d.get('candidate_count_saved') or 0)>0
        add('T1','Nineveh/eBL medical fragment joins',2 if control else 1,control,cand,'Text/sign benchmark only; physical/image join verification still required.',d)
    else:add('T1','Nineveh/eBL medical fragment joins',0,False,False,'No result artifact.')

    # T2 East Asian lost medical texts + Yanxia
    g=find(root,'T2_GANGIBANG_RESULT.json'); y=find(root,'T2_YANXIA_RESULT.json');gd=load(g) if g else None;yd=load(y) if y else None
    control=bool(gd and (gd.get('explicit_source_blocks') or 0)>0)
    ycand=bool(yd and (yd.get('explicit_source_blocks') or 0)>0 and (yd.get('formula_name_candidate_count') or 0)>0)
    stage=2 if control else (1 if yd or gd else 0)
    add('T2','Uibangyuchwi lost texts: Jianqifang control + Yanxia Shengxiaofang',stage,control,ycand,'Yanxia candidates are source-attributed text candidates, not novelty claims.',{'gangibang':gd,'yanxia':yd})

    # T3 DECRYPT
    p=find(root,'T3_DECRYPT_RESULT.json');d=load(p) if p else None
    accessible=bool(d and (d.get('unique_ids_first_pages') or 0)>0);cand=bool(d and len(d.get('medical_alchemical_keyword_records') or [])>0)
    add('T3','Historical medical/alchemical ciphers',1 if accessible else 0,False,cand,'No exact plaintext/key control in this run; public-target triage only.',d or {})

    # T4 u4
    p=find(root,'T4_U4_RESULT.json');d=load(p) if p else None
    control=bool(d and d.get('controls_reproduced'));dec=bool(d and d.get('decisive_candidate'))
    stage=3 if control and dec else (2 if control else (1 if d else 0))
    add('T4','Positive implicational logic u4',stage,control,dec,'Any decisive u4 status remains a candidate until independent checking; timeout is no result.',d or {})

    # T5 RePAIR from access probe
    p=find(root,'T0_ACCESS_RESULT.json');access=load(p) if p else None
    rep=((access or {}).get('tracks') or {}).get('T5_repair') or {}
    if rep.get('status')=='METADATA_ONLY_LEGAL_GATE':add('T5','Pompeii RePAIR fresco joins',0,False,False,'Blocked at custom-license/user-acceptance gate; dataset not downloaded.',rep)
    elif rep:add('T5','Pompeii RePAIR fresco joins',1,False,False,'Metadata/data probe completed.',rep)
    else:add('T5','Pompeii RePAIR fresco joins',0,False,False,'No result artifact.')

    # T6 Vesuvius
    p=find(root,'T6_VESUVIUS_RESULT.json');d=load(p) if p else None
    ok=bool(d and str(d.get('status','')).upper() not in ('FAILED',''))
    add('T6','Vesuvius unopened-scroll reading',1 if ok else 0,False,False,'Current public-volume random access is not ink/text recovery; verified-region smoke remains required.',d or {})

    # T7 CASP17
    p=find(root,'T7_CASP_RESULT.json');d=load(p) if p else None
    ok=bool(d and str(d.get('status','')).upper() not in ('FAILED',''))
    candidate=bool(d and (d.get('frozen_candidates') or d.get('candidate') or d.get('selected_target')))
    add('T7','CASP17 blind structure prediction',1 if ok else 0,False,candidate,'Not an official blind result unless prediction is registered/submitted before truth release.',d or {})

    # T8 B2
    p=find(root,'T8_B2_RESULT.json');d=load(p) if p else None
    witness=bool(d and d.get('verified'))
    stage=4 if witness else (1 if d else 0)
    add('T8','B2[2] subset of Z_100, 13-vs-14',stage,False,witness,'A verified 14-set is decisive. UNKNOWN/near-miss is no negative result.',d or {})

    rows.sort(key=lambda x:(x['score'],x['evidence_stage']),reverse=True)
    result={'scoring_frozen_before_run':True,'formula':'10*evidence_stage + 3*control_pass + 2*candidate_ready + medical_bonus','winner':rows[0]['track'] if rows else None,'ranking':rows,
            'claim_boundary':'Tournament rank is a triage decision, not proof that the top track is solved. Only evidence stage 4 is a verified witness/certificate/official blind result.'}
    (out/'TOURNAMENT_RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    md=['# Cross-domain tournament result','',f"Winner by frozen triage score: **{result['winner']}**",'', '|Rank|Track|Stage|Control|Candidate|Score|Note|','|---:|---|---:|---|---|---:|---|']
    for i,r in enumerate(rows,1):md.append(f"|{i}|{r['track']} {r['title']}|{r['evidence_stage']}|{'Y' if r['control_pass'] else 'N'}|{'Y' if r['candidate_ready'] else 'N'}|{r['score']}|{r['note']}|")
    md += ['', '> No stage-3 candidate is called a discovery. Stage 4 requires the track-specific independent certificate.']
    (out/'TOURNAMENT_RESULTS.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
