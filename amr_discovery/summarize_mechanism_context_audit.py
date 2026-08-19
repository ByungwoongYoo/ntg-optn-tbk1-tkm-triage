#!/usr/bin/env python3
"""Integrate direct sequence recall, AMRFinder, Kleborate and literature audits.

The resulting `advance_to_independent_validation` flag is a computational triage
status, not a novelty or causality claim.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd


def args():
    p=argparse.ArgumentParser();p.add_argument('--classification',required=True);p.add_argument('--support',required=True);p.add_argument('--sequence-dir',required=True);p.add_argument('--amrfinder',default='');p.add_argument('--kleborate',default='');p.add_argument('--literature',default='');p.add_argument('--context-amrfinder',default='');p.add_argument('--out',required=True);return p.parse_args()

def read(path,sep=None):
    p=Path(path) if path else None
    if not p or not p.exists() or p.stat().st_size==0:return pd.DataFrame()
    try:return pd.read_csv(p,sep=sep or ('\t' if p.suffix.lower() in ['.tsv','.txt'] else ','),dtype=str)
    except:return pd.DataFrame()

def assembly_col(df):
    for c in ['assembly_ID','assembly','Assembly','sample','strain','input_file','file']:
        if c in df.columns:return c
    return None

def norm_acc(x):
    s=Path(str(x)).name
    for suf in ['.fna.gz','.fasta.gz','.fa.gz','.fna','.fasta','.fa']:
        if s.endswith(suf):s=s[:-len(suf)];break
    m=re.search(r'(GC[AF]_[0-9]+\.[0-9]+)',s)
    return m.group(1) if m else s

def flag_mcr(df):
    if df.empty:return pd.DataFrame(columns=['assembly_ID','mcr_detected','colistin_keyword_detected'])
    ac=assembly_col(df);text=df.astype(str).agg(' '.join,axis=1).str.lower();x=pd.DataFrame({'assembly_ID':df[ac].map(norm_acc) if ac else 'UNKNOWN','mcr':text.str.contains(r'\bmcr[-_ ]?[0-9]',regex=True),'col':text.str.contains('colistin|polymyxin',regex=True)})
    return x.groupby('assembly_ID',as_index=False).agg(mcr_detected=('mcr','max'),colistin_keyword_detected=('col','max'))
def kleb_flags(df):
    if df.empty:return pd.DataFrame(columns=['assembly_ID','kleborate_known_colistin_signal','kleborate_signal_text'])
    ac=assembly_col(df);cols=[c for c in df.columns if re.search('colistin|polymyxin|mcr|mgrb|pmrb|phoq|crrb',c,re.I)]
    if not cols:
        cols=[c for c in df.columns if re.search('resistance|amr',c,re.I)]
    x=pd.DataFrame({'assembly_ID':df[ac].map(norm_acc) if ac else 'UNKNOWN'})
    x['kleborate_signal_text']=df[cols].fillna('').astype(str).agg(' | '.join,axis=1) if cols else ''
    x['kleborate_known_colistin_signal']=x.kleborate_signal_text.str.lower().str.contains(r'mcr|mgrb|pmrb|phoq|crrb|colistin|polymyxin',regex=True) & ~x.kleborate_signal_text.str.lower().isin(['','-','none','nan'])
    return x.groupby('assembly_ID',as_index=False).agg(kleborate_known_colistin_signal=('kleborate_known_colistin_signal','max'),kleborate_signal_text=('kleborate_signal_text',lambda s:' || '.join(sorted(set(v for v in s if v and v not in ['-','nan'])))))
def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);cls=read(a.classification);support=read(a.support);seqdir=Path(a.sequence_dir);vr=read(seqdir/'TARGETED_VARIANT_DIRECT_RECALL.csv');ur=read(seqdir/'UNITIG_DIRECT_CONTEXT_RECALL.csv');amr=flag_mcr(read(a.amrfinder,sep='\t'));kleb=kleb_flags(read(a.kleborate,sep='\t'));lit=read(a.literature);ctx=read(a.context_amrfinder,sep='\t')
    if cls.empty:
        summary={'status':'NO_CANDIDATE_TO_MECHANISM_AUDIT','n_advance_to_independent_validation':0};(out/'MECHANISM_CONTEXT_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');pd.DataFrame().to_csv(out/'CANDIDATE_MECHANISM_CONTEXT_VERDICT.csv',index=False);print(json.dumps(summary));return
    rows=[]
    for _,c in cls.iterrows():
        cand=str(c.candidate);kind=str(c.candidate_class);carriers=support[support.candidate.astype(str).eq(cand)].copy() if len(support) else pd.DataFrame();ncar=len(carriers)
        if len(carriers):
            carriers['assembly_ID']=carriers.assembly_ID.map(norm_acc);z=carriers[['assembly_ID']].drop_duplicates().merge(amr,on='assembly_ID',how='left').merge(kleb,on='assembly_ID',how='left');mcr_fraction=float(z.mcr_detected.fillna(False).astype(bool).mean()) if len(z) else np.nan;kleb_fraction=float(z.kleborate_known_colistin_signal.fillna(False).astype(bool).mean()) if len(z) else np.nan
        else:mcr_fraction=kleb_fraction=np.nan
        sequence_tested=sequence_confirmed=sequence_edge=0;sequence_gate=False;sequence_note=''
        if kind=='targeted_variant':
            q=vr[vr.candidate.astype(str).eq(cand)] if len(vr) else pd.DataFrame();sequence_tested=len(q);sequence_confirmed=int(q.get('mutation_confirmed',pd.Series(dtype=str)).astype(str).str.lower().eq('true').sum()) if len(q) else 0;rate=sequence_confirmed/sequence_tested if sequence_tested else 0;sequence_gate=sequence_tested>=5 and rate>=.90;sequence_note=f'direct_recall={sequence_confirmed}/{sequence_tested}'
        elif kind=='unitig':
            q=ur[ur.candidate.astype(str).eq(cand)] if len(ur) else pd.DataFrame();sequence_tested=len(q);sequence_confirmed=int(q.get('perfect_full_length',pd.Series(dtype=str)).astype(str).str.lower().eq('true').sum()) if len(q) else 0;edge=pd.to_numeric(q.get('distance_to_nearest_contig_edge',pd.Series(dtype=float)),errors='coerce');sequence_edge=int((edge<500).sum());rate=sequence_confirmed/sequence_tested if sequence_tested else 0;edge_fraction=sequence_edge/sequence_tested if sequence_tested else 1;sequence_gate=sequence_tested>=5 and rate>=.90 and edge_fraction<=.20;sequence_note=f'perfect={sequence_confirmed}/{sequence_tested};edge_lt500={sequence_edge}'
        elif kind=='gene_burden':
            sequence_note='aggregate burden; component-level functional/context audit required';sequence_gate=False
        litrow=lit[lit.candidate.astype(str).eq(cand)].iloc[0] if len(lit) and 'candidate' in lit.columns and any(lit.candidate.astype(str).eq(cand)) else None;exact_hits=int(float(litrow.get('exact_query_hit_count',0))) if litrow is not None and pd.notna(litrow.get('exact_query_hit_count')) else 0
        known_pathway=str(c.get('known_colistin_pathway','')).lower()=='true';strong_known_conf=(np.isfinite(mcr_fraction) and mcr_fraction>.20) or (np.isfinite(kleb_fraction) and kleb_fraction>.50)
        literature_clear=exact_hits==0
        advance=bool(sequence_gate and not strong_known_conf and literature_clear)
        if known_pathway and advance:level='KNOWN_PATHWAY_EXACT_VARIANT_CANDIDATE_REQUIRES_MANUAL_NOVELTY_AND_FUNCTIONAL_TEST'
        elif advance:level='COMPUTATIONAL_ASSOCIATION_CANDIDATE_ADVANCES_TO_INDEPENDENT_DATA'
        elif exact_hits>0:level='AUTOMATED_LITERATURE_HIT_REQUIRES_MANUAL_REVIEW_BEFORE_NOVELTY'
        elif strong_known_conf:level='LIKELY_CONFOUNDED_BY_KNOWN_COLISTIN_SIGNAL'
        elif not sequence_gate:level='FAILED_OR_INCOMPLETE_DIRECT_SEQUENCE_CONTEXT_GATE'
        else:level='DOES_NOT_ADVANCE'
        rows.append({'candidate':cand,'candidate_class':kind,'gene':c.get('gene'),'known_colistin_pathway':known_pathway,'n_carrier_assemblies':ncar,'mcr_carrier_fraction':mcr_fraction,'kleborate_known_signal_fraction':kleb_fraction,'sequence_tested':sequence_tested,'sequence_confirmed':sequence_confirmed,'sequence_near_contig_edge':sequence_edge,'sequence_gate':sequence_gate,'sequence_note':sequence_note,'automated_exact_literature_hits':exact_hits,'automated_literature_clear':literature_clear,'known_mechanism_confounding':strong_known_conf,'advance_to_independent_validation':advance,'audit_classification':level,'claim_boundary':'No causality or novelty claim; automated literature zero-hit is not proof of novelty.'})
    verdict=pd.DataFrame(rows).sort_values(['advance_to_independent_validation','sequence_gate','candidate'],ascending=[False,False,True]);verdict.to_csv(out/'CANDIDATE_MECHANISM_CONTEXT_VERDICT.csv',index=False);advance=verdict[verdict.advance_to_independent_validation].copy();advance.to_csv(out/'CANDIDATES_ADVANCING_TO_INDEPENDENT_VALIDATION.csv',index=False)
    summary={'status':'CANDIDATES_ADVANCE_TO_INDEPENDENT_VALIDATION' if len(advance) else 'NO_CANDIDATE_SURVIVED_MECHANISM_CONTEXT_GATE','n_candidates':len(verdict),'n_sequence_gate':int(verdict.sequence_gate.sum()),'n_known_mechanism_confounding':int(verdict.known_mechanism_confounding.sum()),'n_automated_exact_literature_hits':int((verdict.automated_exact_literature_hits>0).sum()),'n_advance_to_independent_validation':len(advance),'advancing_candidates':advance.candidate.astype(str).tolist(),'claim_boundary':'Advancing candidates remain non-causal genomic associations. Independent phenotype-linked genomes and laboratory validation are required.'}
    (out/'MECHANISM_CONTEXT_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n');report=['# Colistin candidate mechanism and context gate','',f"- Candidates audited: **{len(verdict):,}**",f"- Direct sequence/context gate: **{summary['n_sequence_gate']:,}**",f"- Known-mechanism confounding: **{summary['n_known_mechanism_confounding']:,}**",f"- Advance to independent data: **{len(advance):,}**",'',summary['claim_boundary'],'','## Candidate verdict','',verdict.to_markdown(index=False)];(out/'MECHANISM_CONTEXT_AUDIT_REPORT.md').write_text('\n'.join(report)+'\n');hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'];(out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n');print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
