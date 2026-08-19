#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
from fasta_utils import read_fasta,write_record

def parse_args():
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--metadata',required=True);p.add_argument('--evidence',required=True);p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--ablation',choices=['full','no_longitudinal','no_long','no_consensus'],default='full');return p.parse_args()

def clamp(x):return max(0.0,min(1.0,x))

def main():
    a=parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);cfg=json.load(open(a.config));seqs=dict(read_fasta(a.fasta));meta={r['representative_id']:r for r in csv.DictReader(open(a.metadata),delimiter='\t')}
    if set(seqs)!=set(meta):raise ValueError('FASTA/metadata mismatch')
    ev=defaultdict(dict);samples=set()
    for r in csv.DictReader(open(a.evidence),delimiter='\t'):
        samples.add(r['sample_id']);ev[r['contig_id']][r['sample_id']]={k:float(v or 0) for k,v in r.items() if k not in ('contig_id','sample_id')}
    if not samples:raise ValueError('No samples in evidence')
    rows=[];selected=[]
    for cid,seq in seqs.items():
        m=meta[cid];strong=weak=0;support_values=[];short_b=[];long_b=[];spans=0
        for sid in sorted(samples):
            e=ev.get(cid,{}).get(sid,{});sr=e.get('short_reads',0);sb=e.get('short_breadth',0);sd=e.get('short_depth',0);lb=e.get('long_breadth',0);ld=e.get('long_depth',0);sp=int(e.get('long_spanning_reads',0));spans+=sp
            if a.ablation=='no_long':lb=ld=0;sp=0
            short_b.append(sb);long_b.append(lb);short_signal=sb*clamp(sd/cfg['depth_saturation_short']);long_signal=lb*clamp(ld/cfg['depth_saturation_long']);support_values.append(max(short_signal,long_signal))
            is_strong=(sb>=cfg['strong_breadth'] and sd>=cfg['strong_depth']) or (lb>=cfg['strong_breadth'] and ld>=cfg['strong_depth'])
            is_weak=(sb>=cfg['weak_breadth'] and sr>=cfg['weak_short_reads']) or lb>=cfg['weak_long_breadth'] or sp>=1
            strong+=int(is_strong);weak+=int(is_weak)
        n=len(samples);temporal_score=sum(support_values)/n;temporal_gate=strong>=cfg['minimum_strong_timepoints'] and weak>=cfg['minimum_weak_timepoints']
        if a.ablation=='no_longitudinal':temporal_score=0;temporal_gate=False
        consensus_count=int(m.get('support_assembler_count') or m.get('assembler_count') or 0);consensus=consensus_count>=cfg['minimum_assembler_sources_for_consensus'] and a.ablation!='no_consensus';long_gate=spans>=cfg['minimum_unique_long_spans'] and a.ablation!='no_long';scopes=set(filter(None,m.get('scopes','').split(',')));single='single' in scopes
        feats={'length':min(1,math.log1p(len(seq))/math.log1p(100000)),'source_consensus':min(1,consensus_count/3),'temporal_recurrence':temporal_score,'short_breadth':max(short_b or [0]),'long_breadth':0 if a.ablation=='no_long' else max(long_b or [0]),'long_spanning':0 if a.ablation=='no_long' else min(1,math.log1p(spans)/math.log1p(6))}
        weights=dict(cfg['weights'])
        if a.ablation=='no_longitudinal':weights['temporal_recurrence']=0
        if a.ablation=='no_long':weights['long_breadth']=weights['long_spanning']=0
        if a.ablation=='no_consensus':weights['source_consensus']=0
        z=sum(weights.values());weights={k:v/z for k,v in weights.items()} if z else weights;score=sum(weights[k]*feats[k] for k in weights);nfrac=seq.count('N')/len(seq);basic=len(seq)>=cfg['minimum_contig_length'] and nfrac<=cfg['maximum_n_fraction'] and strong>=1;evidence_gate=consensus or temporal_gate or long_gate or (single and strong>=1);gate=basic and score>=cfg['selection_score_minimum'] and evidence_gate;rescue=gate and not single and not consensus and (temporal_gate or long_gate)
        rows.append({'contig_id':cid,'selected':str(gate).lower(),'rescue':str(rescue).lower(),'score':score,'length':len(seq),'n_fraction':nfrac,'support_assembler_count':consensus_count,'strong_timepoints':strong,'weak_timepoints':weak,'n_timepoints':n,'temporal_gate':str(temporal_gate).lower(),'consensus_gate':str(consensus).lower(),'long_gate':str(long_gate).lower(),'max_short_breadth':max(short_b or [0]),'max_long_breadth':max(long_b or [0]),'long_spanning_reads':spans,'scopes':m.get('scopes',''),'assemblers':m.get('support_assemblers') or m.get('assemblers',''),'ablation':a.ablation,**{f'feature_{k}':v for k,v in feats.items()}});selected.extend([cid] if gate else [])
    rows.sort(key=lambda r:(r['selected']!='true',-r['score'],-r['length'],r['contig_id']))
    with open(out/'LANTERN_SELECTION.tsv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    with open(out/'LANTERN_ASSEMBLY.fasta','w') as f:
        for cid in selected:write_record(f,cid,seqs[cid])
    summary={'version':cfg['version'],'selector':'lantern_select_v2','ablation':a.ablation,'n_input':len(seqs),'n_selected':len(selected),'n_rescue':sum(r['rescue']=='true' for r in rows),'total_selected_bp':sum(len(seqs[x]) for x in selected),'truth_blind':True,'config':cfg};(out/'LANTERN_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
