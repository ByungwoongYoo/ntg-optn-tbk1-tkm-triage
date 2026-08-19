#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
from fasta_utils import read_fasta,write_record

def args():
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--metadata',required=True);p.add_argument('--evidence',required=True);p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--ablation',choices=['full','no_longitudinal','no_long','no_consensus'],default='full');return p.parse_args()

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);cfg=json.load(open(a.config));seqs=dict(read_fasta(a.fasta));meta={}
    with open(a.metadata,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):meta[r['representative_id']]=r
    ev=defaultdict(dict); samples=set()
    with open(a.evidence,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            samples.add(r['sample_id']);ev[r['contig_id']][r['sample_id']]={k:float(v or 0) for k,v in r.items() if k not in ('contig_id','sample_id')}
    w=cfg['weights'].copy()
    if a.ablation=='no_longitudinal': w['temporal_recurrence']=0
    if a.ablation=='no_long': w['long_breadth']=w['long_spanning']=0
    if a.ablation=='no_consensus': w['source_consensus']=0
    sw=sum(w.values()); w={k:v/sw for k,v in w.items()} if sw else w
    rows=[]; selected=[]
    for cid,seq in seqs.items():
        m=meta[cid]; by=ev.get(cid,{})
        supports=[]; sb=[];lb=[];sp=0
        for s in sorted(samples):
            e=by.get(s,{})
            sb.append(e.get('short_breadth',0)); lb.append(e.get('long_breadth',0)); sp+=int(e.get('long_spanning_reads',0))
            breadth=max(e.get('short_breadth',0),0 if a.ablation=='no_long' else e.get('long_breadth',0))
            depth=max(e.get('short_depth',0),0 if a.ablation=='no_long' else e.get('long_depth',0))
            supports.append(breadth>=cfg['support_breadth'] and depth>=cfg['support_depth'])
        t=sum(supports); nsrc=int(m['source_count'])
        feats={'length':min(1,math.log1p(len(seq))/math.log1p(100000)),'source_consensus':min(1,nsrc/3),'temporal_recurrence':t/max(1,len(samples)),'short_breadth':max(sb or [0]),'long_breadth':0 if a.ablation=='no_long' else max(lb or [0]),'long_spanning':0 if a.ablation=='no_long' else min(1,math.log1p(sp)/math.log1p(6))}
        score=sum(w[k]*feats[k] for k in w)
        scopes=set(filter(None,m.get('scopes','').split(','))); single='single' in scopes
        consensus=(nsrc>=cfg['minimum_assembler_sources_for_consensus']) and a.ablation!='no_consensus'
        temporal=(t>=cfg['minimum_rescue_timepoints']) and a.ablation!='no_longitudinal'
        longok=(sp>=cfg['minimum_unique_long_spans']) and a.ablation!='no_long'
        usable_breadths=sb if a.ablation=='no_long' else sb+lb
        basic=len(seq)>=cfg['minimum_contig_length'] and seq.count('N')/len(seq)<=cfg['maximum_n_fraction'] and max(usable_breadths+[0])>=cfg['support_breadth']
        gate=basic and score>=cfg['selection_score_minimum'] and (consensus or temporal or longok or (single and t>=1))
        rescue=gate and (not single) and (temporal or longok) and not consensus
        row={'contig_id':cid,'selected':str(gate).lower(),'rescue':str(rescue).lower(),'score':score,'length':len(seq),'source_count':nsrc,'timepoints_supported':t,'n_timepoints':len(samples),'max_short_breadth':max(sb or [0]),'max_long_breadth':max(lb or [0]),'long_spanning_reads':sp,'scopes':m.get('scopes',''),'assemblers':m.get('assemblers',''),'ablation':a.ablation,**{f'feature_{k}':v for k,v in feats.items()}}
        rows.append(row)
        if gate:selected.append(cid)
    rows.sort(key=lambda r:(r['selected']!='true',-r['score'],-r['length'],r['contig_id']))
    if not rows: raise SystemExit('No representative contigs')
    with open(out/'LANTERN_SELECTION.tsv','w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
    with open(out/'LANTERN_ASSEMBLY.fasta','w') as f:
        for cid in selected:write_record(f,cid,seqs[cid])
    summary={'version':cfg['version'],'ablation':a.ablation,'n_input':len(seqs),'n_selected':len(selected),'n_rescue':sum(r['rescue']=='true' for r in rows),'total_selected_bp':sum(len(seqs[x]) for x in selected),'truth_blind':cfg.get('truth_blind',False),'config':cfg}
    (out/'LANTERN_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
