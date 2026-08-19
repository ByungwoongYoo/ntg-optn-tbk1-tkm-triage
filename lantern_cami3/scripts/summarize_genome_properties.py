#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import Counter
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--recovery',action='append',required=True,help='METHOD=per_genome_recovery.tsv');p.add_argument('--k',type=int,default=15);p.add_argument('--max-kmers',type=int,default=200000);p.add_argument('--out',required=True);return p.parse_args()

def read_fasta(path):
    name=None;seq=[]
    with open(path) as f:
        for line in f:
            line=line.strip().upper()
            if not line:continue
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line)
    if name is not None:yield name,''.join(seq)

def props(path,k,max_kmers):
    total=gc=n=0;seqs=[];sample=[]
    for _,s in read_fasta(path):
        seqs.append(len(s));total+=len(s);gc+=s.count('G')+s.count('C');n+=s.count('N');step=max(1,(max(0,len(s)-k+1)+max_kmers-1)//max_kmers);sample.extend(s[i:i+k] for i in range(0,max(0,len(s)-k+1),step) if 'N' not in s[i:i+k]);sample=sample[:max_kmers]
    counts=Counter(sample);unique=sum(v==1 for v in counts.values());ratio=unique/len(sample) if sample else None;entropy=None
    if sample:
        total_k=len(sample);entropy=-sum((v/total_k)*math.log2(v/total_k) for v in counts.values())/max(1,k*2)
    return {'truth_bp':total,'n_truth_sequences':len(seqs),'mean_truth_sequence_length':sum(seqs)/len(seqs) if seqs else None,'gc_percent':100*gc/(total-n) if total>n else None,'n_fraction':n/total if total else None,'sampled_kmers':len(sample),'unique_sampled_kmer_fraction':ratio,'normalized_kmer_entropy':entropy}

def main():
    a=parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);manifest=list(csv.DictReader(open(a.manifest),delimiter='\t'));pmap={}
    for r in manifest:
        x=props(r['reference_path'],a.k,a.max_kmers);gc=x['gc_percent'];nseq=x['n_truth_sequences'];uniq=x['unique_sampled_kmer_fraction'];gc_class='low_gc' if gc is not None and gc<35 else ('high_gc' if gc is not None and gc>65 else 'mid_gc');frag='highly_fragmented_truth' if nseq>=50 else ('moderately_fragmented_truth' if nseq>=10 else 'low_fragmentation_truth');repeat='repeat_rich_proxy' if uniq is not None and uniq<.80 else 'non_repeat_rich_proxy';pmap[r['genome_id']]={'genome_id':r['genome_id'],**x,'gc_class':gc_class,'fragmentation_class':frag,'repeat_proxy_class':repeat}
    detail=[]
    for spec in a.recovery:
        method,path=spec.split('=',1)
        for r in csv.DictReader(open(path),delimiter='\t'):
            if r['genome_id'] in pmap:detail.append({'method':method,**pmap[r['genome_id']],'recovery_fraction':float(r['recovery_fraction']),'recovered_50':r['recovered_50'],'recovered_90':r['recovered_90']})
    summary=[];groups=[('gc_class',['low_gc','mid_gc','high_gc']),('fragmentation_class',['low_fragmentation_truth','moderately_fragmented_truth','highly_fragmented_truth']),('repeat_proxy_class',['repeat_rich_proxy','non_repeat_rich_proxy'])]
    for method in sorted({r['method'] for r in detail}):
        for field,strata in groups:
            for stratum in strata:
                rows=[r for r in detail if r['method']==method and r[field]==stratum];vals=[r['recovery_fraction'] for r in rows];summary.append({'method':method,'property':field,'stratum':stratum,'n_genomes':len(rows),'mean_recovery':sum(vals)/len(vals) if vals else None,'median_recovery':statistics.median(vals) if vals else None,'recovered_50':sum(r['recovered_50']=='true' for r in rows),'recovered_90':sum(r['recovered_90']=='true' for r in rows)})
    for name,rows in [('GENOME_PROPERTY_DETAIL.tsv',detail),('GENOME_PROPERTY_RECOVERY_SUMMARY.tsv',summary)]:
        fields=list(rows[0]) if rows else ['method','property','stratum','n_genomes','mean_recovery'];
        with open(out/name,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    report={'n_genomes':len(pmap),'k':a.k,'max_sampled_kmers':a.max_kmers,'boundary':'GC and fragmentation strata are direct truth properties. The repeat-rich label is a sampled k-mer uniqueness proxy, not a validated repeat annotation.'};(out/'GENOME_PROPERTY_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
