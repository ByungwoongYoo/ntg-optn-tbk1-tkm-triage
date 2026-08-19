#!/usr/bin/env python3
"""Development-only audit of prefix sampling against the full CAMI read mapping."""
from __future__ import annotations
import argparse,json,math,subprocess,tarfile,gzip,io
from collections import defaultdict
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--url');g.add_argument('--archive')
    p.add_argument('--head-molecules',type=int,required=True);p.add_argument('--paired',action='store_true');p.add_argument('--out',required=True);return p.parse_args()

def source(a):
    if a.archive:return open(a.archive,'rb'),None
    proc=subprocess.Popen(['curl','-fsSL','--retry','8','--retry-all-errors',a.url],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if proc.stdout is None:raise RuntimeError('curl stdout unavailable')
    return proc.stdout,proc

def jsd(p,q):
    keys=set(p)|set(q);m={k:(p.get(k,0)+q.get(k,0))/2 for k in keys}
    def kl(x,y):return sum(v*math.log2(v/y[k]) for k,v in x.items() if v>0 and y.get(k,0)>0)
    return (kl(p,m)+kl(q,m))/2

def main():
    a=parse_args();raw,proc=source(a);full=defaultdict(int);head=defaultdict(int);total=head_total=0;limit=a.head_molecules*(2 if a.paired else 1);member=None
    try:
        with tarfile.open(fileobj=raw,mode='r|gz') as tf:
            for m in tf:
                if not m.isfile() or not m.name.endswith('reads_mapping.tsv.gz'):continue
                member=m.name;ex=tf.extractfile(m)
                if ex is None:raise RuntimeError('mapping member could not be opened')
                with gzip.GzipFile(fileobj=ex) as z,io.TextIOWrapper(z,encoding='utf-8',errors='replace') as f:
                    for line in f:
                        if line.startswith('#'):continue
                        x=line.rstrip().split('\t')
                        if len(x)<2:continue
                        gid=x[1];full[gid]+=1;total+=1
                        if head_total<limit:head[gid]+=1;head_total+=1
                break
    finally:
        try:raw.close()
        except Exception:pass
        if proc is not None:
            proc.terminate()
            try:proc.wait(timeout=15)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
    if member is None or total==0:raise RuntimeError('No mapping records')
    pf={k:v/total for k,v in full.items()};ph={k:v/head_total for k,v in head.items()};keys=set(pf)|set(ph)
    diffs=sorted(((abs(ph.get(k,0)-pf.get(k,0)),k,ph.get(k,0),pf.get(k,0)) for k in keys),reverse=True)
    report={'mapping_member':member,'paired':a.paired,'full_read_records':total,'head_read_records':head_total,'head_molecules':head_total/(2 if a.paired else 1),'full_taxa':len(full),'head_taxa':len(head),'head_taxon_coverage':len(head)/len(full) if full else 0,'jensen_shannon_divergence_bits':jsd(ph,pf),'total_variation_distance':0.5*sum(abs(ph.get(k,0)-pf.get(k,0)) for k in keys),'maximum_absolute_fraction_difference':diffs[0][0] if diffs else 0,'top_absolute_differences':[{'genome_id':k,'head_fraction':h,'full_fraction':f,'absolute_difference':d} for d,k,h,f in diffs[:20]],'development_only':True,'boundary':'This audit may inform the frozen sampling strategy only on the development individual. It must not be repeated to tune held-out samples.'}
    Path(a.out).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
