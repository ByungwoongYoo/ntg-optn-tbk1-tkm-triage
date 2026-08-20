#!/usr/bin/env python3
"""Extract exact genome-level read abundance from a CAMI reads_mapping member.

This is an evaluation-only operation and must run after the truth-blind assembly freeze.
It stops after the mapping member and does not materialize the multi-gigabyte reads archive.
"""
from __future__ import annotations
import argparse,csv,gzip,io,json,subprocess,tarfile,time
from collections import defaultdict
from pathlib import Path

def args():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--url');g.add_argument('--archive');p.add_argument('--sample-id',required=True);p.add_argument('--paired',action='store_true');p.add_argument('--out',required=True);return p.parse_args()

def source(a):
    if a.archive:return open(a.archive,'rb'),None
    proc=subprocess.Popen(['curl','-fsSL','--retry','8','--retry-all-errors',a.url],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if proc.stdout is None:raise RuntimeError('curl stdout unavailable')
    return proc.stdout,proc

def main():
    a=args();out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);raw,proc=source(a);counts=defaultdict(int);tax={};total=0;member_name=None;t0=time.time()
    try:
        with tarfile.open(fileobj=raw,mode='r|gz') as tf:
            for m in tf:
                if not m.isfile() or not m.name.endswith('reads_mapping.tsv.gz'):continue
                member_name=m.name;ex=tf.extractfile(m)
                if ex is None:raise RuntimeError('mapping member could not be opened')
                with gzip.GzipFile(fileobj=ex) as z,io.TextIOWrapper(z,encoding='utf-8',errors='replace') as f:
                    reader=csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t',fieldnames=['anonymous_read_id','genome_id','tax_id','read_id'])
                    for row in reader:
                        gid=(row.get('genome_id') or '').strip()
                        if not gid or gid=='genome_id':continue
                        counts[gid]+=1;tax[gid]=(row.get('tax_id') or '').strip();total+=1
                break
    finally:
        try:raw.close()
        except Exception:pass
        if proc is not None:
            proc.terminate()
            try:proc.wait(timeout=15)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
    if member_name is None or total==0:raise RuntimeError('No populated reads_mapping.tsv.gz member found')
    divisor=2 if a.paired else 1;denom=total/divisor;rows=[]
    for gid,n in sorted(counts.items()):
        molecules=n/divisor;rows.append({'sample_id':a.sample_id,'genome_id':gid,'tax_id':tax.get(gid,''),'read_records':n,'estimated_read_molecules':molecules,'relative_read_fraction':molecules/denom if denom else 0,'relative_read_percent':100*molecules/denom if denom else 0})
    with open(out,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    summary={'sample_id':a.sample_id,'mapping_member':member_name,'paired':a.paired,'total_read_records':total,'estimated_read_molecules':denom,'n_genomes':len(rows),'elapsed_seconds':time.time()-t0,'evaluation_only':True,'boundary':'Read-to-genome mapping and abundance were accessed only after truth-blind assembly outputs and hashes were frozen.'}
    out.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
