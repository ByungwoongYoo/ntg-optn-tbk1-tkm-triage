#!/usr/bin/env python3
"""Stream only a prespecified prefix of a CAMI FASTQ member and terminate transfer early."""
from __future__ import annotations
import argparse,gzip,io,json,subprocess,tarfile,time
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--url');g.add_argument('--archive')
    p.add_argument('--sample-id',required=True);p.add_argument('--mode',choices=['short','long'],required=True);p.add_argument('--molecules',type=int,required=True);p.add_argument('--out-dir',required=True);return p.parse_args()

def source(a):
    if a.archive:return open(a.archive,'rb'),None
    proc=subprocess.Popen(['curl','-fsSL','--retry','8','--retry-all-errors',a.url],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if proc.stdout is None:raise RuntimeError('curl stdout unavailable')
    return proc.stdout,proc

def records(stream):
    while True:
        h=stream.readline()
        if not h:return
        s=stream.readline();p=stream.readline();q=stream.readline()
        if not(s and p and q):raise RuntimeError('truncated FASTQ')
        if not h.startswith('@') or not p.startswith('+'):raise RuntimeError('invalid FASTQ')
        yield h,s,p,q

def main():
    a=parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True);raw,proc=source(a);member=None;kept=0;t0=time.time();paths=[]
    try:
        with tarfile.open(fileobj=raw,mode='r|gz') as tf:
            for m in tf:
                if not m.isfile() or not (m.name.endswith('.fq.gz') or m.name.endswith('.fastq.gz')):continue
                member=m.name;ex=tf.extractfile(m)
                if ex is None:raise RuntimeError('FASTQ member could not be opened')
                with gzip.GzipFile(fileobj=ex) as z,io.TextIOWrapper(z,encoding='utf-8',errors='strict',newline='') as text:
                    if a.mode=='short':
                        p1=out/f'{a.sample_id}_R1.fastq.gz';p2=out/f'{a.sample_id}_R2.fastq.gz';paths=[p1,p2]
                        with gzip.open(p1,'wt',newline='') as o1,gzip.open(p2,'wt',newline='') as o2:
                            it=iter(records(text))
                            for _ in range(a.molecules):
                                try:r1=next(it);r2=next(it)
                                except StopIteration:break
                                id1=r1[0].split()[0].removesuffix('/1');id2=r2[0].split()[0].removesuffix('/2')
                                if id1!=id2:raise RuntimeError(f'pair mismatch {id1} {id2}')
                                o1.writelines(r1);o2.writelines(r2);kept+=1
                    else:
                        p=out/f'{a.sample_id}_long.fastq.gz';paths=[p]
                        with gzip.open(p,'wt',newline='') as o:
                            for i,r in enumerate(records(text)):
                                if i>=a.molecules:break
                                o.writelines(r);kept+=1
                break
    finally:
        try:raw.close()
        except Exception:pass
        if proc is not None:
            proc.terminate()
            try:proc.wait(timeout=15)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
    if member is None or kept==0:raise RuntimeError('No FASTQ records were retained')
    report={'sample_id':a.sample_id,'mode':a.mode,'member':member,'requested_molecules':a.molecules,'kept_molecules':kept,'outputs':[str(x) for x in paths],'elapsed_seconds':time.time()-t0,'selection':'prespecified_archive_prefix','boundary':'Prefix size and representativeness thresholds must be frozen on development data before held-out use.'}
    (out/f'{a.sample_id}_{a.mode}_prefix.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
