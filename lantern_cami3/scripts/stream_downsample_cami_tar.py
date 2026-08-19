#!/usr/bin/env python3
"""Stream a CAMI tar.gz read archive and deterministically downsample FASTQ.

CAMI III Toy short-read archives contain one interleaved FASTQ.GZ member; long-read
archives contain one FASTQ.GZ member. This script also supports archives with separate
R1/R2 members. The outer tar is never materialized.
"""
from __future__ import annotations
import argparse,gzip,hashlib,io,json,re,subprocess,sys,tarfile
from dataclasses import dataclass,asdict
from pathlib import Path
from typing import BinaryIO,Iterator,TextIO

FASTQ_EXT=re.compile(r"\.(?:fastq|fq)(?:\.gz)?$",re.I)
M1=[re.compile(p,re.I) for p in [r"(?:^|[._-])R?1(?:[._-]|$)",r"(?:^|[._-])read1(?:[._-]|$)"]]
M2=[re.compile(p,re.I) for p in [r"(?:^|[._-])R?2(?:[._-]|$)",r"(?:^|[._-])read2(?:[._-]|$)"]]

@dataclass
class MemberStats:
    member:str; role:str; total_records:int=0; kept_records:int=0; selected_hash_xor:int=0; selected_hash_sum:int=0; output:str=""

def parse_args():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--url');g.add_argument('--archive')
    p.add_argument('--sample-id',required=True);p.add_argument('--mode',choices=['short','long'],required=True);p.add_argument('--fraction',type=float,required=True);p.add_argument('--seed',type=int,default=20260819);p.add_argument('--out-dir',required=True);p.add_argument('--curl-retries',type=int,default=8);p.add_argument('--min-kept',type=int,default=1000);return p.parse_args()

def norm_id(h:str)->str:
    t=h.strip().split()[0].lstrip('@');t=re.sub(r'/(?:1|2)$','',t);t=re.sub(r'(?:[._-])R?[12]$','',t,flags=re.I);return t

def hash64(s:str,seed:int)->int:
    return int.from_bytes(hashlib.blake2b(s.encode('utf-8','replace'),digest_size=8,key=seed.to_bytes(8,'little',signed=False)).digest(),'big')

def keep(h:int,f:float)->bool:return f>=1 or h<int(f*(1<<64))

def role(name:str)->str:
    b=Path(name).name
    if any(p.search(b) for p in M1):return 'R1'
    if any(p.search(b) for p in M2):return 'R2'
    return 'unknown'

def open_member(raw:BinaryIO,name:str)->TextIO:
    z=gzip.GzipFile(fileobj=raw,mode='rb') if name.lower().endswith('.gz') else raw
    return io.TextIOWrapper(z,encoding='utf-8',errors='replace',newline='')

def records(h:TextIO)->Iterator[tuple[str,str,str,str]]:
    while True:
        a=h.readline()
        if not a:return
        b=h.readline();c=h.readline();d=h.readline()
        if not (b and c and d):raise ValueError('Truncated FASTQ')
        if not a.startswith('@') or not c.startswith('+'):raise ValueError(f'Bad FASTQ header {a[:80]!r}')
        yield a,b,c,d

def source(args):
    if args.archive:return open(args.archive,'rb'),None
    cmd=['curl','-fL','--retry',str(args.curl_retries),'--retry-all-errors','--connect-timeout','30','--speed-time','120','--speed-limit','1024',args.url]
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=sys.stderr);assert p.stdout is not None;return p.stdout,p

def update(st:MemberStats,hv:int):
    st.kept_records+=1;st.selected_hash_xor^=hv;st.selected_hash_sum=(st.selected_hash_sum+hv)&((1<<64)-1)

def process_long(raw:BinaryIO,name:str,args,out:Path)->MemberStats:
    op=out/f'{args.sample_id}_long.fastq.gz';st=MemberStats(name,'long',output=str(op))
    with open_member(raw,name) as inp,gzip.open(op,'wt',encoding='utf-8',newline='') as dst:
        for rec in records(inp):
            st.total_records+=1;hv=hash64(norm_id(rec[0]),args.seed)
            if keep(hv,args.fraction):dst.writelines(rec);update(st,hv)
    return st

def process_interleaved(raw:BinaryIO,name:str,args,out:Path)->list[MemberStats]:
    o1=out/f'{args.sample_id}_R1.fastq.gz';o2=out/f'{args.sample_id}_R2.fastq.gz'
    s1=MemberStats(name,'R1',output=str(o1));s2=MemberStats(name,'R2',output=str(o2))
    with open_member(raw,name) as inp,gzip.open(o1,'wt',encoding='utf-8',newline='') as a,gzip.open(o2,'wt',encoding='utf-8',newline='') as b:
        it=records(inp)
        while True:
            try:r1=next(it)
            except StopIteration:break
            try:r2=next(it)
            except StopIteration:raise ValueError('Odd number of records in interleaved FASTQ')
            s1.total_records+=1;s2.total_records+=1
            i1=norm_id(r1[0]);i2=norm_id(r2[0])
            if i1!=i2:raise ValueError(f'Interleaved pair IDs differ: {r1[0].strip()} vs {r2[0].strip()}')
            hv=hash64(i1,args.seed)
            if keep(hv,args.fraction):a.writelines(r1);b.writelines(r2);update(s1,hv);update(s2,hv)
    return [s1,s2]

def process_separate(raw:BinaryIO,name:str,r:str,args,out:Path)->MemberStats:
    op=out/f'{args.sample_id}_{r}.fastq.gz';st=MemberStats(name,r,output=str(op))
    with open_member(raw,name) as inp,gzip.open(op,'wt',encoding='utf-8',newline='') as dst:
        for rec in records(inp):
            st.total_records+=1;hv=hash64(norm_id(rec[0]),args.seed)
            if keep(hv,args.fraction):dst.writelines(rec);update(st,hv)
    return st

def main():
    a=parse_args()
    if not 0<a.fraction<=1:raise SystemExit('fraction must be (0,1]')
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    for p in out.glob(f'{a.sample_id}_*.fastq.gz'):p.unlink()
    fh,proc=source(a);stats=[];fastq_members=[]
    try:
        with tarfile.open(fileobj=fh,mode='r|gz') as tf:
            for m in tf:
                if not m.isfile() or not FASTQ_EXT.search(m.name):continue
                raw=tf.extractfile(m)
                if raw is None:raise RuntimeError(m.name)
                fastq_members.append(m.name)
                if a.mode=='long':stats.append(process_long(raw,m.name,a,out))
                else:
                    r=role(m.name)
                    if r=='unknown':stats.extend(process_interleaved(raw,m.name,a,out))
                    else:stats.append(process_separate(raw,m.name,r,a,out))
    finally:fh.close()
    if proc is not None:
        rc=proc.wait()
        if rc:raise SystemExit(f'curl exit {rc}')
    if a.mode=='short':
        by={s.role:s for s in stats}
        if set(by)!={'R1','R2'}:raise SystemExit(f'Expected R1/R2, got {sorted(by)} from {fastq_members}')
        x,y=by['R1'],by['R2']
        if (x.total_records,x.kept_records,x.selected_hash_xor,x.selected_hash_sum)!=(y.total_records,y.kept_records,y.selected_hash_xor,y.selected_hash_sum):raise SystemExit('pair audit failed')
        if x.kept_records<a.min_kept:raise SystemExit(f'too few pairs {x.kept_records}')
    elif sum(s.kept_records for s in stats)<a.min_kept:raise SystemExit('too few long reads')
    summary={'sample_id':a.sample_id,'mode':a.mode,'source':a.url or a.archive,'fraction':a.fraction,'seed':a.seed,'members':[asdict(s) for s in stats],'total_records':sum(s.total_records for s in stats),'kept_records':sum(s.kept_records for s in stats),'pair_integrity':True}
    (out/f'{a.sample_id}_{a.mode}_downsample.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
