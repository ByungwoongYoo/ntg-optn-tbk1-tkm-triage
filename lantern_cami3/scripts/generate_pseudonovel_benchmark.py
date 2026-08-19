#!/usr/bin/env python3
"""Generate a deterministic, reference-withheld longitudinal assembly stress test.

The generated target genomes are mutated descendants of parental sequences. Only the
parents are written to the supplied reference collection, while the target descendants
used for read simulation are withheld. The assembly pipeline is de novo and does not use
this reference collection; the benchmark therefore tests sequence-level recovery of unseen,
strain-related genomes rather than reference-assisted taxonomic assignment.

This is a synthetic stress test, not evidence that a biological taxon is novel.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
from pathlib import Path

DNA = "ACGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--length", type=int, default=160000)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--insert-mean", type=int, default=350)
    p.add_argument("--error-rate", type=float, default=0.005)
    return p.parse_args()


def gc_sequence(length: int, gc: float, rng: random.Random) -> str:
    weights = [("A", (1-gc)/2), ("T", (1-gc)/2), ("G", gc/2), ("C", gc/2)]
    cumulative=[]; total=0.0
    for base,w in weights:
        total += w; cumulative.append((total,base))
    chars=[]
    for _ in range(length):
        x=rng.random()
        chars.append(next(base for threshold,base in cumulative if x<=threshold))
    return "".join(chars)


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    out=list(seq)
    for i,b in enumerate(out):
        if rng.random() < rate:
            out[i]=rng.choice([x for x in DNA if x!=b])
    return "".join(out)


def add_repeat_structure(seq: str, repeat: str, offsets: list[int]) -> str:
    out=list(seq)
    for offset in offsets:
        end=min(len(out),offset+len(repeat)); out[offset:end]=repeat[:end-offset]
    return "".join(out)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def introduce_errors(seq: str, rate: float, rng: random.Random) -> str:
    out=list(seq)
    for i,b in enumerate(out):
        if rng.random()<rate:
            out[i]=rng.choice([x for x in DNA if x!=b])
    return "".join(out)


def write_fasta(path: Path, records: list[tuple[str,str]]) -> None:
    with path.open("w") as fh:
        for name,seq in records:
            fh.write(f">{name}\n")
            for i in range(0,len(seq),80): fh.write(seq[i:i+80]+"\n")


def simulate_pairs(genome_id: str, seq: str, coverage: float, read_length: int, insert_mean: int, error_rate: float, rng: random.Random) -> list[tuple[str,str,str]]:
    n=max(0,round(coverage*len(seq)/(2*read_length)))
    pairs=[]
    for index in range(n):
        insert=max(2*read_length,int(rng.gauss(insert_mean,25)))
        start=rng.randrange(len(seq))
        fragment=(seq+seq)[start:start+insert]
        if len(fragment)<2*read_length: continue
        r1=introduce_errors(fragment[:read_length],error_rate,rng)
        r2=introduce_errors(reverse_complement(fragment[-read_length:]),error_rate,rng)
        pairs.append((f"{genome_id}_{index}",r1,r2))
    return pairs


def write_fastq_pair(prefix: Path, records: list[tuple[str,str,str]], rng: random.Random) -> tuple[int,str,str]:
    rng.shuffle(records)
    r1=prefix.with_name(prefix.name+"_R1.fastq.gz"); r2=prefix.with_name(prefix.name+"_R2.fastq.gz")
    h1=hashlib.sha256(); h2=hashlib.sha256()
    with gzip.open(r1,"wt") as a, gzip.open(r2,"wt") as b:
        for name,s1,s2 in records:
            x=f"@{name}/1\n{s1}\n+\n{'I'*len(s1)}\n"; y=f"@{name}/2\n{s2}\n+\n{'I'*len(s2)}\n"
            a.write(x); b.write(y); h1.update(x.encode()); h2.update(y.encode())
    return len(records),h1.hexdigest(),h2.hexdigest()


def main() -> None:
    args=parse_args(); out=args.out; out.mkdir(parents=True,exist_ok=True)
    rng=random.Random(args.seed)
    gc_values=[0.34,0.43,0.52,0.64]
    divergence=[0.02,0.05,0.08,0.12]
    coverage_profiles=[(5.0,1.2),(1.2,5.0),(3.0,3.0),(1.5,1.5)]
    parents=[]; targets=[]; manifest=[]
    shared_repeat=gc_sequence(2200,0.50,random.Random(args.seed+9999))
    for i,(gc,div,cov) in enumerate(zip(gc_values,divergence,coverage_profiles),1):
        parent=gc_sequence(args.length,gc,random.Random(args.seed+100*i))
        parent=add_repeat_structure(parent,shared_repeat,[20000,70000,125000])
        target=mutate(parent,div,random.Random(args.seed+1000*i))
        target=add_repeat_structure(target,shared_repeat,[25000,76000,130000])
        pid=f"PARENT_{i}"; gid=f"PSEUDONOVEL_{i}"
        parents.append((pid,parent)); targets.append((gid,target))
        manifest.append({'genome_id':gid,'withheld_parent':pid,'designed_divergence':div,'designed_ani':1-div,'gc':gc,'length':len(target),'coverage_t0':cov[0],'coverage_t1':cov[1],'cumulative_coverage':sum(cov),'low_abundance':sum(cov)<=3.1})
    write_fasta(out/'reference_parents_only.fasta',parents)
    write_fasta(out/'withheld_target_truth.fasta',targets)
    with (out/'truth_mapping.tsv').open('w',newline='') as fh:
        fields=['sequence_id','genome_id','original_id','length','sha256','samples']; w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader()
        for gid,seq in targets: w.writerow({'sequence_id':gid,'genome_id':gid,'original_id':gid,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'samples':'T0,T1'})
    with (out/'PSEUDONOVEL_MANIFEST.csv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
    timepoint_records={'T0':[],'T1':[]}
    for row,(gid,seq) in zip(manifest,targets):
        for tp,key in [('T0','coverage_t0'),('T1','coverage_t1')]:
            timepoint_records[tp].extend(simulate_pairs(gid,seq,float(row[key]),args.read_length,args.insert_mean,args.error_rate,random.Random(args.seed+hash((gid,tp))%1000000)))
    read_summary={}
    for tp in ['T0','T1']:
        n,h1,h2=write_fastq_pair(out/tp,timepoint_records[tp],random.Random(args.seed+(0 if tp=='T0' else 1)))
        read_summary[tp]={'pairs':n,'r1_content_sha256':h1,'r2_content_sha256':h2}
    summary={'seed':args.seed,'n_withheld_targets':len(targets),'target_total_bp':sum(len(s) for _,s in targets),'parents_in_reference':[x[0] for x in parents],'targets_absent_from_reference':[x[0] for x in targets],'read_summary':read_summary,'boundary':'Synthetic reference-withheld sequence recovery stress test; not biological novelty evidence.'}
    (out/'BENCHMARK_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
