#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,io,json,re,tarfile
from collections import defaultdict
from pathlib import Path
from fasta_utils import sanitize_seq,safe_id,write_record,seq_sha

def args():
 p=argparse.ArgumentParser();p.add_argument('--gsa-tar',action='append',required=True,help='SAMPLE=tar.gz');p.add_argument('--out',required=True);return p.parse_args()

def genome_from_member(name):
 b=Path(name).name;m=re.match(r'sample\d+_(.+)_gsa\.fasta\.gz$',b);return m.group(1) if m else re.sub(r'_gsa\.fasta\.gz$','',b)

def parse_fasta_text(text):
 name=None;seq=[]
 for line in text:
  line=line.strip()
  if not line:continue
  if line.startswith('>'):
   if name is not None:yield name,''.join(seq)
   name=line[1:].split()[0];seq=[]
  else:seq.append(line)
 if name is not None:yield name,''.join(seq)

def main():
 a=args();out=Path(a.out);refs=out/'references';refs.mkdir(parents=True,exist_ok=True);genomes=defaultdict(dict);presence=defaultdict(set)
 for spec in a.gsa_tar:
  sample,path=spec.split('=',1)
  with tarfile.open(path,'r:gz') as tf:
   for m in tf:
    if not m.isfile() or not m.name.endswith('_gsa.fasta.gz'):continue
    gid=genome_from_member(m.name);raw=tf.extractfile(m)
    if raw is None:continue
    with gzip.GzipFile(fileobj=raw) as z,io.TextIOWrapper(z,encoding='utf-8',errors='replace') as text:
     for orig,seq in parse_fasta_text(text):
      seq=sanitize_seq(seq);sha=seq_sha(seq);genomes[gid][sha]=(orig,seq);presence[gid].add(sample)
 manifest=[]
 with open(out/'combined_truth.fasta','w') as allf,open(out/'truth_mapping.tsv','w',newline='') as mf:
  w=csv.DictWriter(mf,fieldnames=['sequence_id','genome_id','original_id','length','sha256','samples'],delimiter='\t');w.writeheader()
  for gid in sorted(genomes):
   ref=refs/f'{safe_id(gid)}.fasta';total=0
   with open(ref,'w') as gf:
    for i,(sha,(orig,seq)) in enumerate(sorted(genomes[gid].items()),1):
     sid=f'{safe_id(gid)}|{i:05d}';write_record(gf,sid,seq);write_record(allf,sid,seq);w.writerow({'sequence_id':sid,'genome_id':gid,'original_id':orig,'length':len(seq),'sha256':sha,'samples':','.join(sorted(presence[gid]))});total+=len(seq)
   manifest.append({'genome_id':gid,'reference_path':str(ref),'truth_bp':total,'n_truth_sequences':len(genomes[gid]),'samples':','.join(sorted(presence[gid]))})
 with open(out/'genome_reference_manifest.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(manifest[0]) if manifest else ['genome_id','reference_path','truth_bp','n_truth_sequences','samples'],delimiter='\t');w.writeheader();w.writerows(manifest)
 (out/'TRUTH_SUMMARY.json').write_text(json.dumps({'n_genomes':len(genomes),'n_sequences':sum(len(v) for v in genomes.values()),'total_truth_bp':sum(r['truth_bp'] for r in manifest),'samples':sorted({s for v in presence.values() for s in v}),'truth_access_boundary':'Gold truth is used only after truth-blind assembly and selection artifacts are frozen.'},indent=2)+'\n')
if __name__=='__main__':main()
