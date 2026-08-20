#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,re,hashlib
from pathlib import Path
from fasta_utils import read_fasta
ALLOWED=set('ATCGN'); HEADER=re.compile(r'^[A-Za-z0-9\s\[\]_:;,\.\|\-]+$')
def main():
 p=argparse.ArgumentParser();p.add_argument('fasta');p.add_argument('--out',required=True);a=p.parse_args();ids=set();bad=[];n=bp=0
 for name,seq in read_fasta(a.fasta):
  n+=1;bp+=len(seq)
  if name in ids:bad.append(f'duplicate:{name}')
  ids.add(name)
  if not HEADER.match(name):bad.append(f'header:{name}')
  if set(seq)-ALLOWED:bad.append(f'alphabet:{name}:{sorted(set(seq)-ALLOWED)}')
  if not seq:bad.append(f'empty:{name}')
 status={'status':'PASS' if not bad and n else 'FAIL','n_sequences':n,'total_bp':bp,'n_errors':len(bad),'errors':bad[:100],'sha256':hashlib.sha256(Path(a.fasta).read_bytes()).hexdigest()}
 Path(a.out).write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2));raise SystemExit(0 if status['status']=='PASS' else 1)
if __name__=='__main__':main()
