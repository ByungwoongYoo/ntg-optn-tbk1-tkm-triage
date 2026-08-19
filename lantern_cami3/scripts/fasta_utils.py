from __future__ import annotations
import gzip, hashlib, re
from typing import Iterator, TextIO

VALID=set('ACGTN')

def open_text(path, mode='rt'):
    p=str(path)
    return gzip.open(p,mode,encoding=None if 'b' in mode else 'utf-8') if p.endswith('.gz') else open(p,mode,encoding=None if 'b' in mode else 'utf-8')

def read_fasta(path) -> Iterator[tuple[str,str]]:
    with open_text(path,'rt') as f:
        name=None; seq=[]
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None: yield name,''.join(seq)
                name=line[1:].strip().split()[0]; seq=[]
            else: seq.append(line)
        if name is not None: yield name,''.join(seq)

def sanitize_seq(seq:str)->str:
    s=seq.upper().replace('U','T')
    return ''.join(c if c in VALID else 'N' for c in s)

def safe_id(text:str)->str:
    s=re.sub(r'[^A-Za-z0-9_.:|,;\[\]-]+','_',text)
    return s[:180] or 'sequence'

def write_record(f:TextIO,name:str,seq:str,width:int=80):
    f.write('>'+safe_id(name)+'\n')
    for i in range(0,len(seq),width): f.write(seq[i:i+width]+'\n')

def seq_sha(seq:str)->str: return hashlib.sha256(seq.encode()).hexdigest()

def revcomp(seq:str)->str: return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
