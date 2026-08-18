#!/usr/bin/env python3
"""Convert T2461 A24 models to CASP17-compliant TS text templates.

Fixes confidence values to the required 0-100 B-factor scale, emits 80-column
PDB ATOM/TER records, adds STOICH A24 and chain-qualified PARENT records, and
builds both single-model and five-model files. AUTHOR remains a conspicuous
placeholder until a real CASP17 regular-group registration code is supplied.
"""
import argparse,hashlib,json
from pathlib import Path
PARENTS={1:'3VCD_A',2:'6FDB_A',3:'3VCD_A',4:'3VCD_A',5:'6FDB_A'}
def sh(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom80(line):
 x=float(line[30:38]);y=float(line[38:46]);z=float(line[46:54]);occ=float(line[54:60]);bf=float(line[60:66]);bf=bf*100 if bf<=1.01 else bf
 out=f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}{line[66:78]}"
 return out[:80].ljust(80)
def ter80(line):return line[:80].ljust(80)
def model_lines(path):
 out=[]
 for l in Path(path).read_text().splitlines():
  if l.startswith('ATOM'):out.append(atom80(l))
  elif l.startswith('TER'):out.append(ter80(l))
 return out
def header():return ['PFRMAT TS','TARGET T2461','AUTHOR REPLACE_WITH_REGISTERED_CASP_GROUP_CODE','REMARK PRE-SUBMISSION TEMPLATE; NOT OFFICIALLY SUBMITTED.','METHOD Public ESMFold v1 full-length monomer prediction.','METHOD Foldseek identified PduT-based cage templates; target threaded rigidly.','METHOD A24 geometry from public RCSB 3VCD/6FDB biological assemblies.']
def block(n,lines):return [f'MODEL {n}', 'STOICH A24', f'PARENT {PARENTS[n]}']+lines+['END']
def validate(lines,n):
 atoms=[x for x in lines if x.startswith('ATOM')];ters=[x for x in lines if x.startswith('TER')];chains={x[21] for x in atoms};b=[float(x[60:66]) for x in atoms]
 return {'model':n,'atoms':len(atoms),'ter_records':len(ters),'chains':len(chains),'all_coordinate_records_80_columns':all(len(x)==80 for x in atoms+ters),'bfactor_min':min(b),'bfactor_max':max(b),'bfactor_unique':len(set(b)),'bfactor_percentage_scale':max(b)>1 and max(b)<=100 and min(b)>=0}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--models-dir',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();src=Path(a.models_dir);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 allout=header();audit=[]
 for n in range(1,6):
  lines=model_lines(src/f'T2461_model{n}.pdb');v=validate(lines,n);assert v['atoms']==33528 and v['ter_records']==24 and v['chains']==24 and v['all_coordinate_records_80_columns'] and v['bfactor_percentage_scale'],v
  p=out/f'T2461_model{n}_CASP17_AUTHOR_CODE_REQUIRED.txt';content=header()+block(n,lines);p.write_text('\n'.join(content)+'\n');v['file']=p.name;v['sha256']=sh(p);v['parent']=PARENTS[n];audit.append(v);allout+=block(n,lines)
 combo=out/'T2461_models1-5_CASP17_AUTHOR_CODE_REQUIRED.txt';combo.write_text('\n'.join(allout)+'\n')
 report={'schema':'casp17-t2461-ts-format-audit-v1','target':'T2461','stoichiometry':'A24','official_submission':False,'credential_gate':'Replace AUTHOR with a real 12-digit registered CASP17 regular-group code, log in, and submit through the regular CASP model submission form.','models':audit,'combined_file':combo.name,'combined_sha256':sh(combo),'combined_model_count':5,'claim_boundary':'Submission-ready format template; not accepted or officially submitted.'}
 (out/'CASP17_FORMAT_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
