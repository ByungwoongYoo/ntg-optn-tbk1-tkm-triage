#!/usr/bin/env python3
"""Build public-data A24 candidate models for CASP17 T2461.

The full-length public ESMFold monomer is placed on each subunit of public
3VCD/6FDB biological assemblies. Residues 9-191 of T2461 align without gaps
to template residues 1-183 (Foldseek: 63.3% identity to 3VCD). No unreleased
experimental target coordinates are used.
"""
import argparse,csv,hashlib,json,math,re
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import cKDTree
SEQ='MHHHHHHGSMAIGIVELSSIAMGLKLADEMLKAADVKLLVSRPILPGKFLIILGGETEAIRKAIAVATEAAGSKLVRSALIEDIHPSVLPAISGINPVEERQAVGIVETESLEAAILAANAAVKGSNVTLVRIRMLSGITGKCYIVVAGDVDDVALAVVVAAEVAASRGKLIYAALIPRPHPAIWPLIVEG'
CID='ABCDEFGHIJKLMNOPQRSTUVWX'
AA={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
def sh(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def monomer(path):
 atoms=[]; ca={}; rn={}
 for l in Path(path).read_text(errors='replace').splitlines():
  if not l.startswith('ATOM') or l[16:17] not in (' ','A'):continue
  n=l[12:16]; r=int(l[22:26]); x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
  a={'name':n,'resname':l[17:20].strip(),'resid':r,'xyz':x,'occ':float(l[54:60]),'bf':float(l[60:66]),'el':l[76:78].strip() or re.sub('[^A-Za-z]','',n)[:1]}
  atoms.append(a);rn[r]=a['resname'];
  if n.strip()=='CA':ca[r]=x
 s=''.join(AA.get(rn.get(i),'X') for i in range(1,192))
 assert s==SEQ and len(ca)==191,(s,len(ca))
 return atoms,ca
def chains(cif):
 st=MMCIFParser(QUIET=True).get_structure('t',cif);m=next(st.get_models());out=[]
 for c in m:
  rs=[r for r in c if is_aa(r,standard=False) and 'CA' in r]
  if len(rs)>=183:out.append((str(c.id),rs,np.array([r['CA'].coord for r in rs])))
 assert len(out)==24,(cif,len(out))
 return out
def fit(P,Q):
 pc=P.mean(0);qc=Q.mean(0);u,s,vt=np.linalg.svd((P-pc).T@(Q-qc));R=vt.T@u.T
 if np.linalg.det(R)<0:vt[-1]*=-1;R=vt.T@u.T
 t=qc-R@pc;z=(R@P.T).T+t
 return R,t,float(np.sqrt(np.mean(np.sum((z-Q)**2,axis=1))))
def build(atoms,ca,cif,delta):
 out=[]; rms=[]; cents=[]
 for k,(tc,rs,Qall) in enumerate(chains(cif)):
  Q=Qall[:183];P=np.array([ca[i] for i in range(9,192)]);R,t,r=fit(P,Q);rms.append(r)
  aa=[]
  for a in atoms:aa.append({**a,'chain':CID[k],'xyz':R@a['xyz']+t})
  cents.append(np.mean([a['xyz'] for a in aa if a['name'].strip()=='CA'],0));out.append(aa)
 center=np.mean(cents,0)
 for aa,c in zip(out,cents):
  v=c-center;v=v/np.linalg.norm(v)*delta
  for a in aa:a['xyz']=a['xyz']+v
 return out,float(np.mean(rms))
def metrics(ch):
 xyz=[np.array([a['xyz'] for a in c]) for c in ch];sev=c15=c20=ct=0;mind=999.;pairs=0
 for i in range(24):
  ti=cKDTree(xyz[i])
  for j in range(i+1,24):
   d=ti.sparse_distance_matrix(cKDTree(xyz[j]),4.5,output_type='coo_matrix').data
   if len(d):mind=min(mind,float(d.min()));sev+=int((d<1.2).sum());c15+=int((d<1.5).sum());c20+=int((d<2).sum());ct+=len(d);pairs+=1
 cent=np.array([np.mean([a['xyz'] for a in c if a['name'].strip()=='CA'],0) for c in ch]);rad=np.linalg.norm(cent-cent.mean(0),axis=1)
 return {'severe_lt1.2':sev,'clash_lt1.5':c15,'close_lt2.0':c20,'contacts_lt4.5':ct,'contacting_pairs':pairs,'min_interchain_A':mind,'radius_mean_A':float(rad.mean()),'radius_sd_A':float(rad.std())}
def write(ch,p):
 s=1
 with open(p,'w') as f:
  f.write('HEADER    CASP17 PREDICTION\nTITLE     T2461 DESIGNED PROTEIN CAGE A24 CANDIDATE\n')
  for c in ch:
   for a in c:
    x=a['xyz'];f.write(f"ATOM  {s:5d} {a['name']}{' ':1s}{a['resname']:>3s} {a['chain']}{a['resid']:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{a['occ']:6.2f}{a['bf']:6.2f}          {a['el']:>2s}\n");s+=1
   f.write(f"TER   {s:5d}      GLY {c[0]['chain']} 191\n");s+=1
  f.write('END\n')
def submission(p,out,n,parent,delta):
 lines=[x for x in Path(p).read_text().splitlines() if x.startswith(('ATOM','TER'))]
 with open(out,'w') as f:
  f.write('PFRMAT TS\nTARGET T2461\nAUTHOR REPLACE_WITH_REGISTERED_CASP_GROUP_CODE\n')
  f.write('REMARK NOT OFFICIALLY SUBMITTED; REPLACE AUTHOR BEFORE UPLOAD.\n')
  f.write('METHOD Public ESMFold v1 full-length monomer prediction.\n')
  f.write(f'METHOD Rigid placement on PDB {parent} biological assembly 1, A24.\n')
  f.write(f'METHOD Radial interface alternative {delta:.2f} Angstrom.\nMODEL {n}\nPARENT {parent}\n')
  f.write('\n'.join(lines)+'\nEND\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--monomer',required=True);ap.add_argument('--template',action='append',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
 at,ca=monomer(a.monomer);all=[]
 for spec in a.template:
  pid,cif=spec.split('=',1)
  for d in (0.,.25,.5,.75,1.):
   try:
    ch,r=build(at,ca,cif,d);m=metrics(ch);p=o/f'candidate_{pid}_{d:.2f}.pdb';write(ch,p)
    rec={'template':pid.upper(),'delta_A':d,'fit_rmsd_A':r,'file':p.name,'sha256':sh(p),'metrics':m,'chains':24,'residues_per_chain':191,'atoms':sum(map(len,ch))};rec['rank']=[m['severe_lt1.2'],m['clash_lt1.5'],m['close_lt2.0'],-m['contacts_lt4.5'],d];all.append(rec)
   except Exception as e:all.append({'template':pid.upper(),'delta_A':d,'error':repr(e)})
 good=[x for x in all if 'error' not in x];assert good,all
 # top 3 from 3VCD and top 2 from 6FDB, then clash-rank globally
 pick=[]
 for pid,n in [('3VCD',3),('6FDB',2)]:pick += sorted([x for x in good if x['template']==pid],key=lambda z:z['rank'])[:n]
 pick=sorted(pick,key=lambda z:z['rank'])[:5]
 for i,x in enumerate(pick,1):
  src=o/x['file'];dst=o/f'T2461_model{i}.pdb';dst.write_bytes(src.read_bytes());x['model']=i;x['selected_file']=dst.name;x['selected_sha256']=sh(dst);submission(dst,o/f'T2461_model{i}_submission_TEMPLATE.txt',i,x['template'],x['delta_A'])
 summary={'schema':'casp17-t2461-a24-v1','target':'T2461','sequence':SEQ,'stoichiometry':'A24','monomer_sha256':sh(a.monomer),'candidates':all,'selected':pick,'claim_boundary':'Candidate public-data predictions; not experimental structures or official CASP submissions.'}
 (o/'MODEL_BUILD_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 with open(o/'MODEL_VALIDATION.csv','w',newline='') as f:
  w=csv.writer(f);w.writerow(['model','template','delta_A','fit_rmsd_A','severe_lt1.2','clash_lt1.5','close_lt2.0','contacts_lt4.5','sha256'])
  for x in pick:w.writerow([x['model'],x['template'],x['delta_A'],x['fit_rmsd_A'],x['metrics']['severe_lt1.2'],x['metrics']['clash_lt1.5'],x['metrics']['close_lt2.0'],x['metrics']['contacts_lt4.5'],x['selected_sha256']])
 print(json.dumps({'selected':pick},indent=2))
if __name__=='__main__':main()
