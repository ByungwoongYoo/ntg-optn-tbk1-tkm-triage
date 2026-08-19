#!/usr/bin/env python3
from __future__ import annotations
import argparse,hashlib,json,shutil,time
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--method',required=True);p.add_argument('--scope',choices=['global_cross_sample','participant_cross_sample','single_sample'],required=True);p.add_argument('--data-mode',choices=['short','long','hybrid'],required=True);p.add_argument('--samples',required=True);p.add_argument('--versions',required=True);p.add_argument('--parameters',required=True);p.add_argument('--git-commit',required=True);p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);src=Path(a.fasta);dest=out/'assembly.fasta';shutil.copy2(src,dest)
    metadata={'method_name':a.method,'assembly_scope':a.scope,'data_mode':a.data_mode,'sample_ids':[x.strip() for x in a.samples.split(',') if x.strip()],'git_commit':a.git_commit,'software_versions_file':Path(a.versions).name,'parameters_file':Path(a.parameters).name,'assembly_sha256':sha(dest),'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'claim_boundary':'Submission-ready files only. No official CAMI score or rank exists until the portal accepts and evaluates the submission.'}
    shutil.copy2(a.versions,out/Path(a.versions).name);shutil.copy2(a.parameters,out/Path(a.parameters).name);(out/'submission_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n');(out/'README.md').write_text(f"# CAMI III assembly submission bundle\n\n- Method: {a.method}\n- Scope: {a.scope}\n- Data mode: {a.data_mode}\n- Samples: {a.samples}\n- Git commit: {a.git_commit}\n- FASTA SHA-256: {metadata['assembly_sha256']}\n\nThe FASTA must be uploaded through the CAMI portal by the human participant. Portal acceptance and evaluation are external actions and are not claimed here.\n")
    files=[]
    for x in sorted(out.iterdir()):
        if x.is_file():files.append({'path':x.name,'bytes':x.stat().st_size,'sha256':sha(x)})
    (out/'SHA256SUMS.json').write_text(json.dumps(files,indent=2)+'\n');print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()
