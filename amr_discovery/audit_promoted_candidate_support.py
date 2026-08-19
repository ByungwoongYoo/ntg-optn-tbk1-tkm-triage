#!/usr/bin/env python3
"""Reconstruct support and classify cross-method-promoted colistin candidates.

The output is a pre-novelty audit. It distinguishes candidates located in already
known colistin-resistance pathways from sequence markers outside those pathways.
It does not assert that any exact variant is novel or causal.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd

KNOWN_PATHWAY_GENES={
    'mcr','mgrb','phop','phoq','pmra','pmrb','pmrd','crra','crrb','crrc',
    'arna','arnb','arnc','arnd','arne','arnf','arnt','epta','ugd','pbgp',
    'lpxa','lpxc','lpxd','lpxm','pagp','acrb','rama','oqxa','oqxb',
}


def args():
    p=argparse.ArgumentParser(); p.add_argument('--consolidation',required=True); p.add_argument('--route-artifacts',required=True); p.add_argument('--upstream-root',required=True); p.add_argument('--out',required=True); return p.parse_args()

def norm_gene(value):
    return re.sub(r'[^a-z0-9]','',str(value).lower())

def candidate_gene(candidate,row):
    for key in ['gene','burden_id','feature','variant']:
        if key in row and pd.notna(row[key]):
            token=str(row[key]).split('|',1)[0].split(':',1)[0]
            if token:return token
    return str(candidate).split('|',1)[0].split(':',1)[0]

def find_rtab(root:Path,preferred='all_variants.Rtab'):
    hits=list(root.rglob(preferred));
    if hits:return hits[0]
    hits=list(root.rglob('*.Rtab'))+list(root.rglob('*.rtab'))
    return hits[0] if hits else None

def read_rtab(path):
    x=pd.read_csv(path,sep='\t',index_col=0); x.index=x.index.astype(str); x.columns=x.columns.astype(str); return x.apply(pd.to_numeric,errors='coerce').fillna(0).astype(np.uint8)

def find_component_features(candidate,raw_rows):
    sub=raw_rows[raw_rows.candidate.astype(str).eq(str(candidate))] if len(raw_rows) else pd.DataFrame()
    components=[]
    for _,r in sub.iterrows():
        if 'component_features' in r.index and pd.notna(r.component_features): components += [v for v in str(r.component_features).split(';') if v]
    return sorted(set(components))

def find_unitig_matrix(route_root,candidate):
    for p in route_root.rglob('*.rtab'):
        try:x=pd.read_csv(p,sep='\t',index_col=0,nrows=5)
        except:continue
        if str(candidate) in x.index.astype(str):return p
    for p in route_root.rglob('*.Rtab'):
        try:x=pd.read_csv(p,sep='\t',index_col=0,nrows=5)
        except:continue
        if str(candidate) in x.index.astype(str):return p
    return None

def main():
    a=args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); cons=Path(a.consolidation); route_root=Path(a.route_artifacts); upstream=Path(a.upstream_root)
    promoted_path=cons/'CROSS_METHOD_PROMOTED_CANDIDATES.csv'; raw_path=cons/'STRICT_ROUTE_ROWS.csv'
    promoted=pd.read_csv(promoted_path,dtype=str) if promoted_path.exists() and promoted_path.stat().st_size else pd.DataFrame()
    raw=pd.read_csv(raw_path,dtype=str) if raw_path.exists() and raw_path.stat().st_size else pd.DataFrame()
    manifest_hits=list(upstream.rglob('gwas_sample_manifest.csv')); rtab_hits=list(upstream.rglob('all_variants.Rtab'))
    if not manifest_hits or not rtab_hits:raise SystemExit('Frozen upstream manifest or targeted Rtab missing')
    manifest=pd.read_csv(manifest_hits[0],dtype={'assembly_ID':str}).drop_duplicates('assembly_ID'); targeted=read_rtab(rtab_hits[0])
    if promoted.empty:
        summary={'status':'NO_CROSS_METHOD_PROMOTED_CANDIDATE_TO_AUDIT','n_candidates':0,'claim_boundary':'No marker advances to novelty or context audit.'}
        pd.DataFrame().to_csv(out/'CANDIDATE_SUPPORT_MANIFEST.csv',index=False); pd.DataFrame().to_csv(out/'CANDIDATE_CLASSIFICATION.csv',index=False); (out/'CANDIDATE_SUPPORT_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2)); return
    classification=[]; support_rows=[]; fasta=[]
    for _,candrow in promoted.iterrows():
        candidate=str(candrow.candidate); cls=str(candrow.get('candidate_class','unknown')); route_rows=raw[raw.candidate.astype(str).eq(candidate)] if len(raw) else pd.DataFrame()
        gene=candidate_gene(candidate,route_rows.iloc[0].to_dict() if len(route_rows) else candrow.to_dict()); gene_norm=norm_gene(gene)
        known_pathway=gene_norm in KNOWN_PATHWAY_GENES or any(gene_norm.startswith(g) for g in ['mcr'])
        vector=None; component_features=[]; matrix_source=''
        if cls=='targeted_variant':
            if candidate in targeted.index: vector=targeted.loc[candidate].to_numpy(np.uint8); matrix_source=str(rtab_hits[0])
        elif cls=='gene_burden':
            component_features=find_component_features(candidate,raw)
            valid=[f for f in component_features if f in targeted.index]
            if valid: vector=(targeted.loc[valid].sum(axis=0)>0).astype(np.uint8).to_numpy(); matrix_source=str(rtab_hits[0])
        elif cls=='unitig':
            sequence=candidate if re.fullmatch('[ACGTNacgtn]+',candidate) else None
            if not sequence and len(route_rows):
                for col in ['canonical_sequence','variant']:
                    if col in route_rows.columns:
                        vals=[str(v) for v in route_rows[col].dropna() if re.fullmatch('[ACGTNacgtn]+',str(v))]
                        if vals:sequence=vals[0];break
            if sequence:
                fasta.append(f'>{candidate}\n{sequence.upper()}\n')
                unitig_matrix=find_unitig_matrix(route_root,sequence) or find_unitig_matrix(route_root,candidate)
                if unitig_matrix:
                    ux=read_rtab(unitig_matrix); key=sequence if sequence in ux.index else candidate
                    if key in ux.index:
                        # Map Rtab columns to frozen assembly IDs by exact/base accession.
                        exact={str(x):str(x) for x in manifest.assembly_ID}; base={str(x).split('.')[0]:str(x) for x in manifest.assembly_ID}
                        mapped=[]
                        for c in ux.columns:
                            token=Path(str(c)).name
                            for suffix in ['.fna.gz','.fasta.gz','.fa.gz','.fna','.fasta','.fa']:
                                if token.endswith(suffix):token=token[:-len(suffix)];break
                            mapped.append(exact.get(token,base.get(token.split('.')[0],token)))
                        ux.columns=mapped; common=[c for c in manifest.assembly_ID if c in ux.columns]
                        vector=pd.Series(0,index=manifest.assembly_ID,dtype=np.uint8); vector.loc[common]=ux.loc[key,common].to_numpy(np.uint8); vector=vector.to_numpy(); matrix_source=str(unitig_matrix)
        if vector is None:
            classification.append({'candidate':candidate,'candidate_class':cls,'gene':gene,'known_colistin_pathway':known_pathway,'support_reconstructed':False,'classification':'SUPPORT_MATRIX_NOT_RECONSTRUCTED','routes':candrow.get('routes')}); continue
        if len(vector)!=len(manifest):
            # Targeted matrix columns are frozen manifest IDs but may differ in order.
            if cls in {'targeted_variant','gene_burden'}:
                series=pd.Series(vector,index=targeted.columns); vector=series.reindex(manifest.assembly_ID).fillna(0).astype(np.uint8).to_numpy()
        present=vector.astype(bool); sub=manifest.loc[present].copy(); sub['candidate']=candidate; sub['candidate_class']=cls; sub['matrix_source']=matrix_source; support_rows.append(sub)
        r=int(sub.phenotype.astype(str).eq('R').sum()); s=int(sub.phenotype.astype(str).eq('S').sum()); sources=int(sub.get('source_group',pd.Series(dtype=str)).nunique(dropna=True)); countries=int(sub.get('ISO_country_code',pd.Series(dtype=str)).nunique(dropna=True))
        if known_pathway: label='KNOWN_COLISTIN_PATHWAY_EXACT_VARIANT_NOVELTY_UNRESOLVED'
        elif cls=='gene_burden': label='OUTSIDE_OR_UNRESOLVED_PATHWAY_GENE_BURDEN'
        elif cls=='unitig': label='UNITIG_CONTEXT_REQUIRED_BEFORE_PATHWAY_CLASSIFICATION'
        else: label='OUTSIDE_KNOWN_PATHWAY_TARGETED_FEATURE'
        classification.append({'candidate':candidate,'candidate_class':cls,'gene':gene,'known_colistin_pathway':known_pathway,'support_reconstructed':True,'n_present':len(sub),'n_present_R':r,'n_present_S':s,'n_source_groups':sources,'n_countries':countries,'component_features':';'.join(component_features),'classification':label,'routes':candrow.get('routes'),'lineage_supported':candrow.get('lineage_supported')})
    classifications=pd.DataFrame(classification); classifications.to_csv(out/'CANDIDATE_CLASSIFICATION.csv',index=False)
    support=pd.concat(support_rows,ignore_index=True) if support_rows else pd.DataFrame(); support.to_csv(out/'CANDIDATE_SUPPORT_MANIFEST.csv',index=False)
    if fasta:(out/'PROMOTED_UNITIGS.fasta').write_text(''.join(fasta))
    # Representative accessions are selected by candidate, phenotype and source group without sequence-context inspection.
    reps=[]
    if len(support):
        for (candidate,phenotype),z in support.groupby(['candidate','phenotype'],dropna=False):
            z=z.sort_values(['source_group','assembly_ID'] if 'source_group' in z.columns else ['assembly_ID'])
            if 'source_group' in z.columns:z=z.groupby('source_group',dropna=False).head(1)
            reps.append(z.head(20))
    representative=pd.concat(reps,ignore_index=True) if reps else pd.DataFrame(); representative.to_csv(out/'REPRESENTATIVE_CONTEXT_ACCESSIONS.csv',index=False)
    if len(representative):
        (out/'REPRESENTATIVE_ACCESSIONS.txt').write_text('\n'.join(representative.assembly_ID.dropna().astype(str).drop_duplicates())+'\n')
    summary={'status':'PROMOTED_CANDIDATES_REQUIRE_CONTEXT_AND_NOVELTY_AUDIT','n_candidates':len(classifications),'n_support_reconstructed':int(classifications.support_reconstructed.eq(True).sum()),'n_known_pathway_candidates':int(classifications.known_colistin_pathway.eq(True).sum()),'n_unitig_candidates':int(classifications.candidate_class.eq('unitig').sum()),'n_representative_accessions':len(representative),'claim_boundary':'Known-pathway candidates cannot be called new resistance genes. Unitigs require genomic-context mapping. No candidate is causal or novel at this stage.'}
    (out/'CANDIDATE_SUPPORT_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    report=['# Promoted candidate support and pathway audit','',f"- Promoted candidates: **{len(classifications):,}**",f"- Support reconstructed: **{summary['n_support_reconstructed']:,}**",f"- Located in known colistin pathways: **{summary['n_known_pathway_candidates']:,}**",f"- Unitig candidates requiring context: **{summary['n_unitig_candidates']:,}**",'',summary['claim_boundary'],'','## Candidate classification','',classifications.to_markdown(index=False)]
    (out/'CANDIDATE_SUPPORT_AUDIT_REPORT.md').write_text('\n'.join(report)+'\n')
    hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt']; (out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n'); print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
