#!/usr/bin/env python3
"""Independently verify targeted substitutions and unitig genomic contexts.

Targeted substitutions are recalled from tblastn alignments of pinned reference
proteins against representative assemblies. Unitig hits are recalled with blastn
and flanking sequence is extracted. This checks sequence support and assembly
context only; it does not establish phenotype causality.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,subprocess,tempfile
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def args():
    p=argparse.ArgumentParser();p.add_argument('--classification',required=True);p.add_argument('--support',required=True);p.add_argument('--representatives',required=True);p.add_argument('--assemblies-dir',required=True);p.add_argument('--reference-gbff',required=True);p.add_argument('--unitigs',default='');p.add_argument('--out',required=True);p.add_argument('--flank',type=int,default=5000);return p.parse_args()

def norm_gene(x):return re.sub(r'[^a-z0-9]','',str(x).lower())

def parse_mutation(candidate):
    m=re.search(r'([A-Z*])([0-9]+)([A-Z*])',str(candidate).upper());return (m.group(1),int(m.group(2)),m.group(3)) if m else None

def reference_proteins(gbff):
    records=list(SeqIO.parse(gbff,'genbank'));out={}
    for rec in records:
        for f in rec.features:
            if f.type!='CDS':continue
            gene=(f.qualifiers.get('gene') or f.qualifiers.get('locus_tag') or [''])[0]
            translation=(f.qualifiers.get('translation') or [''])[0]
            if gene and translation:out[norm_gene(gene)]={'gene':gene,'protein':translation,'locus_tag':(f.qualifiers.get('locus_tag') or [''])[0]}
    return out

def find_assembly(root,acc):
    root=Path(root);hits=list(root.glob(f'{acc}*.fna'))+list(root.glob(f'{acc.split(".")[0]}*.fna'));return hits[0] if hits else None

def run(cmd):
    return subprocess.run(cmd,text=True,capture_output=True,check=False)

def tblastn_recall(protein,gene,candidate,assembly,tmpdir):
    mut=parse_mutation(candidate)
    q=Path(tmpdir)/f'{re.sub(r"[^A-Za-z0-9_.-]","_",candidate)}.faa';q.write_text(f'>{gene}\n{protein}\n')
    fmt='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq'
    r=run(['tblastn','-query',str(q),'-subject',str(assembly),'-evalue','1e-20','-max_target_seqs','5','-outfmt',fmt])
    if r.returncode!=0:return {'status':'TBLASTN_ERROR','stderr':r.stderr[:2000]}
    rows=[]
    for line in r.stdout.splitlines():
        f=line.split('\t');
        if len(f)<12:continue
        rows.append({'qseqid':f[0],'sseqid':f[1],'pident':float(f[2]),'length':int(f[3]),'qstart':int(f[4]),'qend':int(f[5]),'sstart':int(f[6]),'send':int(f[7]),'evalue':float(f[8]),'bitscore':float(f[9]),'qseq':f[10],'sseq':f[11]})
    if not rows:return {'status':'NO_HIT'}
    hit=max(rows,key=lambda x:x['bitscore']);result={'status':'HIT','hit':{k:v for k,v in hit.items() if k not in ['qseq','sseq']}}
    if mut:
        ref,pos,alt=mut;query_position=hit['qstart'];observed=None;aligned_ref=None
        for qa,sa in zip(hit['qseq'],hit['sseq']):
            if qa!='-':
                if query_position==pos:aligned_ref=qa;observed=sa;break
                query_position+=1
        result.update({'expected_reference_residue':ref,'reference_position':pos,'expected_alternate_residue':alt,'aligned_reference_residue':aligned_ref,'observed_assembly_residue':observed,'mutation_confirmed':observed==alt and aligned_ref==ref})
    return result

def load_unitigs(path):
    if not path or not Path(path).exists():return {}
    return {rec.id:str(rec.seq).upper() for rec in SeqIO.parse(path,'fasta')}

def blast_unitig(candidate,seq,assembly,outdir,flank):
    q=Path(outdir)/'query.fasta';q.write_text(f'>{candidate}\n{seq}\n')
    fmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    r=run(['blastn','-task','blastn-short','-query',str(q),'-subject',str(assembly),'-evalue','1e-10','-max_target_seqs','20','-outfmt',fmt])
    if r.returncode!=0:return {'status':'BLASTN_ERROR','stderr':r.stderr[:2000]},None
    hits=[]
    for line in r.stdout.splitlines():
        f=line.split('\t');
        if len(f)<12:continue
        hits.append({'qseqid':f[0],'sseqid':f[1],'pident':float(f[2]),'length':int(f[3]),'mismatch':int(f[4]),'gapopen':int(f[5]),'qstart':int(f[6]),'qend':int(f[7]),'sstart':int(f[8]),'send':int(f[9]),'evalue':float(f[10]),'bitscore':float(f[11])})
    if not hits:return {'status':'NO_HIT'},None
    best=max(hits,key=lambda x:(x['bitscore'],x['length'],x['pident']));records=SeqIO.to_dict(SeqIO.parse(assembly,'fasta'));rec=records.get(best['sseqid'])
    if rec is None:return {'status':'HIT_CONTIG_NOT_FOUND','best':best},None
    lo=min(best['sstart'],best['send'])-1;hi=max(best['sstart'],best['send']);start=max(0,lo-flank);end=min(len(rec.seq),hi+flank);context=rec[start:end]
    edge_distance=min(lo,len(rec.seq)-hi);reverse=best['sstart']>best['send']
    if reverse:context=context.reverse_complement()
    header=f'{candidate}|{Path(assembly).stem}|{best["sseqid"]}:{start+1}-{end}|reverse={int(reverse)}'
    return {'status':'HIT','best':best,'contig_length':len(rec.seq),'context_start':start+1,'context_end':end,'distance_to_nearest_contig_edge':edge_distance,'context_length':len(context),'reverse_complemented':reverse,'perfect_full_length':best['pident']==100.0 and best['length']==len(seq)},(header,str(context))

def main():
    a=args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);classification=pd.read_csv(a.classification,dtype=str) if Path(a.classification).exists() else pd.DataFrame();support=pd.read_csv(a.support,dtype=str) if Path(a.support).exists() else pd.DataFrame();reps=pd.read_csv(a.representatives,dtype=str) if Path(a.representatives).exists() else pd.DataFrame();refs=reference_proteins(a.reference_gbff);unitigs=load_unitigs(a.unitigs);variant_rows=[];unitig_rows=[];contexts=[]
    if classification.empty:
        summary={'status':'NO_CANDIDATE_TO_SEQUENCE_VERIFY','n_targeted_calls':0,'n_unitig_contexts':0};(out/'SEQUENCE_VERIFICATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));return
    with tempfile.TemporaryDirectory() as td:
        for _,c in classification.iterrows():
            cand=str(c.candidate);cls=str(c.candidate_class);gene=str(c.get('gene',''));candidate_reps=reps[reps.candidate.astype(str).eq(cand)] if len(reps) and 'candidate' in reps.columns else pd.DataFrame()
            if cls=='targeted_variant':
                ref=refs.get(norm_gene(gene))
                for _,s in candidate_reps.iterrows():
                    acc=str(s.assembly_ID);assembly=find_assembly(a.assemblies_dir,acc)
                    if not ref or assembly is None:res={'status':'REFERENCE_OR_ASSEMBLY_MISSING'}
                    else:res=tblastn_recall(ref['protein'],gene,cand,assembly,td)
                    variant_rows.append({'candidate':cand,'gene':gene,'assembly_ID':acc,'phenotype':s.get('phenotype'),'source_group':s.get('source_group'),'country':s.get('ISO_country_code'),**res})
            elif cls=='unitig':
                seq=unitigs.get(cand) or (cand.upper() if re.fullmatch('[ACGTNacgtn]+',cand) else None)
                for _,s in candidate_reps.iterrows():
                    acc=str(s.assembly_ID);assembly=find_assembly(a.assemblies_dir,acc);cdir=out/'unitig_tmp'/re.sub(r'[^A-Za-z0-9_.-]','_',cand)/acc;cdir.mkdir(parents=True,exist_ok=True)
                    if not seq or assembly is None:res,context={'status':'SEQUENCE_OR_ASSEMBLY_MISSING'},None
                    else:res,context=blast_unitig(cand,seq,assembly,cdir,a.flank)
                    unitig_rows.append({'candidate':cand,'assembly_ID':acc,'phenotype':s.get('phenotype'),'source_group':s.get('source_group'),'country':s.get('ISO_country_code'),**res})
                    if context:contexts.append(context)
    pd.DataFrame(variant_rows).to_csv(out/'TARGETED_VARIANT_DIRECT_RECALL.csv',index=False);pd.DataFrame(unitig_rows).to_csv(out/'UNITIG_DIRECT_CONTEXT_RECALL.csv',index=False)
    if contexts:
        with open(out/'UNITIG_FLANKING_CONTEXTS.fasta','w') as f:
            for h,s in contexts:f.write(f'>{h}\n{s}\n')
    if (out/'unitig_tmp').exists():
        import shutil;shutil.rmtree(out/'unitig_tmp')
    vr=pd.DataFrame(variant_rows);ur=pd.DataFrame(unitig_rows)
    summary={'status':'SEQUENCE_SUPPORT_RECALLED_REQUIRES_MECHANISM_INTERPRETATION','n_targeted_calls':len(vr),'n_targeted_confirmed':int(vr.get('mutation_confirmed',pd.Series(dtype=bool)).eq(True).sum()) if len(vr) else 0,'n_targeted_ambiguous_or_failed':int(len(vr)-vr.get('mutation_confirmed',pd.Series([False]*len(vr))).eq(True).sum()) if len(vr) else 0,'n_unitig_calls':len(ur),'n_unitig_hits':int(ur.get('status',pd.Series(dtype=str)).eq('HIT').sum()) if len(ur) else 0,'n_unitig_perfect_full_length':int(ur.get('perfect_full_length',pd.Series(dtype=bool)).eq(True).sum()) if len(ur) else 0,'n_unitig_contexts':len(contexts),'n_unitig_near_contig_edge_lt_500':int((pd.to_numeric(ur.get('distance_to_nearest_contig_edge',pd.Series(dtype=float)),errors='coerce')<500).sum()) if len(ur) else 0,'claim_boundary':'Direct sequence recall checks assembly support only. It does not establish resistance causality, biological novelty, or clinical validity.'}
    (out/'SEQUENCE_VERIFICATION_SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    report=['# Candidate direct sequence verification','',f"- Targeted calls: **{summary['n_targeted_calls']:,}**; confirmed: **{summary['n_targeted_confirmed']:,}**",f"- Unitig calls: **{summary['n_unitig_calls']:,}**; hits: **{summary['n_unitig_hits']:,}**",f"- Unitig contexts recovered: **{summary['n_unitig_contexts']:,}**",f"- Unitig hits <500 bp from contig edge: **{summary['n_unitig_near_contig_edge_lt_500']:,}**",'',summary['claim_boundary']]
    (out/'SEQUENCE_VERIFICATION_REPORT.md').write_text('\n'.join(report)+'\n');hashes=[f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}" for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'];(out/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n');print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__':main()
