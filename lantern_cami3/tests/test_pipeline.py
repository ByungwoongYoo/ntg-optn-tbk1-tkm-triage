#!/usr/bin/env python3
from __future__ import annotations
import csv,gzip,json,os,random,subprocess,sys,tarfile,tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/'lantern_cami3'/'scripts'
CONFIG=ROOT/'lantern_cami3'/'config'/'lantern_frozen_v1.json'
ENV={**os.environ,'PYTHONPATH':str(SCRIPTS)}

def run(*args,cwd=None):
    print('+',' '.join(map(str,args)))
    subprocess.run([str(x) for x in args],check=True,cwd=cwd,env=ENV)

def dna(n,seed):
    r=random.Random(seed);return ''.join(r.choice('ACGT') for _ in range(n))

def mutate(s,rate,seed):
    r=random.Random(seed);x=list(s);bases='ACGT'
    for i in range(len(x)):
        if r.random()<rate:x[i]=r.choice([b for b in bases if b!=x[i]])
    return ''.join(x)

def write_fasta(path,records):
    with open(path,'w') as f:
        for name,seq in records:
            f.write(f'>{name}\n')
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')

def fastq_record(name,seq):return f'@{name}\n{seq}\n+\n'+('I'*len(seq))+'\n'

def make_read_archives(work):
    short_member=work/'anonymous_short.fq.gz'
    with gzip.open(short_member,'wt') as f:
        for i in range(100):
            f.write(fastq_record(f'S0R{i}/1',dna(80,1000+i)))
            f.write(fastq_record(f'S0R{i}/2',dna(80,2000+i)))
    short_tar=work/'short.tar.gz'
    with tarfile.open(short_tar,'w:gz') as tf:tf.add(short_member,arcname='sample_0_reads/anonymous_reads.fq.gz')
    long_member=work/'anonymous_long.fq.gz'
    with gzip.open(long_member,'wt') as f:
        for i in range(50):f.write(fastq_record(f'S0R{i}',dna(1200,3000+i)))
    long_tar=work/'long.tar.gz'
    with tarfile.open(long_tar,'w:gz') as tf:tf.add(long_member,arcname='sample_0_reads/anonymous_reads.fq.gz')
    out=work/'reads';out.mkdir()
    run(sys.executable,SCRIPTS/'stream_downsample_cami_tar.py','--archive',short_tar,'--sample-id','sample0','--mode','short','--fraction','1','--out-dir',out,'--min-kept','10')
    run(sys.executable,SCRIPTS/'stream_downsample_cami_tar.py','--archive',long_tar,'--sample-id','sample0','--mode','long','--fraction','1','--out-dir',out,'--min-kept','10')
    sj=json.loads((out/'sample0_short_downsample.json').read_text());lj=json.loads((out/'sample0_long_downsample.json').read_text())
    assert sj['members'][0]['kept_records']==100 and sj['members'][1]['kept_records']==100
    assert lj['members'][0]['kept_records']==50 and sj['pair_integrity']
    assert (out/'sample0_R1.fastq.gz').is_file() and (out/'sample0_R2.fastq.gz').is_file()
    assert (out/'sample0_long.fastq.gz').is_file()

def make_assemblies(work):
    a=dna(2400,1);b=dna(2200,2);c=dna(1800,3);d=dna(1700,4);e=dna(2100,5);noise=dna(900,6)
    assemblies=work/'assemblies';assemblies.mkdir()
    write_fasta(assemblies/'single0.fa',[('A0',a),('C0',c),('short_noise',noise)])
    write_fasta(assemblies/'single1.fa',[('A1',mutate(a,.008,7))])
    write_fasta(assemblies/'co.fa',[('Aco',a),('Bco',b),('Eco',e)])
    write_fasta(assemblies/'long.fa',[('Blong',mutate(b,.005,8)),('Dlong',d)])
    manifest=work/'sources.tsv'
    with open(manifest,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['source_id','assembler','scope','mode'],delimiter='\t');w.writeheader();w.writerows([
            {'source_id':'s0','assembler':'megahit','scope':'single','mode':'short'},
            {'source_id':'s1','assembler':'megahit','scope':'single','mode':'short'},
            {'source_id':'co','assembler':'metaspades','scope':'longitudinal','mode':'short'},
            {'source_id':'lr','assembler':'flye','scope':'longitudinal','mode':'long'}])
    prep=work/'prep'
    run(sys.executable,SCRIPTS/'prepare_candidates.py','--assembly',f's0={assemblies/"single0.fa"}','--assembly',f's1={assemblies/"single1.fa"}','--assembly',f'co={assemblies/"co.fa"}','--assembly',f'lr={assemblies/"long.fa"}','--source-manifest',manifest,'--min-length','1000','--out',prep)
    seqs={}
    name=None;buf=[]
    for line in (prep/'combined_candidates.fasta').read_text().splitlines():
        if line.startswith('>'):
            if name:seqs[name]=''.join(buf)
            name=line[1:];buf=[]
        else:buf.append(line)
    if name:seqs[name]=''.join(buf)
    paf=work/'self.paf'
    with open(paf,'w') as f:
        ids=list(seqs)
        for q in ids:
            for t in ids:
                if q==t:continue
                l=min(len(seqs[q]),len(seqs[t]));matches=sum(x==y for x,y in zip(seqs[q][:l],seqs[t][:l]));ident=matches/l
                if ident>=.97:f.write(f'{q}\t{len(seqs[q])}\t0\t{l}\t+\t{t}\t{len(seqs[t])}\t0\t{l}\t{matches}\t{l}\t60\n')
    clusters=work/'clusters'
    run(sys.executable,SCRIPTS/'cluster_candidates.py','--fasta',prep/'combined_candidates.fasta','--metadata',prep/'candidates.tsv','--paf',paf,'--out',clusters)
    return clusters,a,b,c,d,e

def evidence_and_selection(work,clusters,e_seq):
    reps=list(csv.DictReader(open(clusters/'representative_metadata.tsv'),delimiter='\t'))
    members=list(csv.DictReader(open(clusters/'cluster_members.tsv'),delimiter='\t'))
    sources={r['representative_id']:set(r['source_ids'].split(',')) for r in reps}
    eco_candidate=next(r['candidate_id'] for r in members if r['source_id']=='co' and r['original_id']=='Eco')
    e_rep=next(r['representative_id'] for r in members if r['candidate_id']==eco_candidate)
    ev=work/'evidence.tsv';fields=['contig_id','sample_id','short_reads','short_breadth','short_depth','long_reads','long_breadth','long_depth','long_spanning_reads']
    with open(ev,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for r in reps:
            cid=r['representative_id'];src=sources[cid]
            for sample in ('t0','t1'):
                base={'contig_id':cid,'sample_id':sample,'short_reads':100,'short_breadth':.9,'short_depth':5,'long_reads':0,'long_breadth':0,'long_depth':0,'long_spanning_reads':0}
                if src=={'lr'}:
                    base.update({'short_breadth':.2,'short_depth':.2,'long_reads':10,'long_breadth':.9,'long_depth':2,'long_spanning_reads':2 if sample=='t0' else 0})
                if cid==e_rep:
                    base.update({'short_breadth':.8,'short_depth':1.2,'long_breadth':0,'long_depth':0,'long_spanning_reads':0})
                w.writerow(base)
    full=work/'select_full';nol=work/'select_no_longitudinal'
    run(sys.executable,SCRIPTS/'lantern_select.py','--fasta',clusters/'cluster_representatives.fasta','--metadata',clusters/'representative_metadata.tsv','--evidence',ev,'--config',CONFIG,'--out',full,'--ablation','full')
    run(sys.executable,SCRIPTS/'lantern_select.py','--fasta',clusters/'cluster_representatives.fasta','--metadata',clusters/'representative_metadata.tsv','--evidence',ev,'--config',CONFIG,'--out',nol,'--ablation','no_longitudinal')
    full_rows={r['contig_id']:r for r in csv.DictReader(open(full/'LANTERN_SELECTION.tsv'),delimiter='\t')}
    no_rows={r['contig_id']:r for r in csv.DictReader(open(nol/'LANTERN_SELECTION.tsv'),delimiter='\t')}
    assert full_rows[e_rep]['selected']=='true' and full_rows[e_rep]['rescue']=='true'
    assert no_rows[e_rep]['selected']=='false', (full_rows[e_rep],no_rows[e_rep])
    run(sys.executable,SCRIPTS/'verify_cami_assembly.py',full/'LANTERN_ASSEMBLY.fasta','--out',full/'VALIDATION.json')
    return full

def scaffold_test(work):
    x=dna(2000,20);y=dna(1900,21);fa=work/'bridge_input.fa';write_fasta(fa,[('x',x),('y',y)])
    paf=work/'bridges.paf'
    with open(paf,'w') as f:
        for i in range(3):
            q=f'read{i}';f.write(f'{q}\t4000\t0\t1800\t+\tx\t2000\t200\t2000\t1780\t1800\t60\n');f.write(f'{q}\t4000\t2000\t3800\t+\ty\t1900\t0\t1800\t1780\t1800\t60\n')
    out=work/'scaffold'
    run(sys.executable,SCRIPTS/'conservative_scaffold.py','--fasta',fa,'--paf',paf,'--config',CONFIG,'--out',out)
    summary=json.loads((out/'SCAFFOLD_SUMMARY.json').read_text());assert summary['n_chosen_edges']==1
    run(sys.executable,SCRIPTS/'verify_cami_assembly.py',out/'LANTERN_SCAFFOLDED.fasta','--out',out/'VALIDATION.json')

def gold_test(work,assembly):
    seqs=[];name=None;buf=[]
    for line in assembly.read_text().splitlines():
        if line.startswith('>'):
            if name:seqs.append((name,''.join(buf)))
            name=line[1:];buf=[]
        else:buf.append(line)
    if name:seqs.append((name,''.join(buf)))
    truth=work/'truth.fa';mapping=work/'truth_mapping.tsv'
    chosen=seqs[:2];write_fasta(truth,[('g1|00001',chosen[0][1]),('g2|00001',chosen[1][1])])
    with open(mapping,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sequence_id','genome_id','original_id','length','sha256','samples'],delimiter='\t');w.writeheader();
        for i,(_,s) in enumerate(chosen,1):w.writerow({'sequence_id':f'g{i}|00001','genome_id':f'g{i}','original_id':f'o{i}','length':len(s),'sha256':'x','samples':'t0,t1'})
    paf=work/'gold.paf'
    with open(paf,'w') as f:
        for i,(q,s) in enumerate(chosen,1):f.write(f'{q}\t{len(s)}\t0\t{len(s)}\t+\tg{i}|00001\t{len(s)}\t0\t{len(s)}\t{len(s)}\t{len(s)}\t60\n')
    out=work/'gold_eval';run(sys.executable,SCRIPTS/'evaluate_gold_coverage.py','--paf',paf,'--truth-mapping',mapping,'--assembly-fasta',assembly,'--min-alignment','500','--out',out)
    summary=json.loads((out/'GOLD_COVERAGE_SUMMARY.json').read_text());assert summary['genomes_recovered_90']==2
    report=work/'report.tsv';report.write_text('Assembly\tmethod\nGenome fraction (%)\t12.5\nDuplication ratio\t1.02\n# misassemblies\t2\nN50\t2000\n')
    run(sys.executable,SCRIPTS/'parse_metaquast.py','--report',report,'--method','fixture','--out',work/'parsed.json')
    assert json.loads((work/'parsed.json').read_text())['genome_fraction']==12.5

def main():
    with tempfile.TemporaryDirectory(prefix='lantern-test-') as td:
        work=Path(td);make_read_archives(work);clusters,*_=make_assemblies(work);full=evidence_and_selection(work,clusters,None);scaffold_test(work);gold_test(work,full/'LANTERN_ASSEMBLY.fasta')
    print('LANTERN_SYNTHETIC_INTEGRATION_PASS')
if __name__=='__main__':main()
