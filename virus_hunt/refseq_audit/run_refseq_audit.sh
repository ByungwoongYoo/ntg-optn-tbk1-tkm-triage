#!/usr/bin/env bash
set -euo pipefail

PROT_QUERY="${1:?candidate protein FASTA}"
NT_QUERY="${2:?candidate contig FASTA}"
OUT="${3:?output directory}"
mkdir -p "$OUT/db" "$OUT/results" "$OUT/trees"

# RefSeq release 236 (2026-07-10) currently contains one viral protein part and
# one viral genomic FASTA part. Use an NCBI public mirror and exact filenames;
# directory-index HTML differs across servers and is not a reliable API.
BASE="https://ftp.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
PROT_FILES=("viral.1.protein.faa.gz")
NT_FILES=("viral.1.1.genomic.fna.gz")
printf '%s\n' "${PROT_FILES[@]}" > "$OUT/db/protein_files.txt"
printf '%s\n' "${NT_FILES[@]}" > "$OUT/db/nucleotide_files.txt"
printf '%s\n' "$BASE" > "$OUT/db/source_base_url.txt"

for f in "${PROT_FILES[@]}" "${NT_FILES[@]}"; do
  aria2c --continue=true --allow-overwrite=true --auto-file-renaming=false \
    --file-allocation=none --max-connection-per-server=4 --split=4 \
    --max-tries=20 --retry-wait=5 --dir "$OUT/db" --out "$f" "$BASE$f"
done

pigz -dc "$OUT"/db/viral.*.protein.faa.gz > "$OUT/db/refseq_viral_proteins.faa"
pigz -dc "$OUT"/db/viral.*.genomic.fna.gz > "$OUT/db/refseq_viral_genomes.fna"
makeblastdb -in "$OUT/db/refseq_viral_proteins.faa" -dbtype prot -parse_seqids \
  -out "$OUT/db/refseq_viral_prot" > "$OUT/results/makeblastdb_prot.log" 2>&1
makeblastdb -in "$OUT/db/refseq_viral_genomes.fna" -dbtype nucl -parse_seqids \
  -out "$OUT/db/refseq_viral_nt" > "$OUT/results/makeblastdb_nt.log" 2>&1

FMT='6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore qcovs stitle'
blastp -query "$PROT_QUERY" -db "$OUT/db/refseq_viral_prot" \
  -evalue 1e-5 -max_target_seqs 100 -seg yes -comp_based_stats 2 \
  -outfmt "$FMT" -out "$OUT/results/refseq_viral_blastp.tsv"
blastn -query "$NT_QUERY" -db "$OUT/db/refseq_viral_nt" \
  -evalue 1e-10 -max_target_seqs 100 -dust yes -soft_masking true \
  -outfmt "$FMT" -out "$OUT/results/refseq_viral_blastn.tsv"

python3 - "$OUT/results/refseq_viral_blastp.tsv" "$OUT/results/refseq_summary.tsv" <<'PY'
import csv,sys
src,dst=sys.argv[1:]
fields=['qseqid','saccver','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','qcovs','stitle']
rows=[]
with open(src) as f:
 for vals in csv.reader(f,delimiter='\t'):
  if len(vals)>=len(fields): rows.append(dict(zip(fields,vals[:len(fields)])))
by={}
for r in rows: by.setdefault(r['qseqid'],[]).append(r)
queries=['PNX_Duplo_A','PNX_Duplo_B','PNX_Picorna_A','PNX_Picorna_B']
with open(dst,'w',newline='') as f:
 outfields=['query','hit_count','top_accession','top_identity','top_qcov','top_evalue','top_bitscore','top_title','near_identity_ge90_qcov80']
 w=csv.DictWriter(f,fieldnames=outfields,delimiter='\t');w.writeheader()
 for q in queries:
  rs=by.get(q,[])
  if not rs:
   w.writerow({'query':q,'hit_count':0,'near_identity_ge90_qcov80':'false'});continue
  top=max(rs,key=lambda x:float(x['bitscore']))
  near=any(float(x['pident'])>=90 and float(x['qcovs'])>=80 for x in rs)
  w.writerow({'query':q,'hit_count':len(rs),'top_accession':top['saccver'],'top_identity':top['pident'],
              'top_qcov':top['qcovs'],'top_evalue':top['evalue'],'top_bitscore':top['bitscore'],
              'top_title':top['stitle'],'near_identity_ge90_qcov80':str(near).lower()})
PY

python3 - "$OUT/results/refseq_viral_blastp.tsv" "$OUT/trees" <<'PY'
import csv,pathlib,sys
src,out=sys.argv[1:];out=pathlib.Path(out)
fields=['qseqid','saccver','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','qcovs','stitle']
by={}
with open(src) as f:
 for vals in csv.reader(f,delimiter='\t'):
  if len(vals)<len(fields):continue
  r=dict(zip(fields,vals[:len(fields)]));by.setdefault(r['qseqid'],[]).append(r)
for q,rs in by.items():
 rs.sort(key=lambda x:float(x['bitscore']),reverse=True)
 ids=[]
 for r in rs:
  if r['saccver'] not in ids: ids.append(r['saccver'])
  if len(ids)>=20: break
 (out/f'{q}.top_ids.txt').write_text('\n'.join(ids)+'\n')
PY

for ids in "$OUT"/trees/*.top_ids.txt; do
  [[ -s "$ids" ]] || continue
  base="$(basename "$ids")"
  q="${base%.top_ids.txt}"
  blastdbcmd -db "$OUT/db/refseq_viral_prot" -entry_batch "$ids" \
    > "$OUT/trees/${q}_homologs.faa" 2> "$OUT/trees/${q}_blastdbcmd.log" || true
  python3 - "$PROT_QUERY" "$q" "$OUT/trees/${q}_query.faa" <<'PY'
from Bio import SeqIO
import sys
src,q,dst=sys.argv[1:]
found=[]
for r in SeqIO.parse(src,'fasta'):
 if r.id==q: found=[r];break
if not found: raise SystemExit(f'query {q} missing')
SeqIO.write(found,dst,'fasta')
PY
  cat "$OUT/trees/${q}_query.faa" "$OUT/trees/${q}_homologs.faa" > "$OUT/trees/${q}_combined.faa"
  mafft --auto "$OUT/trees/${q}_combined.faa" > "$OUT/trees/${q}.aln.faa" 2> "$OUT/trees/${q}.mafft.log"
  FastTree -wag "$OUT/trees/${q}.aln.faa" > "$OUT/trees/${q}.tree.nwk" 2> "$OUT/trees/${q}.fasttree.log"
done

{
  echo "generated_utc=$(date -u +%FT%TZ)"
  echo "refseq_release=236"
  echo "refseq_release_date=2026-07-10"
  echo "source_base_url=$BASE"
  blastp -version
  blastn -version
  mafft --version 2>&1 | head -n2
  FastTree 2>&1 | head -n3 || true
  du -sh "$OUT/db"/*
} > "$OUT/results/PROVENANCE.txt"

sha256sum "$OUT"/db/viral.*.gz > "$OUT/results/DOWNLOADED_DATABASE_SHA256.txt"
rm -rf "$OUT/db"
find "$OUT" -type f ! -name SHA256SUMS_COMPLETE.txt -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS_COMPLETE.txt"
