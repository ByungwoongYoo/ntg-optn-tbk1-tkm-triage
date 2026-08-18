#!/usr/bin/env bash
set -euo pipefail

PROT_QUERY="${1:?candidate protein FASTA}"
NT_QUERY="${2:?candidate contig FASTA}"
OUT="${3:?output directory}"
mkdir -p "$OUT/db" "$OUT/results" "$OUT/trees"

BASE="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
curl -fsSL --retry 10 "$BASE" -o "$OUT/db/viral_index.html"
mapfile -t PROT_FILES < <(grep -oE 'href="viral\.[0-9]+\.protein\.faa\.gz"' "$OUT/db/viral_index.html" | cut -d'"' -f2 | sort -u)
mapfile -t NT_FILES < <(grep -oE 'href="viral\.[0-9]+\.genomic\.fna\.gz"' "$OUT/db/viral_index.html" | cut -d'"' -f2 | sort -u)
if [[ ${#PROT_FILES[@]} -eq 0 || ${#NT_FILES[@]} -eq 0 ]]; then
  echo "Could not parse RefSeq viral release index" >&2
  exit 2
fi
printf '%s\n' "${PROT_FILES[@]}" > "$OUT/db/protein_files.txt"
printf '%s\n' "${NT_FILES[@]}" > "$OUT/db/nucleotide_files.txt"

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
  if len(vals)>=len(fields):rows.append(dict(zip(fields,vals[:len(fields)])))
by={}
for r in rows:
 q=r['qseqid']; by.setdefault(q,[]).append(r)
with open(dst,'w',newline='') as f:
 outfields=['query','hit_count','top_accession','top_identity','top_qcov','top_evalue','top_bitscore','top_title','near_identity_ge90_qcov80']
 w=csv.DictWriter(f,fieldnames=outfields,delimiter='\t');w.writeheader()
 for q,rs in sorted(by.items()):
  top=max(rs,key=lambda x:float(x['bitscore']))
  near=any(float(x['pident'])>=90 and float(x['qcovs'])>=80 for x in rs)
  w.writerow({'query':q,'hit_count':len(rs),'top_accession':top['saccver'],'top_identity':top['pident'],
              'top_qcov':top['qcovs'],'top_evalue':top['evalue'],'top_bitscore':top['bitscore'],
              'top_title':top['stitle'],'near_identity_ge90_qcov80':str(near).lower()})
PY

# Build one compact query-plus-RefSeq-homolog tree per candidate.
python3 - "$PROT_QUERY" "$OUT/results/refseq_viral_blastp.tsv" "$OUT/trees/top_ids" <<'PY'
import csv,pathlib,sys
query,tsv,prefix=sys.argv[1:]
fields=['qseqid','saccver','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','qcovs','stitle']
by={}
with open(tsv) as f:
 for vals in csv.reader(f,delimiter='\t'):
  if len(vals)<len(fields):continue
  r=dict(zip(fields,vals[:len(fields)]));by.setdefault(r['qseqid'],[]).append(r)
for q,rs in by.items():
 rs.sort(key=lambda x:float(x['bitscore']),reverse=True)
 ids=[]
 for r in rs:
  if r['saccver'] not in ids:ids.append(r['saccver'])
  if len(ids)>=20:break
 pathlib.Path(f'{prefix}_{q}.txt').write_text('\n'.join(ids)+'\n')
PY

for ids in "$OUT"/trees/top_ids_*.txt; do
  [[ -s "$ids" ]] || continue
  q="${ids##*_}"; q="${q%.txt}"
  blastdbcmd -db "$OUT/db/refseq_viral_prot" -entry_batch "$ids" \
    > "$OUT/trees/${q}_homologs.faa" 2> "$OUT/trees/${q}_blastdbcmd.log" || true
  python3 - "$PROT_QUERY" "$q" "$OUT/trees/${q}_query.faa" <<'PY'
from Bio import SeqIO
import sys
src,q,dst=sys.argv[1:]
for r in SeqIO.parse(src,'fasta'):
 if r.id==q:
  SeqIO.write([r],dst,'fasta');break
PY
  cat "$OUT/trees/${q}_query.faa" "$OUT/trees/${q}_homologs.faa" > "$OUT/trees/${q}_combined.faa"
  mafft --auto "$OUT/trees/${q}_combined.faa" > "$OUT/trees/${q}.aln.faa" 2> "$OUT/trees/${q}.mafft.log"
  FastTree -wag "$OUT/trees/${q}.aln.faa" > "$OUT/trees/${q}.tree.nwk" 2> "$OUT/trees/${q}.fasttree.log"
done

{
  date -u +%FT%TZ
  blastp -version
  blastn -version
  mafft --version 2>&1 | head -n2
  FastTree 2>&1 | head -n3 || true
  du -sh "$OUT/db"/*
} > "$OUT/results/PROVENANCE.txt"

# Large downloaded databases are reproducible from the preserved index/file list
# and are omitted from the final evidence artifact. Keep their source hashes.
sha256sum "$OUT"/db/viral.*.gz > "$OUT/results/DOWNLOADED_DATABASE_SHA256.txt"
rm -rf "$OUT/db"
find "$OUT" -type f ! -name SHA256SUMS_COMPLETE.txt -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS_COMPLETE.txt"
