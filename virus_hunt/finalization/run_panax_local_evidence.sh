#!/usr/bin/env bash
# Reproducible local evidence gate using the current official RefSeq viral,
# Pfam-A, and UniVec_Core releases available at run time.
set -Eeuo pipefail

QUERY_ROOT="${1:?query directory required}"
CURATED_MANIFEST="${2:?curated reference manifest required}"
OUT="${3:?output directory required}"
CACHE="${4:?cache/work directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CURRENT_HIT_FASTA="$SCRIPT_DIR/current_nr_top_hit_proteins.faa"
CURRENT_HIT_MANIFEST="$SCRIPT_DIR/current_nr_top_hit_proteins.tsv"
CURRENT_HIT_PANEL_SHA256="$SCRIPT_DIR/current_nr_top_hit_proteins.sha256"
mkdir -p "$OUT" "$CACHE"
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi

for command in curl python md5sum sha256sum gzip makeblastdb blastp blastn \
               dustmasker segmasker hmmpress hmmscan; do
  command -v "$command" >/dev/null || { echo "missing required command: $command" >&2; exit 2; }
done
for fixed_input in "$CURRENT_HIT_FASTA" "$CURRENT_HIT_MANIFEST" "$CURRENT_HIT_PANEL_SHA256"; do
  [[ -s "$fixed_input" ]] || { echo "missing fixed current-hit reference input: $fixed_input" >&2; exit 2; }
done
(cd "$SCRIPT_DIR" && sha256sum -c "$(basename "$CURRENT_HIT_PANEL_SHA256")")

download() {
  local url="$1" destination="$2"
  local tmp="${destination}.part"
  rm -f "$tmp"
  curl --fail --location --silent --show-error --retry 5 --retry-all-errors \
    --connect-timeout 30 --max-time 3600 --output "$tmp" "$url"
  [[ -s "$tmp" ]] || { echo "empty download: $url" >&2; return 1; }
  mv "$tmp" "$destination"
}

date -u +%FT%TZ > "$OUT/STARTED_UTC.txt"
{
  printf 'curl\t'; curl --version | sed -n '1p'
  printf 'python\t'; python --version 2>&1
  printf 'blastp\t'; blastp -version | sed -n '1p'
  printf 'dustmasker\t'; dustmasker -version | sed -n '1p'
  printf 'segmasker\t'; segmasker -version | sed -n '1p'
  printf 'hmmscan\t'; hmmscan -h | sed -n '2p'
} > "$OUT/TOOL_VERSIONS.txt"

REFSEQ_CATALOG_ROOT='https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog'
REFSEQ_RELEASE_ROOT='https://ftp.ncbi.nlm.nih.gov/refseq/release/viral'
PFAM_ROOT='https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release'
UNIVEC_URL='https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'

# Resolve the highest current RefSeq release from the official catalog index;
# do not silently reuse a historical hard-coded mirror.
download "$REFSEQ_CATALOG_ROOT/" "$CACHE/refseq_catalog_index.html"
release="$(python - "$CACHE/refseq_catalog_index.html" <<'PY'
from pathlib import Path
import re,sys
versions={int(x) for x in re.findall(r'release(\d+)\.files\.installed',Path(sys.argv[1]).read_text(errors='replace'))}
if not versions: raise SystemExit('no RefSeq release catalog found in official index')
print(max(versions))
PY
)"
catalog="$CACHE/release${release}.files.installed"
download "$REFSEQ_CATALOG_ROOT/release${release}.files.installed" "$catalog"
catalog_sha="$(sha256sum "$catalog" | awk '{print $1}')"
cp "$catalog" "$OUT/release${release}.files.installed"
printf '%s  %s\n' "$catalog_sha" "release${release}.files.installed" > "$OUT/REFSEQ_RELEASE_CATALOG_SHA256.txt"

python - "$catalog" "$release" "$REFSEQ_RELEASE_ROOT" "$OUT/REFSEQ_SELECTED_FILES.tsv" <<'PY'
from pathlib import Path
import re,sys
catalog,release,root,out=sys.argv[1:]
selected=[]
for line in Path(catalog).read_text().splitlines():
    fields=line.split('\t')
    if len(fields)!=2: continue
    md5,name=fields
    kind=None
    if re.fullmatch(r'viral\.\d+\.protein\.faa\.gz',name): kind='protein'
    if re.fullmatch(r'viral\.\d+\.\d+\.genomic\.fna\.gz',name): kind='genomic'
    if kind: selected.append((kind,md5,name,f'{root}/{name}'))
if not any(x[0]=='protein' for x in selected) or not any(x[0]=='genomic' for x in selected):
    raise SystemExit(f'incomplete viral file selection: {selected}')
with open(out,'w') as f:
    f.write('refseq_release\ttype\tofficial_md5\tfilename\turl\n')
    for kind,md5,name,url in selected: f.write(f'{release}\t{kind}\t{md5}\t{name}\t{url}\n')
PY

mkdir -p "$CACHE/refseq_gz"
printf 'filename\tbytes\tofficial_md5\tobserved_md5\tobserved_sha256\n' > "$OUT/REFSEQ_DOWNLOADED_FILES.tsv"
while IFS=$'\t' read -r selected_release kind expected_md5 filename url; do
  [[ "$selected_release" == "refseq_release" ]] && continue
  destination="$CACHE/refseq_gz/$filename"
  download "$url" "$destination"
  echo "$expected_md5  $destination" | md5sum -c -
  gzip -t "$destination"
  printf '%s\t%s\t%s\t%s\t%s\n' "$filename" "$(stat -c '%s' "$destination")" "$expected_md5" \
    "$(md5sum "$destination" | awk '{print $1}')" "$(sha256sum "$destination" | awk '{print $1}')" \
    >> "$OUT/REFSEQ_DOWNLOADED_FILES.tsv"
done < "$OUT/REFSEQ_SELECTED_FILES.tsv"

mapfile -t protein_gz < <(awk -F '\t' 'NR>1 && $2=="protein" {print $4}' "$OUT/REFSEQ_SELECTED_FILES.tsv" | sort)
mapfile -t genomic_gz < <(awk -F '\t' 'NR>1 && $2=="genomic" {print $4}' "$OUT/REFSEQ_SELECTED_FILES.tsv" | sort)
(( ${#protein_gz[@]} > 0 && ${#genomic_gz[@]} > 0 ))
protein_paths=(); genomic_paths=()
for filename in "${protein_gz[@]}"; do protein_paths+=("$CACHE/refseq_gz/$filename"); done
for filename in "${genomic_gz[@]}"; do genomic_paths+=("$CACHE/refseq_gz/$filename"); done
gzip -cd "${protein_paths[@]}" > "$CACHE/refseq_viral.protein.faa"
gzip -cd "${genomic_paths[@]}" > "$CACHE/refseq_viral.genomic.fna"

python - "$CACHE/refseq_viral.protein.faa" "$CACHE/refseq_viral.genomic.fna" "$OUT/REFSEQ_FASTA_STATS.json" <<'PY'
from pathlib import Path
import hashlib,json,sys
out={}
for label,name in zip(('protein','genomic'),sys.argv[1:3]):
    p=Path(name); records=0; digest=hashlib.sha256()
    with p.open('rb') as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b'>'): records += 1
    if records<1000 or p.stat().st_size<1_000_000:
        raise SystemExit(f'implausibly small current RefSeq viral {label} FASTA: records={records}, bytes={p.stat().st_size}')
    out[label]={'records':records,'bytes':p.stat().st_size,'sha256':digest.hexdigest()}
Path(sys.argv[3]).write_text(json.dumps(out,indent=2)+'\n')
PY

makeblastdb -in "$CACHE/refseq_viral.protein.faa" -dbtype prot -parse_seqids \
  -out "$CACHE/refseq_viral_prot" > "$OUT/MAKEBLASTDB_PROTEIN.log" 2>&1
makeblastdb -in "$CACHE/refseq_viral.genomic.fna" -dbtype nucl -parse_seqids \
  -out "$CACHE/refseq_viral_nt" > "$OUT/MAKEBLASTDB_NUCLEOTIDE.log" 2>&1

blast_fmt='6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore qcovs staxids sscinames stitle sseq'
blastp -query "$QUERY_ROOT/panax_candidates_plus_controls_orfs.faa" -db "$CACHE/refseq_viral_prot" \
  -evalue 1e-5 -max_target_seqs 100 -seg yes -comp_based_stats 2 -outfmt "$blast_fmt" \
  -out "$OUT/LOCAL_REFSEQ_VIRAL_BLASTP.tsv"
blastn -query "$QUERY_ROOT/panax_candidates_plus_controls_contigs.fna" -db "$CACHE/refseq_viral_nt" \
  -task blastn -evalue 1e-5 -max_target_seqs 100 -dust yes -soft_masking true -outfmt "$blast_fmt" \
  -out "$OUT/LOCAL_REFSEQ_VIRAL_BLASTN.tsv"

mkdir -p "$CACHE/references"
python "$SCRIPT_DIR/audit_panax_local_evidence.py" extract-references \
  --manifest "$CURATED_MANIFEST" \
  --current-panel-manifest "$CURRENT_HIT_MANIFEST" \
  --protein-fasta "$CACHE/refseq_viral.protein.faa" \
  --protein-fasta "$CURRENT_HIT_FASTA" \
  --out "$CACHE/references"
cp "$CACHE/references/CURATED_REFERENCE_FULL.faa" "$OUT/CURATED_REFERENCE_FULL.faa"
cp "$CACHE/references/CURATED_REFERENCE_PROVENANCE.json" "$OUT/CURATED_REFERENCE_PROVENANCE.json"
cp "$CACHE/references/CURRENT_NR_REFERENCE_CONTRACT.tsv" "$OUT/CURRENT_NR_REFERENCE_CONTRACT.tsv"
cp "$CACHE/references/CURRENT_NR_REFERENCE_CONTRACT.json" "$OUT/CURRENT_NR_REFERENCE_CONTRACT.json"
# Preserve the original basenames so the copied SHA-256 manifest remains
# directly checkable inside the evidence artifact.
cp "$CURRENT_HIT_FASTA" "$OUT/$(basename "$CURRENT_HIT_FASTA")"
cp "$CURRENT_HIT_MANIFEST" "$OUT/$(basename "$CURRENT_HIT_MANIFEST")"
cp "$CURRENT_HIT_PANEL_SHA256" "$OUT/$(basename "$CURRENT_HIT_PANEL_SHA256")"

# UniVec has no checksum sidecar in its current public directory, so retain the
# HTTPS headers and an observed SHA-256 rather than inventing an official hash.
curl --fail --location --silent --show-error --retry 5 --retry-all-errors \
  --connect-timeout 30 --max-time 600 --dump-header "$OUT/UNIVEC_HTTP_HEADERS.txt" \
  --output "$CACHE/UniVec_Core.part" "$UNIVEC_URL"
[[ -s "$CACHE/UniVec_Core.part" ]]
mv "$CACHE/UniVec_Core.part" "$CACHE/UniVec_Core"
cp "$CACHE/UniVec_Core" "$OUT/UniVec_Core"
printf '%s  UniVec_Core\n' "$(sha256sum "$CACHE/UniVec_Core" | awk '{print $1}')" > "$OUT/UNIVEC_SHA256.txt"
makeblastdb -in "$CACHE/UniVec_Core" -dbtype nucl -parse_seqids -out "$CACHE/univec_core" \
  > "$OUT/MAKEBLASTDB_UNIVEC.log" 2>&1
blastn -task blastn -query "$QUERY_ROOT/panax_three_contigs.fna" -db "$CACHE/univec_core" \
  -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true \
  -evalue 700 -searchsp 1750000000000 -max_target_seqs 1000 \
  -outfmt '6 qseqid saccver pident length qlen qstart qend sstart send evalue score bitscore qseq sseq' \
  -out "$OUT/UNIVEC_RAW.tsv"

dustmasker -in "$QUERY_ROOT/panax_three_contigs.fna" -outfmt fasta -out "$OUT/DUST_MASKED.fna"
segmasker -in "$QUERY_ROOT/panax_three_partial_orfs.faa" -outfmt fasta -out "$OUT/SEG_MASKED.faa"

download "$PFAM_ROOT/md5_checksums" "$OUT/PFAM_MD5_CHECKSUMS.txt"
download "$PFAM_ROOT/Pfam-A.hmm.gz" "$CACHE/Pfam-A.hmm.gz"
pfam_md5="$(python - "$OUT/PFAM_MD5_CHECKSUMS.txt" <<'PY'
from pathlib import Path
import os,re,sys
matches=[]
for line in Path(sys.argv[1]).read_text(errors='replace').splitlines():
    fields=line.replace('*',' ').split()
    if len(fields)>=2 and os.path.basename(fields[-1])=='Pfam-A.hmm.gz' and re.fullmatch(r'[0-9a-fA-F]{32}',fields[0]):
        matches.append(fields[0].lower())
if len(set(matches))!=1: raise SystemExit(f'could not resolve one official Pfam-A.hmm.gz MD5: {matches}')
print(matches[0])
PY
)"
echo "$pfam_md5  $CACHE/Pfam-A.hmm.gz" | md5sum -c -
gzip -t "$CACHE/Pfam-A.hmm.gz"
gzip -cd "$CACHE/Pfam-A.hmm.gz" > "$CACHE/Pfam-A.hmm"
hmmpress -f "$CACHE/Pfam-A.hmm" > "$OUT/HMMPRESS.log" 2>&1
cp "$QUERY_ROOT/panax_three_partial_orfs.faa" "$CACHE/CANDIDATES_AND_REFERENCES.faa"
sed -n '/^>/,$p' "$CACHE/references/CURATED_REFERENCE_FULL.faa" >> "$CACHE/CANDIDATES_AND_REFERENCES.faa"
hmmscan --cpu 2 --cut_ga --noali --domtblout "$OUT/PFAM_DOMAINS.domtblout" \
  "$CACHE/Pfam-A.hmm" "$CACHE/CANDIDATES_AND_REFERENCES.faa" \
  > "$OUT/PFAM_HMMSCAN.txt" 2> "$OUT/PFAM_HMMSCAN.stderr.txt"

printf 'database\trelease_or_snapshot\turl\tofficial_integrity\tobserved_sha256\n' > "$OUT/DATABASE_PROVENANCE.tsv"
printf 'RefSeq_viral\t%s\t%s/release%s.files.installed\tofficial per-file MD5 in catalog\t%s\n' \
  "$release" "$REFSEQ_CATALOG_ROOT" "$release" "$catalog_sha" >> "$OUT/DATABASE_PROVENANCE.tsv"
printf 'Pfam-A\tcurrent_at_run_time\t%s/Pfam-A.hmm.gz\tofficial MD5 %s\t%s\n' \
  "$PFAM_ROOT" "$pfam_md5" "$(sha256sum "$CACHE/Pfam-A.hmm.gz" | awk '{print $1}')" >> "$OUT/DATABASE_PROVENANCE.tsv"
printf 'UniVec_Core\tcurrent_at_run_time\t%s\tno checksum sidecar; HTTPS headers retained\t%s\n' \
  "$UNIVEC_URL" "$(sha256sum "$CACHE/UniVec_Core" | awk '{print $1}')" >> "$OUT/DATABASE_PROVENANCE.tsv"

date -u +%FT%TZ > "$OUT/FINISHED_UTC.txt"
python "$SCRIPT_DIR/audit_panax_local_evidence.py" finalize \
  --query-root "$QUERY_ROOT" \
  --reference-full "$CACHE/references/CURATED_REFERENCE_FULL.faa" \
  --reference-provenance "$CACHE/references/CURATED_REFERENCE_PROVENANCE.tsv" \
  --current-panel-contract "$CACHE/references/CURRENT_NR_REFERENCE_CONTRACT.tsv" \
  --domtbl "$OUT/PFAM_DOMAINS.domtblout" \
  --refseq-blastp "$OUT/LOCAL_REFSEQ_VIRAL_BLASTP.tsv" \
  --refseq-blastn "$OUT/LOCAL_REFSEQ_VIRAL_BLASTN.tsv" \
  --univec "$OUT/UNIVEC_RAW.tsv" --dust-masked "$OUT/DUST_MASKED.fna" \
  --seg-masked "$OUT/SEG_MASKED.faa" --database-provenance "$OUT/DATABASE_PROVENANCE.tsv" \
  --out "$OUT"
