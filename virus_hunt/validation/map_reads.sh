#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?run accession required}"
REFS="${2:?reference FASTA required}"
OUT="${3:?output directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT" "$RUN.fastq"

META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,base_count,read_count,first_public,last_updated,country,collection_date,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
META_FALLBACK_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,base_count,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
if ! curl -fL --retry 8 --retry-delay 5 "$META_URL" -o "$OUT/ena_run_metadata.tsv"; then
  curl -fL --retry 8 --retry-delay 5 "$META_FALLBACK_URL" -o "$OUT/ena_run_metadata.tsv"
fi

python3 - "$OUT/ena_run_metadata.tsv" "$RUN.fastq/download_plan.tsv" <<'PY'
import csv,sys
src,dst=sys.argv[1:]
rows=list(csv.DictReader(open(src),delimiter='\t'))
if not rows:
 raise SystemExit('ENA returned no metadata row')
r=rows[0]
urls=(r.get('fastq_ftp') or '').split(';') if r.get('fastq_ftp') else []
md5s=(r.get('fastq_md5') or '').split(';') if r.get('fastq_md5') else []
bytes_=(r.get('fastq_bytes') or '').split(';') if r.get('fastq_bytes') else []
with open(dst,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['url','md5','bytes'])
 for i,u in enumerate(urls):
  if not u: continue
  if not u.startswith(('http://','https://')): u='https://'+u
  w.writerow([u,md5s[i] if i<len(md5s) else '',bytes_[i] if i<len(bytes_) else ''])
PY

mapfile -t URLS < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f1)
mapfile -t MD5S < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f2)
if [[ ${#URLS[@]} -eq 0 ]]; then
  echo "No ENA FASTQ URLs returned for $RUN" >&2
  exit 2
fi
FILES=()
for i in "${!URLS[@]}"; do
  url="${URLS[$i]}"
  file="$RUN.fastq/$(basename "$url")"
  curl -fL --retry 8 --retry-delay 10 --continue-at - "$url" -o "$file"
  if [[ -n "${MD5S[$i]:-}" ]]; then
    echo "${MD5S[$i]}  $file" | md5sum -c -
  fi
  FILES+=("$file")
done

minimap2 --version > "$OUT/minimap2_version.txt"
samtools --version > "$OUT/samtools_version.txt"
if [[ ${#FILES[@]} -ge 2 ]]; then
  minimap2 -ax sr --secondary=no -t 4 "$REFS" \
    <(pigz -dc "${FILES[0]}") <(pigz -dc "${FILES[1]}") \
    2> "$OUT/minimap2.stderr.txt" \
    | samtools view -b -F 4 - \
    | samtools sort -@ 2 -o "$OUT/mapped.bam" -
else
  minimap2 -ax sr --secondary=no -t 4 "$REFS" \
    <(pigz -dc "${FILES[0]}") \
    2> "$OUT/minimap2.stderr.txt" \
    | samtools view -b -F 4 - \
    | samtools sort -@ 2 -o "$OUT/mapped.bam" -
fi
samtools index "$OUT/mapped.bam"
samtools flagstat "$OUT/mapped.bam" > "$OUT/flagstat.txt"
samtools idxstats "$OUT/mapped.bam" > "$OUT/idxstats.tsv"
samtools coverage "$OUT/mapped.bam" > "$OUT/samtools_coverage.tsv"
python3 "$SCRIPT_DIR/summarize_read_support.py" \
  --run "$RUN" --bam "$OUT/mapped.bam" --refs "$REFS" \
  --metadata "$OUT/ena_run_metadata.tsv" --output "$OUT/read_support.tsv"

rm -rf "$RUN.fastq" "$OUT/mapped.bam" "$OUT/mapped.bam.bai"
find "$OUT" -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS.txt"
