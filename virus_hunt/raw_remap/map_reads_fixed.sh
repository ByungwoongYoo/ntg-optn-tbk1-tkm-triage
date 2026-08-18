#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?run accession required}"
REFS="${2:?reference FASTA required}"
OUT="${3:?output directory required}"
SUMMARIZER="${4:?summarize_read_support.py required}"
mkdir -p "$OUT" "$RUN.fastq"

META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,base_count,read_count,first_public,last_updated,country,collection_date,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
META_FALLBACK_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,base_count,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
if ! curl -fL --retry 12 --retry-all-errors --retry-delay 5 --connect-timeout 30 "$META_URL" -o "$OUT/ena_run_metadata.tsv"; then
  curl -fL --retry 12 --retry-all-errors --retry-delay 5 --connect-timeout 30 "$META_FALLBACK_URL" -o "$OUT/ena_run_metadata.tsv"
fi

# ENA TSV responses may carry CRLF endings. Every scalar is stripped before
# byte/MD5 validation so identical numeric values do not compare unequal.
python3 - "$OUT/ena_run_metadata.tsv" "$RUN.fastq/download_plan.tsv" <<'PY'
import csv,sys
src,dst=sys.argv[1:]
rows=list(csv.DictReader(open(src,newline=''),delimiter='\t'))
if not rows:
 raise SystemExit('ENA returned no metadata row')
r={k:(v or '').strip() for k,v in rows[0].items()}
urls=[x.strip() for x in r.get('fastq_ftp','').split(';') if x.strip()]
md5s=[x.strip() for x in r.get('fastq_md5','').split(';')]
bytes_=[x.strip() for x in r.get('fastq_bytes','').split(';')]
with open(dst,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['url','md5','bytes'])
 for i,u in enumerate(urls):
  if not u.startswith(('http://','https://')): u='https://'+u
  w.writerow([u,md5s[i] if i<len(md5s) else '',bytes_[i] if i<len(bytes_) else ''])
PY

mapfile -t URLS < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f1 | tr -d '\r')
mapfile -t MD5S < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f2 | tr -d '\r')
mapfile -t EXPECTED_BYTES < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f3 | tr -d '\r')
if [[ ${#URLS[@]} -eq 0 ]]; then
  echo "No ENA FASTQ URLs returned for $RUN" >&2
  exit 2
fi

FILES=()
for i in "${!URLS[@]}"; do
  url="${URLS[$i]}"; base="$(basename "$url")"; file="$RUN.fastq/$base"
  echo "[$(date -u +%FT%TZ)] downloading $url" | tee -a "$OUT/download.log"
  aria2c \
    --continue=true --allow-overwrite=true --auto-file-renaming=false \
    --file-allocation=none --max-connection-per-server=4 --split=4 \
    --min-split-size=32M --max-tries=0 --retry-wait=5 \
    --connect-timeout=30 --timeout=60 --lowest-speed-limit=1K \
    --console-log-level=notice --summary-interval=30 \
    --dir "$RUN.fastq" --out "$base" "$url" 2>&1 | tee -a "$OUT/download.log"

  [[ -s "$file" ]] || { echo "Downloaded file absent or empty: $file" >&2; exit 3; }
  actual_bytes="$(stat -c '%s' "$file" | tr -d '[:space:]')"
  expected_bytes="$(printf '%s' "${EXPECTED_BYTES[$i]:-}" | tr -d '[:space:]')"
  if [[ -n "$expected_bytes" && "$actual_bytes" != "$expected_bytes" ]]; then
    echo "Byte-count mismatch for $file: expected=$expected_bytes actual=$actual_bytes" >&2
    exit 4
  fi
  expected_md5="$(printf '%s' "${MD5S[$i]:-}" | tr -d '[:space:]')"
  if [[ -n "$expected_md5" ]]; then
    echo "$expected_md5  $file" | md5sum -c -
  fi
  FILES+=("$file")
done

printf '%s\n' "${FILES[@]}" > "$OUT/downloaded_fastq_files.txt"
df -h . > "$OUT/disk_before_mapping.txt"
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
python3 "$SUMMARIZER" \
  --run "$RUN" --bam "$OUT/mapped.bam" --refs "$REFS" \
  --metadata "$OUT/ena_run_metadata.tsv" --output "$OUT/read_support.tsv"

python3 - "$OUT/read_support.tsv" "$OUT/RAW_MAPPING_STATUS.json" "$RUN" <<'PY'
import csv,json,sys
from datetime import datetime,timezone
src,out,run=sys.argv[1:]
rows=list(csv.DictReader(open(src),delimiter='\t'))
status={'run':run,'success':True,'exit_code':0,'lineages_reported':len(rows),
        'generated_utc':datetime.now(timezone.utc).isoformat()}
open(out,'w').write(json.dumps(status,indent=2)+'\n')
PY

rm -rf "$RUN.fastq" "$OUT/mapped.bam" "$OUT/mapped.bam.bai"
find "$OUT" -type f ! -name SHA256SUMS_COMPLETE.txt -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS_COMPLETE.txt"
