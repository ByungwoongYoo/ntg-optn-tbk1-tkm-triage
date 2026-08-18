#!/usr/bin/env bash
set -euo pipefail

RUN="${1:?run accession required}"
REFS="${2:?reference FASTA required}"
OUT="${3:?output directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT" "$RUN.fastq"

META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,base_count,read_count,first_public,last_updated,country,collection_date,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
META_FALLBACK_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,scientific_name,sample_title,library_strategy,library_source,library_selection,library_layout,base_count,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
if ! curl -fL --retry 12 --retry-all-errors --retry-delay 5 --connect-timeout 30 "$META_URL" -o "$OUT/ena_run_metadata.tsv"; then
  curl -fL --retry 12 --retry-all-errors --retry-delay 5 --connect-timeout 30 "$META_FALLBACK_URL" -o "$OUT/ena_run_metadata.tsv"
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
mapfile -t EXPECTED_BYTES < <(tail -n +2 "$RUN.fastq/download_plan.tsv" | cut -f3)
if [[ ${#URLS[@]} -eq 0 ]]; then
  echo "No ENA FASTQ URLs returned for $RUN" >&2
  exit 2
fi

FILES=()
for i in "${!URLS[@]}"; do
  url="${URLS[$i]}"
  base="$(basename "$url")"
  file="$RUN.fastq/$base"
  echo "[$(date -u +%FT%TZ)] downloading $url" | tee -a "$OUT/download.log"

  # ENA occasionally closes long transfers. aria2 resumes partial byte ranges and
  # retries indefinitely within the enclosing GitHub job timeout. A curl fallback
  # is retained for environments where aria2 is unavailable.
  if command -v aria2c >/dev/null 2>&1; then
    aria2c \
      --continue=true \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --file-allocation=none \
      --max-connection-per-server=4 \
      --split=4 \
      --min-split-size=32M \
      --max-tries=0 \
      --retry-wait=5 \
      --connect-timeout=30 \
      --timeout=60 \
      --lowest-speed-limit=1K \
      --console-log-level=notice \
      --summary-interval=30 \
      --dir "$RUN.fastq" \
      --out "$base" \
      "$url" 2>&1 | tee -a "$OUT/download.log"
  else
    curl -fL \
      --retry 60 --retry-all-errors --retry-delay 5 \
      --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
      --continue-at - "$url" -o "$file" 2>&1 | tee -a "$OUT/download.log"
  fi

  if [[ ! -s "$file" ]]; then
    echo "Downloaded file is absent or empty: $file" >&2
    exit 3
  fi
  actual_bytes="$(stat -c '%s' "$file")"
  expected_bytes="${EXPECTED_BYTES[$i]:-}"
  if [[ -n "$expected_bytes" && "$actual_bytes" != "$expected_bytes" ]]; then
    echo "Byte-count mismatch for $file: expected=$expected_bytes actual=$actual_bytes" >&2
    exit 4
  fi
  if [[ -n "${MD5S[$i]:-}" ]]; then
    echo "${MD5S[$i]}  $file" | md5sum -c -
  fi
  FILES+=("$file")
done

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
python3 "$SCRIPT_DIR/summarize_read_support.py" \
  --run "$RUN" --bam "$OUT/mapped.bam" --refs "$REFS" \
  --metadata "$OUT/ena_run_metadata.tsv" --output "$OUT/read_support.tsv"

rm -rf "$RUN.fastq" "$OUT/mapped.bam" "$OUT/mapped.bam.bai"
find "$OUT" -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS.txt"
