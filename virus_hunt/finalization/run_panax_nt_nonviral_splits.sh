#!/usr/bin/env bash
# Three candidate+both-control requests for the NCBI nt nonviral partition.
set -Eeuo pipefail

QUERY_ROOT="${1:?query directory required}"
OUT="${2:?output directory required}"
REMOTE_RUNNER="${3:?remote runner required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CANDIDATE_FASTA="$QUERY_ROOT/panax_three_contigs.fna"
PANAX_CONTROL="$SCRIPT_DIR/remote_partition_controls.fna"
NONPANAX_CONTROL="$SCRIPT_DIR/remote_nonpanax_control.fna"

[[ -d "$OUT" ]] || { echo "output directory missing: $OUT" >&2; exit 2; }
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi
date -u +%FT%TZ > "$OUT/STARTED_UTC.txt"
cp "$PANAX_CONTROL" "$NONPANAX_CONTROL" "$OUT/"
candidate_payload="$OUT/CANDIDATE_QUERIES.fna"
cp "$CANDIDATE_FASTA" "$candidate_payload"
cat "$CANDIDATE_FASTA" "$PANAX_CONTROL" "$NONPANAX_CONTROL" \
  > "$OUT/SEARCH_QUERIES.fna"
mkdir -p "$OUT/SPLITS"

split_budget="${PANAX_NT_NONVIRAL_SPLIT_BUDGET_SECONDS:-6000}"
[[ "$split_budget" =~ ^(0|[1-9][0-9]{2,4})$ ]] || {
  echo "PANAX_NT_NONVIRAL_SPLIT_BUDGET_SECONDS must be a canonical integer" >&2
  exit 2
}
split_budget=$((10#$split_budget))
(( split_budget >= 300 && split_budget <= 6000 )) || {
  echo "PANAX_NT_NONVIRAL_SPLIT_BUDGET_SECONDS must be 300 through 6000" >&2
  exit 2
}
printf '%s\n' "$((split_budget * 3))" > "$OUT/SEARCH_BUDGET_SECONDS.txt"

split_rc=0
for candidate in PNX_Picorna_A1 PNX_Picorna_A2 PNX_Picorna_B; do
  child_out="$OUT/SPLITS/$candidate"
  PANAX_NONVIRAL_SPLIT_CANDIDATE="$candidate" \
  PANAX_REMOTE_SEARCH_BUDGET_SECONDS="$split_budget" \
    bash "$REMOTE_RUNNER" nt_nonviral "$QUERY_ROOT" "$child_out" || split_rc=1
done

aggregate_rc=0
python "$SCRIPT_DIR/aggregate_panax_nonviral_splits.py" \
  --mode nt_nonviral \
  --out "$OUT" \
  --candidate-fasta "$candidate_payload" \
  --control-fasta "$PANAX_CONTROL" \
  --control-fasta "$NONPANAX_CONTROL" \
  --query-prefix "$OUT" || aggregate_rc=$?
date -u +%FT%TZ > "$OUT/FINISHED_UTC.txt"
(
  cd "$OUT"
  find . -type f ! -path './SHA256SUMS.txt' -print0 | sort -z | \
    xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
if (( split_rc != 0 || aggregate_rc != 0 )); then
  exit 1
fi
python - "$OUT/SEARCH_STATUS.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1]))['technical_complete'] else 1)
PY
