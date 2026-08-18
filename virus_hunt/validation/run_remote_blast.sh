#!/usr/bin/env bash
# Current-database audit for four candidate RdRP proteins/contigs.
# Failure is recorded and causes a conservative incomplete verdict; it is never
# converted into evidence of novelty.
set -uo pipefail

PROT="${1:?candidate protein FASTA required}"
NT="${2:?candidate nucleotide FASTA required}"
OUT="${3:?output directory required}"
mkdir -p "$OUT"

FMT='6 qseqid saccver pident length qlen slen qcovs evalue bitscore staxids sscinames stitle'
: > "$OUT/blastp_viral.tsv"
: > "$OUT/blastp_nonviral.tsv"
: > "$OUT/blastp_nr.tsv"
: > "$OUT/blastn_nt.tsv"
: > "$OUT/REMOTE_BLAST_COMMANDS.txt"

run_with_retry() {
  local label="$1" output="$2"; shift 2
  local log="$OUT/${label}.log"
  local ok=0
  : > "$log"
  printf '%q ' "$@" >> "$OUT/REMOTE_BLAST_COMMANDS.txt"
  printf '\n' >> "$OUT/REMOTE_BLAST_COMMANDS.txt"
  for attempt in 1 2; do
    echo "[$(date -u +%FT%TZ)] $label attempt $attempt" | tee -a "$log"
    if timeout 70m "$@" > "$output" 2>> "$log"; then
      # A clean empty file is a valid no-hit result; explicit errors are not.
      if ! grep -Eqi '(^|[[:space:]])(error|failed|exception|cannot|timed out)([[:space:]:]|$)' "$log"; then
        ok=1
        break
      fi
    fi
    sleep $((attempt * 30))
  done
  echo "$ok" > "$OUT/${label}.success"
}

blastp -version > "$OUT/blastp_version.txt" 2>&1 || true
blastn -version > "$OUT/blastn_version.txt" 2>&1 || true

run_with_retry blastp_viral "$OUT/blastp_viral.tsv" \
  blastp -query "$PROT" -db nr -remote \
  -entrez_query 'txid10239[ORGN]' -evalue 1e-5 -max_target_seqs 30 \
  -seg yes -comp_based_stats 2 -outfmt "$FMT"

# Explicit nonviral search to identify a stronger cellular or endogenous match.
run_with_retry blastp_nonviral "$OUT/blastp_nonviral.tsv" \
  blastp -query "$PROT" -db nr -remote \
  -entrez_query 'NOT txid10239[ORGN]' -evalue 1e-5 -max_target_seqs 30 \
  -seg yes -comp_based_stats 2 -outfmt "$FMT"

# Near-identical nucleotide audit.  No hit is not, by itself, proof of novelty.
run_with_retry blastn_nt "$OUT/blastn_nt.tsv" \
  blastn -query "$NT" -db nt -remote -task megablast \
  -evalue 1e-10 -max_target_seqs 30 -dust yes -outfmt "$FMT"

# The combined protein table is reproducibly derived from the two disjoint
# organism searches; its success requires both component searches.
cat "$OUT/blastp_viral.tsv" "$OUT/blastp_nonviral.tsv" > "$OUT/blastp_nr.tsv"
if [[ "$(cat "$OUT/blastp_viral.success" 2>/dev/null)" == 1 && \
      "$(cat "$OUT/blastp_nonviral.success" 2>/dev/null)" == 1 ]]; then
  echo 1 > "$OUT/blastp_nr.success"
else
  echo 0 > "$OUT/blastp_nr.success"
fi

python3 - "$OUT" <<'PY'
import json, pathlib, sys
from datetime import datetime, timezone
out=pathlib.Path(sys.argv[1])
labels=['blastp_viral','blastp_nonviral','blastp_nr','blastn_nt']
status={
 'generated_utc':datetime.now(timezone.utc).isoformat(),
 'database_note':'NCBI remote BLAST databases queried at workflow runtime; empty output is not itself evidence of novelty.',
 'blastp_nr_derivation':'concatenation of disjoint viral and nonviral NR searches',
 'searches':{}
}
for label in labels:
    success=(out/f'{label}.success').read_text().strip()=='1' if (out/f'{label}.success').exists() else False
    result=out/f'{label}.tsv'
    status['searches'][label]={
      'success':success,
      'result_bytes':result.stat().st_size if result.exists() else 0,
      'result_lines':sum(1 for _ in result.open(errors='ignore')) if result.exists() else 0,
      'log':f'{label}.log' if (out/f'{label}.log').exists() else 'derived'
    }
status['all_required_searches_succeeded']=all(status['searches'][x]['success'] for x in labels)
(out/'REMOTE_BLAST_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps(status,indent=2))
PY

exit 0
