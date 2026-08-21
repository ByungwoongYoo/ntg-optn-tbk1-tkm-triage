#!/usr/bin/env bash
# Current NCBI remote-database audit for immutable Panax A1/A2/B queries.
# A zero-hit result is accepted only when BLAST produced a valid archive,
# recovered every exact query, and reported nonzero statistics. Entrez-
# partitioned modes also carry a same-request positive control so an empty
# candidate result cannot be confused with a silently invalid remote search.
set -Eeuo pipefail

MODE="${1:?search mode required}"
QUERY_ROOT="${2:?query directory required}"
OUT="${3:?output directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
partition_controls=()
mkdir -p "$OUT"
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi

# Once the empty output directory is under our control, every nonzero exit must
# leave a checksummed diagnostic artifact. This includes failures before the
# remote-search loop, such as a missing tool/query, bad environment value,
# malformed control manifest, or an internal status-builder error.
finalization_complete=0
failure_line=unknown
failure_command=explicit_exit
capture_error() {
  failure_line="$1"
  failure_command="$2"
}
finalize_incomplete() {
  local rc="$1"
  trap - ERR EXIT
  set +e
  [[ -f "$OUT/STARTED_UTC.txt" ]] || date -u +%FT%TZ > "$OUT/STARTED_UTC.txt"
  printf '0\n' > "$OUT/SEARCH_SUCCESS.txt"
  [[ -f "$OUT/ATTEMPT_COUNT.txt" ]] || printf '0\n' > "$OUT/ATTEMPT_COUNT.txt"
  [[ -f "$OUT/SUCCESS_ATTEMPT.txt" ]] || printf '0\n' > "$OUT/SUCCESS_ATTEMPT.txt"
  [[ -f "$OUT/SEARCH_BUDGET_SECONDS.txt" ]] || printf '0\n' > "$OUT/SEARCH_BUDGET_SECONDS.txt"
  printf 'preflight_or_internal_failure\n' > "$OUT/TERMINATION_REASON.txt"
  [[ -f "$OUT/HITS.tsv" ]] || : > "$OUT/HITS.tsv"
  [[ -f "$OUT/STDOUT.txt" ]] || : > "$OUT/STDOUT.txt"
  {
    [[ -f "$OUT/STDERR.txt" ]] && cat "$OUT/STDERR.txt"
    printf 'runner_exit_code=%s\nfailure_line=%s\nfailure_command=%s\n' \
      "$rc" "$failure_line" "$failure_command"
  } > "$OUT/STDERR.finalization.tmp"
  mv "$OUT/STDERR.finalization.tmp" "$OUT/STDERR.txt"
  [[ -f "$OUT/REMOTE_ATTEMPTS.tsv" ]] || \
    printf 'attempt\tstart_utc\tend_utc\tbackoff_before_seconds\tattempt_timeout_seconds\tblast_rc\tjson_formatter_rc\ttsv_formatter_rc\tvalidator_rc\tfailure_stage\tfailure_class\tretryable\tresult_archive_bytes\tresult_archive_sha256\n' \
      > "$OUT/REMOTE_ATTEMPTS.tsv"
  python - "$OUT/SEARCH_STATUS.json" "$MODE" "$rc" "$failure_line" "$failure_command" <<'PY'
from datetime import datetime, timezone
from pathlib import Path
import json, sys
path, mode, rc, line, command = sys.argv[1:]
payload = {
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "mode": mode,
    "command_completed_successfully": False,
    "result_archive_valid": False,
    "technical_complete": False,
    "failure_stage": "preflight_or_internal",
    "runner_exit_code": int(rc),
    "failure_line": line,
    "failure_command": command,
    "query_ids": [],
    "validation_control_ids": [],
    "validation_control_results": {},
    "per_query": {},
    "interpretation_boundary": (
        "No biological or sequence-level inference is permitted from this "
        "incomplete diagnostic artifact."
    ),
}
Path(path).write_text(json.dumps(payload, indent=2) + "\n")
PY
  date -u +%FT%TZ > "$OUT/FINISHED_UTC.txt"
  (
    cd "$OUT" || exit 0
    find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | \
      xargs -0 sha256sum > SHA256SUMS.txt
  )
}
on_exit() {
  local rc="$1"
  if (( rc != 0 && finalization_complete == 0 )); then
    finalize_incomplete "$rc"
  fi
}
trap 'capture_error "$LINENO" "$BASH_COMMAND"' ERR
trap 'on_exit "$?"' EXIT

# Standard-task nr/nt coverage is split into explicit viral and indexed-
# complement Entrez partitions. The complement uses an all-record left operand
# before NOT; the unfiltered remote service can emit zero-statistic archives that
# look successful but contain no usable search result.
case "$MODE" in
  protein_viral)
    program=blastp; query="$QUERY_ROOT/panax_candidates_plus_controls_orfs.faa"; database=nr
    extra=(-entrez_query 'txid10239[ORGN]' -seg yes -comp_based_stats 2)
    ;;
  protein_nonviral)
    program=blastp; candidate_query="$QUERY_ROOT/panax_three_partial_orfs.faa"
    partition_controls+=("$SCRIPT_DIR/remote_partition_controls.faa")
    query="$OUT/SEARCH_QUERIES.faa"; database=nr
    extra=(-entrez_query 'all[filter] NOT txid10239[ORGN]' -seg yes -comp_based_stats 2)
    ;;
  protein_tsa)
    program=blastp; query="$QUERY_ROOT/panax_three_partial_orfs.faa"; database=tsa_nr
    extra=(-seg yes -comp_based_stats 2)
    ;;
  protein_environmental)
    program=blastp; query="$QUERY_ROOT/panax_three_partial_orfs.faa"; database=env_nr
    extra=(-seg yes -comp_based_stats 2)
    ;;
  nt_viral)
    program=blastn; query="$QUERY_ROOT/panax_candidates_plus_controls_contigs.fna"; database=nt
    extra=(-task blastn -entrez_query 'txid10239[ORGN]' -dust yes -soft_masking true)
    ;;
  nt_nonviral)
    program=blastn; candidate_query="$QUERY_ROOT/panax_three_contigs.fna"
    partition_controls+=("$SCRIPT_DIR/remote_partition_controls.fna")
    partition_controls+=("$SCRIPT_DIR/remote_nonpanax_control.fna")
    query="$OUT/SEARCH_QUERIES.fna"; database=nt
    extra=(-task blastn -entrez_query 'all[filter] NOT txid10239[ORGN]' -dust yes -soft_masking true)
    ;;
  nt_megablast)
    program=blastn; query="$QUERY_ROOT/panax_candidates_plus_controls_contigs.fna"; database=nt
    extra=(-task megablast -dust yes -soft_masking true)
    ;;
  nt_panax)
    program=blastn; candidate_query="$QUERY_ROOT/panax_three_contigs.fna"
    partition_controls+=("$SCRIPT_DIR/remote_partition_controls.fna")
    query="$OUT/SEARCH_QUERIES.fna"; database=nt
    extra=(-task blastn -entrez_query 'txid44586[ORGN]' -dust yes -soft_masking true)
    ;;
  nt_tsa)
    program=blastn; query="$QUERY_ROOT/panax_three_contigs.fna"; database=tsa_nt
    extra=(-task blastn -dust yes -soft_masking true)
    ;;
  *)
    echo "unknown mode: $MODE" >&2
    exit 2
    ;;
esac

if (( ${#partition_controls[@]} > 0 )); then
  control_manifest="$SCRIPT_DIR/remote_partition_controls.json"
  python - "$control_manifest" "$MODE" "${partition_controls[@]}" <<'PY'
from pathlib import Path
import hashlib,json,sys
manifest_path=Path(sys.argv[1]); mode=sys.argv[2]; fasta_paths=[Path(x) for x in sys.argv[3:]]
manifest=json.loads(manifest_path.read_text())
controls=[
    row for row in manifest.get('controls',[])
    if mode in row.get('required_modes',[])
]
if len(controls) != len(fasta_paths) or {
    row.get('fasta_file') for row in controls
} != {path.name for path in fasta_paths}:
    raise SystemExit(
        f'partition-control file set mismatch for {mode}: '
        f'controls={[row.get("fasta_file") for row in controls]}, '
        f'files={[path.name for path in fasta_paths]}'
    )
rows_by_file={row['fasta_file']:row for row in controls}
if len(rows_by_file) != len(controls):
    raise SystemExit(f'duplicate partition-control FASTA contract for {mode}')
for fasta in fasta_paths:
    row=rows_by_file[fasta.name]
    if hashlib.sha256(fasta.read_bytes()).hexdigest() != row.get('fasta_file_sha256'):
        raise SystemExit(f'control FASTA hash mismatch for {mode}: {fasta.name}')
    records={}; name=None
    for raw in fasta.read_text().splitlines():
        line=raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in records:
                raise SystemExit(f'duplicate partition control: {name}')
            records[name]=''
        elif name is None:
            raise SystemExit('sequence before partition-control header')
        else:
            records[name]+=line.upper()
    if set(records) != {row['control']}:
        raise SystemExit(f'partition-control ID mismatch: {sorted(records)}')
    seq=records[row['control']]
    if len(seq) != int(row['length']):
        raise SystemExit(f'partition-control length mismatch: {row["control"]}')
    if hashlib.sha256(seq.encode()).hexdigest() != row['sequence_sha256']:
        raise SystemExit(f'partition-control sequence hash mismatch: {row["control"]}')
PY
  cat "$candidate_query" "${partition_controls[@]}" > "$query"
  for partition_control in "${partition_controls[@]}"; do
    cp "$partition_control" "$OUT/"
  done
  cp "$control_manifest" "$OUT/REMOTE_PARTITION_CONTROLS.json"
fi

# Emit immutable expected query metadata before making a network request.
python - "$query" "$program" "$MODE" "$OUT/EXPECTED_QUERIES.json" <<'PY'
from pathlib import Path
import hashlib,json,re,sys
p=Path(sys.argv[1]); program=sys.argv[2]; mode=sys.argv[3]; out=Path(sys.argv[4])
records={}; name=None
for raw in p.read_text().splitlines():
    line=raw.strip()
    if not line: continue
    if line.startswith('>'):
        name=line[1:].split()[0]
        if name in records: raise SystemExit(f'duplicate query id: {name}')
        records[name]=''
    elif name is None:
        raise SystemExit(f'sequence before first FASTA header: {p}')
    else:
        records[name]+=line.upper()
expected={'PNX_Picorna_A1','PNX_Picorna_A2','PNX_Picorna_B'}
controls={'PNX_Duplo_A_control','PNX_Duplo_B_control'}
mode_controls={
    'protein_viral': controls,
    'nt_viral': controls,
    'nt_megablast': controls,
    'protein_nonviral': {'PNX_Panax_L2_control'},
    'nt_nonviral': {
        'PNX_Panax_cpDNA_control', 'PNX_NonPanax_mtDNA_control',
    },
    'nt_panax': {'PNX_Panax_cpDNA_control'},
}
validation_controls=mode_controls.get(mode,set())
exact_control_expectations={
    'protein_nonviral': {
        'PNX_Panax_L2_control': {
            'expected_accession': 'YP_009121238.1',
            'min_query_coverage': 99.0,
            'min_identity': 99.0,
        },
    },
    'nt_nonviral': {
        'PNX_Panax_cpDNA_control': {
            'expected_accession': 'NC_026447.1',
            'min_query_coverage': 99.0,
            'min_identity': 99.0,
        },
        'PNX_NonPanax_mtDNA_control': {
            'expected_accession': 'NC_012920.1',
            'min_query_coverage': 99.0,
            'min_identity': 99.0,
        },
    },
    'nt_panax': {
        'PNX_Panax_cpDNA_control': {
            'expected_accession': 'NC_026447.1',
            'min_query_coverage': 99.0,
            'min_identity': 99.0,
        },
    },
}
required=expected|validation_controls
if set(records)!=required:
    raise SystemExit(f'query set mismatch: observed={sorted(records)}, required={sorted(required)}')
allowed=r'[ACGTN]+' if program=='blastn' else r'[ABCDEFGHIKLMNPQRSTVWXYZ]+'
for name,seq in records.items():
    if not seq or not re.fullmatch(allowed,seq):
        raise SystemExit(f'empty or invalid query sequence: {name}')
payload={'query_file':str(p),'query_file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
         'candidate_ids':sorted(expected),'validation_control_ids':sorted(validation_controls),
         'validation_controls':[
             {'id':control,**exact_control_expectations.get(mode,{}).get(control,{})}
             for control in sorted(validation_controls)
         ],
         'queries':[{'id':k,'length':len(v),'sequence_sha256':hashlib.sha256(v.encode()).hexdigest()}
                    for k,v in records.items()]}
out.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(payload,indent=2))
PY

date -u +%FT%TZ > "$OUT/STARTED_UTC.txt"
"$program" -version > "$OUT/TOOL_VERSION.txt" 2>&1
blast_formatter -version >> "$OUT/TOOL_VERSION.txt" 2>&1
sha256sum "$query" > "$OUT/QUERY_SHA256.txt"

# ASN.1 archive output permits lossless conversion to both machine-readable
# BLAST JSON and a compact TSV without repeating the remote search.
cmd=("$program" -remote -query "$query" -db "$database" -evalue 1e-5 \
     -max_target_seqs 100 -max_hsps 1 -outfmt 11 -out "$OUT/RESULTS.asn" "${extra[@]}")
printf '%q ' "${cmd[@]}" > "$OUT/COMMAND.txt"
printf '\n' >> "$OUT/COMMAND.txt"
: > "$OUT/STDOUT.txt"
: > "$OUT/STDERR.txt"

validate_remote_archive() {
  python "$SCRIPT_DIR/validate_panax_remote_archive.py" "$@"
}

success=0
attempts=0
success_attempt=0
# A short immediate retry is insufficient for a shared remote service.  The
# 2026-08-21 audit saw six otherwise unrelated matrix cells return the same
# Blast4 transport exception in a 20-minute window, while sibling modes later
# completed.  Keep the scientific gate fail-closed, but retry the transport on
# a bounded minute-scale schedule that stays well below NCBI's request-rate
# guidance.  The schedule is deterministic so the provenance is auditable.
max_attempts="${PANAX_REMOTE_MAX_ATTEMPTS:-8}"
[[ "$max_attempts" =~ ^[1-8]$ ]] || {
  echo "PANAX_REMOTE_MAX_ATTEMPTS must be an integer from 1 through 8" >&2
  exit 2
}
attempt_timeout_seconds="${PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS:-6300}"
[[ "$attempt_timeout_seconds" =~ ^(0|[1-9][0-9]{0,4})$ ]] || {
  echo "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS must be a canonical integer from 60 through 6300" >&2
  exit 2
}
attempt_timeout_seconds=$((10#$attempt_timeout_seconds))
(( attempt_timeout_seconds >= 60 && attempt_timeout_seconds <= 6300 )) || {
  echo "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS must be an integer from 60 through 6300" >&2
  exit 2
}
# GitHub-hosted jobs have a hard six-hour ceiling.  Stop network work after a
# five-hour internal budget so logs, status, checksums, and the failure artifact
# are always finalized before the outer 360-minute job timeout.
search_budget_seconds="${PANAX_REMOTE_SEARCH_BUDGET_SECONDS:-18000}"
[[ "$search_budget_seconds" =~ ^(0|[1-9][0-9]{0,4})$ ]] || {
  echo "PANAX_REMOTE_SEARCH_BUDGET_SECONDS must be a canonical integer from 300 through 19800" >&2
  exit 2
}
search_budget_seconds=$((10#$search_budget_seconds))
(( search_budget_seconds >= 300 && search_budget_seconds <= 19800 )) || {
  echo "PANAX_REMOTE_SEARCH_BUDGET_SECONDS must be an integer from 300 through 19800" >&2
  exit 2
}
minimum_attempt_seconds=60
backoff_seconds=(0 120 300 600 900 1200 1800 2700)
search_start_epoch="$(date +%s)"
search_deadline_epoch=$((search_start_epoch + search_budget_seconds))
termination_reason=max_attempts_exhausted
printf '%s\n' "$search_budget_seconds" > "$OUT/SEARCH_BUDGET_SECONDS.txt"
printf 'attempt\tstart_utc\tend_utc\tbackoff_before_seconds\tattempt_timeout_seconds\tblast_rc\tjson_formatter_rc\ttsv_formatter_rc\tvalidator_rc\tfailure_stage\tfailure_class\tretryable\tresult_archive_bytes\tresult_archive_sha256\n' \
  > "$OUT/REMOTE_ATTEMPTS.tsv"
for attempt in $(seq 1 "$max_attempts"); do
  attempt_stdout="$OUT/STDOUT.attempt${attempt}.txt"
  attempt_stderr="$OUT/STDERR.attempt${attempt}.txt"
  backoff="${backoff_seconds[$((attempt-1))]}"

  now_epoch="$(date +%s)"
  remaining_seconds=$((search_deadline_epoch - now_epoch))
  if (( remaining_seconds < backoff + minimum_attempt_seconds )); then
    termination_reason=search_budget_exhausted_before_attempt
    break
  fi
  if (( backoff > 0 )); then
    sleep "$backoff"
  fi
  now_epoch="$(date +%s)"
  remaining_seconds=$((search_deadline_epoch - now_epoch))
  current_timeout_seconds="$attempt_timeout_seconds"
  if (( current_timeout_seconds > remaining_seconds )); then
    current_timeout_seconds="$remaining_seconds"
  fi
  if (( current_timeout_seconds < minimum_attempt_seconds )); then
    termination_reason=search_budget_exhausted_before_attempt
    break
  fi

  attempts="$attempt"
  start_utc="$(date -u +%FT%TZ)"
  printf '[%s] attempt=%s backoff_before_seconds=%s timeout_seconds=%s\n' \
    "$start_utc" "$attempt" "$backoff" "$current_timeout_seconds" > "$attempt_stderr"
  : > "$attempt_stdout"
  rm -f "$OUT/RESULTS.asn" "$OUT/RESULTS.json" "$OUT/HITS.tsv"
  blast_rc=0
  json_formatter_rc=-1
  tsv_formatter_rc=-1
  validator_rc=-1
  failure_stage=none
  failure_class=none
  retryable=1

  timeout --signal=TERM --kill-after=60s "${current_timeout_seconds}s" "${cmd[@]}" \
    >> "$attempt_stdout" 2>> "$attempt_stderr" || blast_rc=$?
  if (( blast_rc != 0 )); then
    failure_stage=blast
    if (( blast_rc == 124 || blast_rc == 137 )); then
      failure_class=transient_timeout
    elif grep -Eiq \
      'unknown argument|invalid argument|command line argument error|query is empty|fasta-reader|cannot open|no such file|BLAST query/options error' \
      "$attempt_stderr"; then
      failure_class=deterministic_input_error
      retryable=0
    elif grep -Eiq \
      'Blast4-request|CRPCClientException|connection stream is in bad state|timed out|timeout|temporar|try again|connection (reset|closed|refused)|service unavailable|HTTP[^0-9]*(429|5[0-9][0-9])' \
      "$attempt_stderr"; then
      failure_class=transient_remote_transport
    else
      failure_class=remote_process_error
    fi
  elif [[ ! -s "$OUT/RESULTS.asn" ]]; then
    failure_stage=blast
    failure_class=missing_remote_archive
  else
    json_formatter_rc=0
    blast_formatter -archive "$OUT/RESULTS.asn" -outfmt 15 -out "$OUT/RESULTS.json" \
      >> "$attempt_stdout" 2>> "$attempt_stderr" || json_formatter_rc=$?
    if (( json_formatter_rc != 0 )) || [[ ! -s "$OUT/RESULTS.json" ]]; then
      failure_stage=json_formatter
      failure_class=structural_remote_archive
    else
      tsv_formatter_rc=0
      blast_formatter -archive "$OUT/RESULTS.asn" \
        -outfmt '6 qseqid saccver sallacc sallseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs staxids sscinames stitle qseq sseq' \
        -out "$OUT/HITS.tsv" >> "$attempt_stdout" 2>> "$attempt_stderr" || tsv_formatter_rc=$?
      if (( tsv_formatter_rc != 0 )); then
        failure_stage=tsv_formatter
        failure_class=structural_remote_archive
      else
        validator_rc=0
        validate_remote_archive \
          "$OUT/RESULTS.json" "$OUT/EXPECTED_QUERIES.json" "$MODE" "$OUT/HITS.tsv" \
          >> "$attempt_stdout" 2>> "$attempt_stderr" || validator_rc=$?
        if (( validator_rc != 0 )); then
          failure_stage=validator
          case "$validator_rc" in
            20)
              failure_class=structural_remote_archive
              ;;
            21)
              failure_class=deterministic_control_failure
              retryable=0
              ;;
            *)
              # An unexpected validator exit is a local/code failure, not a
              # remote no-hit result.  Preserve it fail-closed without burning
              # the search budget on identical local failures.
              failure_class=validator_internal_error
              retryable=0
              ;;
          esac
        fi
      fi
    fi
  fi

  archive_bytes=0
  archive_sha256=NA
  [[ -f "$OUT/RESULTS.asn" ]] && archive_bytes=$(stat -c '%s' "$OUT/RESULTS.asn")
  [[ -f "$OUT/RESULTS.asn" ]] && archive_sha256=$(sha256sum "$OUT/RESULTS.asn" | cut -d' ' -f1)
  if (( blast_rc == 0 && json_formatter_rc == 0 && tsv_formatter_rc == 0 && validator_rc == 0 )); then
    success=1
    success_attempt="$attempt"
    retryable=0
    termination_reason=success
  fi
  end_utc="$(date -u +%FT%TZ)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$attempt" "$start_utc" "$end_utc" "$backoff" "$current_timeout_seconds" \
    "$blast_rc" "$json_formatter_rc" "$tsv_formatter_rc" "$validator_rc" \
    "$failure_stage" "$failure_class" "$retryable" "$archive_bytes" "$archive_sha256" \
    >> "$OUT/REMOTE_ATTEMPTS.tsv"
  if (( success == 1 )); then
    break
  fi
  if (( retryable == 0 )); then
    termination_reason=nonretryable_failure
    break
  fi
done
printf '%s\n' "$termination_reason" > "$OUT/TERMINATION_REASON.txt"
for log in "$OUT"/STDOUT.attempt*.txt; do [[ -f "$log" ]] && { printf '===== %s =====\n' "$(basename "$log")"; cat "$log"; }; done > "$OUT/STDOUT.txt"
for log in "$OUT"/STDERR.attempt*.txt; do [[ -f "$log" ]] && { printf '===== %s =====\n' "$(basename "$log")"; cat "$log"; }; done > "$OUT/STDERR.txt"
printf '%s\n' "$success" > "$OUT/SEARCH_SUCCESS.txt"
printf '%s\n' "$attempts" > "$OUT/ATTEMPT_COUNT.txt"
printf '%s\n' "$success_attempt" > "$OUT/SUCCESS_ATTEMPT.txt"
[[ -f "$OUT/HITS.tsv" ]] || : > "$OUT/HITS.tsv"

python - "$MODE" "$query" "$database" "$OUT" <<'PY'
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,math,sys
mode,query,database,out=sys.argv[1:]; out=Path(out); query=Path(query)
expected=json.loads((out/'EXPECTED_QUERIES.json').read_text())
expected_lengths={x['id']:x['length'] for x in expected['queries']}
validation_controls=expected.get('validation_control_ids',[])
validation_control_specs=expected.get('validation_controls',[])
fields=['qseqid','saccver','sallacc','sallseqid','pident','length','qlen','slen','qstart','qend','sstart','send',
        'evalue','bitscore','qcovs','staxids','sscinames','stitle','qseq','sseq']
rows=[]
tsv_validation_errors=[]
for line_number,values in enumerate(
    csv.reader((out/'HITS.tsv').open(errors='replace'),delimiter='\t'), 1
):
    if not values: continue
    if len(values) != len(fields):
        tsv_validation_errors.append(
            f'line {line_number}: expected {len(fields)} fields, observed {len(values)}'
        )
        continue
    row=dict(zip(fields,values))
    try:
        int(row['qlen'])
        int(row['slen'])
        int(row['length'])
        for key in ('pident','evalue','bitscore','qcovs'):
            value=float(row[key])
            if not math.isfinite(value):
                raise ValueError(f'{key} is not finite')
    except (TypeError,ValueError) as exc:
        tsv_validation_errors.append(f'line {line_number}: malformed numeric field: {exc}')
        continue
    rows.append(row)
with (out/'REMOTE_ATTEMPTS.tsv').open(errors='replace') as handle:
    attempt_history=list(csv.DictReader(handle,delimiter='\t'))
success=(out/'SEARCH_SUCCESS.txt').read_text().strip()=='1'
archive_structurally_valid=success or any(
    row.get('blast_rc')=='0'
    and row.get('json_formatter_rc')=='0'
    and row.get('tsv_formatter_rc')=='0'
    and row.get('validator_rc') in {'0','21'}
    for row in attempt_history
)
json_queries={}
json_error=''
if archive_structurally_valid:
    try:
        payload=json.loads((out/'RESULTS.json').read_text())
        reports=payload.get('BlastOutput2',[])
        for item in reports:
            search=item['report']['results']['search']
            title=str(search.get('query_title','')).strip()
            qid=title.split()[0] if title else ''
            qlen=int(search.get('query_len',0))
            if not qid or qid in json_queries:
                raise ValueError(f'missing or duplicated JSON query title: {title!r}')
            json_queries[qid]=qlen
    except Exception as exc:
        json_error=f'{type(exc).__name__}: {exc}'
valid_json=(not json_error and json_queries==expected_lengths)
unexpected=sorted({r['qseqid'] for r in rows}-set(expected_lengths))
bad_qlen=sorted({r['qseqid'] for r in rows
                 if r['qseqid'] in expected_lengths and int(r['qlen'])!=expected_lengths[r['qseqid']]})
success_attempt=(out/'SUCCESS_ATTEMPT.txt').read_text().strip()
stderr_path=out/f'STDERR.attempt{success_attempt}.txt' if success_attempt!='0' else out/'STDERR.txt'
stderr=stderr_path.read_text(errors='replace')
fatal_markers=[x for x in ('Query is Empty','BLAST Database error','Error:','FATAL') if x.lower() in stderr.lower()]
command_valid=bool(
    success and valid_json and not tsv_validation_errors
    and not unexpected and not bad_qlen and not fatal_markers
)
# A failed attempt may still leave a partial or even superficially plausible
# formatted table. Preserve its raw row count for diagnostics, but never publish
# hit/control annotations from rows that did not pass the archive validator.
trusted_rows=rows if command_valid else []
per_query={}
for candidate in expected_lengths:
    hits=[r for r in trusted_rows if r['qseqid']==candidate]
    hits.sort(key=lambda r:float(r['bitscore']),reverse=True)
    distinct={}
    for hit in hits:
        distinct.setdefault(hit['saccver'],hit)
    hits=list(distinct.values())
    top=hits[0] if hits else None
    per_query[candidate]={
        'query_length':expected_lengths[candidate], 'hit_count':len(hits),
        'near_identical_qcov80_pident90_count':sum(float(r['qcovs'])>=80 and float(r['pident'])>=90 for r in hits),
        'near_identical_qcov80_pident95_count':sum(float(r['qcovs'])>=80 and float(r['pident'])>=95 for r in hits),
        'top_hit':None if top is None else {
            k:top[k] for k in fields if k not in {'qseq','sseq','sscinames'}
        },
    }
control_results={}
def exact_accessions(row):
    accessions={row['saccver']}
    accessions.update(token for token in row.get('sallacc','').split(';') if token)
    for identifier in row['sallseqid'].split(';'):
        accessions.update(token for token in identifier.split('|') if token)
    return accessions
for spec in validation_control_specs:
    control=spec['id']
    accession=spec.get('expected_accession')
    matches=[
        row for row in trusted_rows
        if row['qseqid']==control
        and (not accession or accession in exact_accessions(row))
        and float(row['pident'])>=float(spec.get('min_identity',0))
        and float(row['qcovs'])>=float(spec.get('min_query_coverage',0))
    ]
    control_results[control]={
        'expected_accession':accession,
        'min_query_coverage':spec.get('min_query_coverage'),
        'min_identity':spec.get('min_identity'),
        'validated_accessions':sorted(
            {row['saccver'] for row in matches}
            | ({accession} if matches and accession else set())
        )[:10],
        'validated':bool(matches),
    }
status={
    'generated_utc':datetime.now(timezone.utc).isoformat(), 'mode':mode, 'database':database,
    'query_file':str(query), 'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),
    'query_count':len(expected_lengths), 'query_ids':list(expected_lengths),
    'validation_control_ids':validation_controls,
    'validation_control_results':control_results,
    'expected_query_lengths':expected_lengths, 'json_query_lengths':json_queries,
    'command_completed_successfully':success,
    'result_archive_valid':bool(archive_structurally_valid and valid_json),
    'attempt_count':len(attempt_history), 'attempt_history':attempt_history,
    'termination_reason':(out/'TERMINATION_REASON.txt').read_text().strip(),
    'search_budget_seconds':int((out/'SEARCH_BUDGET_SECONDS.txt').read_text().strip()),
    'fatal_stderr_markers':fatal_markers, 'unexpected_result_query_ids':unexpected,
    'result_query_length_mismatches':bad_qlen, 'technical_complete':command_valid,
    'result_row_count':len(trusted_rows),
    'unvalidated_diagnostic_row_count':len(rows) if not command_valid else 0,
    'per_query':per_query,
    'annotation_validation':(
        'subject accession/version, title, and taxonomy ID are JSON/TSV-bound; '
        'sscinames is retained only in raw HITS.tsv and is not asserted'
        if command_valid else
        'no hit or control annotation is asserted because archive validation failed'
    ),
    'interpretation_boundary':'An empty hit table is evidence only when technical_complete is true; no-hit or divergence does not establish a new taxon.',
}
if json_error: status['json_validation_error']=json_error
if tsv_validation_errors: status['tsv_validation_errors']=tsv_validation_errors
(out/'SEARCH_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps(status,indent=2))
PY

date -u +%FT%TZ > "$OUT/FINISHED_UTC.txt"
(
  cd "$OUT"
  find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
)
# Preserve the complete failure artifact; the workflow uploads it with
# `if: always()` and then fails this matrix cell when technical_complete=false.
finalization_complete=1
python - "$OUT/SEARCH_STATUS.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1]))['technical_complete'] else 1)
PY
