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
mkdir -p "$OUT"
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi

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
    partition_control="$SCRIPT_DIR/remote_partition_controls.faa"
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
    partition_control="$SCRIPT_DIR/remote_partition_controls.fna"
    query="$OUT/SEARCH_QUERIES.fna"; database=nt
    extra=(-task blastn -entrez_query 'all[filter] NOT txid10239[ORGN]' -dust yes -soft_masking true)
    ;;
  nt_megablast)
    program=blastn; query="$QUERY_ROOT/panax_candidates_plus_controls_contigs.fna"; database=nt
    extra=(-task megablast -dust yes -soft_masking true)
    ;;
  nt_panax)
    program=blastn; candidate_query="$QUERY_ROOT/panax_three_contigs.fna"
    partition_control="$SCRIPT_DIR/remote_partition_controls.fna"
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

if [[ -n "${partition_control:-}" ]]; then
  control_manifest="$SCRIPT_DIR/remote_partition_controls.json"
  python - "$control_manifest" "$partition_control" "$MODE" <<'PY'
from pathlib import Path
import hashlib,json,sys
manifest_path=Path(sys.argv[1]); fasta_path=Path(sys.argv[2]); mode=sys.argv[3]
manifest=json.loads(manifest_path.read_text())
controls=[
    row for row in manifest.get('controls',[])
    if mode in row.get('required_modes',[])
]
if len(controls) != 1:
    raise SystemExit(f'expected one partition control for {mode}, found {len(controls)}')
row=controls[0]
fasta=Path(fasta_path)
if fasta.name != row.get('fasta_file'):
    raise SystemExit(f'control FASTA mismatch for {mode}: {fasta.name}')
if hashlib.sha256(fasta.read_bytes()).hexdigest() != row.get('fasta_file_sha256'):
    raise SystemExit(f'control FASTA hash mismatch for {mode}')
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
  cat "$candidate_query" "$partition_control" > "$query"
  cp "$partition_control" "$OUT/"
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
    'nt_nonviral': {'PNX_Panax_cpDNA_control'},
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
  local result_json="$1"
  local expected_json="$2"
  local mode="$3"
  local hits_tsv="$4"
  python - "$result_json" "$expected_json" "$mode" "$hits_tsv" <<'PY'
from pathlib import Path
import json,math,sys
result_path=Path(sys.argv[1]); expected_path=Path(sys.argv[2]); mode=sys.argv[3]
hits_path=Path(sys.argv[4])
payload=json.loads(result_path.read_text())
expected_payload=json.loads(expected_path.read_text())
expected_lengths={x['id']:int(x['length']) for x in expected_payload['queries']}
control_specs=expected_payload.get('validation_controls',[])
reports=payload.get('BlastOutput2',[])
errors=[]; observed={}; hit_counts={}
if not isinstance(reports,list) or not reports:
    errors.append('missing BlastOutput2 reports')
else:
    for index,item in enumerate(reports,1):
        try:
            search=item['report']['results']['search']
        except Exception as exc:
            errors.append(f'report {index} has no search payload: {exc}')
            continue
        title=str(search.get('query_title','')).strip()
        qid=title.split()[0] if title else ''
        try:
            qlen=int(search.get('query_len',0))
        except Exception:
            qlen=0
        if not qid:
            errors.append(f'report {index} has no query title')
            continue
        if qid in observed:
            errors.append(f'duplicate query report: {qid}')
            continue
        observed[qid]=qlen
        stat=search.get('stat') or {}
        for key in ('db_num','db_len'):
            try:
                value=int(stat.get(key,0))
            except Exception:
                value=0
            if value <= 0:
                errors.append(f'{qid} has invalid database statistic: {key}={stat.get(key)!r}')
        bad_stats=[]
        for key in ('kappa','lambda','entropy'):
            try:
                value=float(stat.get(key,0))
            except Exception:
                value=0.0
            if not math.isfinite(value) or value <= 0:
                bad_stats.append(f'{key}={stat.get(key)!r}')
        if bad_stats:
            errors.append(f'{qid} has invalid result statistics: {", ".join(bad_stats)}')
        hits=search.get('hits',[])
        if not isinstance(hits,list):
            errors.append(f'{qid} has a malformed hit list')
            hits=[]
        if not hits and search.get('message') != 'No hits found':
            errors.append(f'{qid} has no hits without the expected completion message')
        hit_counts[qid]=len(hits)
if observed != expected_lengths:
    errors.append(f'query reports mismatch: observed={observed}, expected={expected_lengths}')
hit_rows=[]
hit_fields=['qseqid','saccver','sallacc','sallseqid','pident','length','qlen','slen',
            'qstart','qend','sstart','send','evalue','bitscore','qcovs','staxids',
            'sscinames','stitle','sseq']
for raw in hits_path.read_text(errors='replace').splitlines():
    values=raw.split('\t')
    if len(values) != len(hit_fields):
        errors.append(f'malformed BLAST result row with {len(values)} fields')
        continue
    row=dict(zip(hit_fields,values))
    try:
        row['pident']=float(row['pident'])
        row['qcovs']=float(row['qcovs'])
    except ValueError:
        errors.append('malformed numeric value in BLAST result row')
        continue
    hit_rows.append(row)
def exact_accessions(row):
    accessions={row['saccver']}
    for identifier in row['sallseqid'].split(';'):
        accessions.update(token for token in identifier.split('|') if token)
    return accessions
for spec in control_specs:
    control=spec['id']
    if hit_counts.get(control,0) < 1:
        errors.append(f'positive control has no hit: {control}')
        continue
    accession=spec.get('expected_accession')
    if accession:
        matches=[
            row for row in hit_rows
            if row['qseqid']==control
            and accession in exact_accessions(row)
            and row['pident']>=float(spec['min_identity'])
            and row['qcovs']>=float(spec['min_query_coverage'])
        ]
        if not matches:
            errors.append(
                f'positive control did not recover a near-exact match: {control}'
            )
if errors:
    print('remote archive validation failed:', file=sys.stderr)
    for error in errors:
        print(f'- {error}', file=sys.stderr)
    raise SystemExit(1)
print(f'validated {len(observed)} query reports with nonzero statistics')
PY
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
[[ "$attempt_timeout_seconds" =~ ^[0-9]+$ ]] && \
  (( attempt_timeout_seconds >= 60 && attempt_timeout_seconds <= 6300 )) || {
  echo "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS must be an integer from 60 through 6300" >&2
  exit 2
}
# GitHub-hosted jobs have a hard six-hour ceiling.  Stop network work after a
# five-hour internal budget so logs, status, checksums, and the failure artifact
# are always finalized before the outer 360-minute job timeout.
search_budget_seconds="${PANAX_REMOTE_SEARCH_BUDGET_SECONDS:-18000}"
[[ "$search_budget_seconds" =~ ^[0-9]+$ ]] && \
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
        -outfmt '6 qseqid saccver sallacc sallseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs staxids sscinames stitle sseq' \
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
          if grep -Eiq 'positive control' "$attempt_stderr"; then
            failure_class=deterministic_control_failure
            retryable=0
          else
            failure_class=structural_remote_archive
          fi
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
import csv,hashlib,json,sys
mode,query,database,out=sys.argv[1:]; out=Path(out); query=Path(query)
expected=json.loads((out/'EXPECTED_QUERIES.json').read_text())
expected_lengths={x['id']:x['length'] for x in expected['queries']}
validation_controls=expected.get('validation_control_ids',[])
validation_control_specs=expected.get('validation_controls',[])
fields=['qseqid','saccver','sallacc','sallseqid','pident','length','qlen','slen','qstart','qend','sstart','send',
        'evalue','bitscore','qcovs','staxids','sscinames','stitle','sseq']
rows=[]
for values in csv.reader((out/'HITS.tsv').open(errors='replace'),delimiter='\t'):
    if not values: continue
    if len(values) != len(fields):
        raise SystemExit(f'malformed BLAST result row with {len(values)} fields')
    rows.append(dict(zip(fields,values)))
with (out/'REMOTE_ATTEMPTS.tsv').open(errors='replace') as handle:
    attempt_history=list(csv.DictReader(handle,delimiter='\t'))
success=(out/'SEARCH_SUCCESS.txt').read_text().strip()=='1'
json_queries={}
json_error=''
if success:
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
command_valid=bool(success and valid_json and not unexpected and not bad_qlen and not fatal_markers)
per_query={}
for candidate in expected_lengths:
    hits=[r for r in rows if r['qseqid']==candidate]
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
        'top_hit':None if top is None else {k:top[k] for k in fields if k!='sseq'},
    }
control_results={}
def exact_accessions(row):
    accessions={row['saccver']}
    for identifier in row['sallseqid'].split(';'):
        accessions.update(token for token in identifier.split('|') if token)
    return accessions
for spec in validation_control_specs:
    control=spec['id']
    accession=spec.get('expected_accession')
    matches=[
        row for row in rows
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
    'command_completed_successfully':success, 'result_archive_valid':valid_json,
    'attempt_count':len(attempt_history), 'attempt_history':attempt_history,
    'termination_reason':(out/'TERMINATION_REASON.txt').read_text().strip(),
    'search_budget_seconds':int((out/'SEARCH_BUDGET_SECONDS.txt').read_text().strip()),
    'fatal_stderr_markers':fatal_markers, 'unexpected_result_query_ids':unexpected,
    'result_query_length_mismatches':bad_qlen, 'technical_complete':command_valid,
    'result_row_count':len(rows), 'per_query':per_query,
    'interpretation_boundary':'An empty hit table is evidence only when technical_complete is true; no-hit or divergence does not establish a new taxon.',
}
if json_error: status['json_validation_error']=json_error
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
python - "$OUT/SEARCH_STATUS.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1]))['technical_complete'] else 1)
PY
