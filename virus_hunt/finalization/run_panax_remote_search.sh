#!/usr/bin/env bash
# Current NCBI remote-database audit for immutable Panax A1/A2/B queries.
# A zero-hit result is accepted only when BLAST produced a valid archive,
# recovered every exact query, and reported nonzero statistics. Modes with
# embedded positive controls must also recover both controls.
set -Eeuo pipefail

MODE="${1:?search mode required}"
QUERY_ROOT="${2:?query directory required}"
OUT="${3:?output directory required}"
mkdir -p "$OUT"
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi

# Standard-task nr/nt coverage is split into explicit viral and nonviral
# Entrez partitions. The unfiltered remote service can emit zero-statistic
# archives that look successful but contain no usable search result.
case "$MODE" in
  protein_viral)
    program=blastp; query="$QUERY_ROOT/panax_candidates_plus_controls_orfs.faa"; database=nr
    extra=(-entrez_query 'txid10239[ORGN]' -seg yes -comp_based_stats 2)
    ;;
  protein_nonviral)
    program=blastp; query="$QUERY_ROOT/panax_three_partial_orfs.faa"; database=nr
    extra=(-entrez_query 'NOT txid10239[ORGN]' -seg yes -comp_based_stats 2)
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
    program=blastn; query="$QUERY_ROOT/panax_three_contigs.fna"; database=nt
    extra=(-task blastn -entrez_query 'NOT txid10239[ORGN]' -dust yes -soft_masking true)
    ;;
  nt_megablast)
    program=blastn; query="$QUERY_ROOT/panax_candidates_plus_controls_contigs.fna"; database=nt
    extra=(-task megablast -dust yes -soft_masking true)
    ;;
  nt_panax)
    program=blastn; query="$QUERY_ROOT/panax_three_contigs.fna"; database=nt
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

# Emit immutable expected query metadata before making a network request.
python - "$query" "$program" "$OUT/EXPECTED_QUERIES.json" <<'PY'
from pathlib import Path
import hashlib,json,re,sys
p=Path(sys.argv[1]); program=sys.argv[2]; out=Path(sys.argv[3])
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
required=(expected|controls) if p.name.startswith('panax_candidates_plus_controls_') else expected
if set(records)!=required:
    raise SystemExit(f'query set mismatch: observed={sorted(records)}, required={sorted(required)}')
allowed=r'[ACGTN]+' if program=='blastn' else r'[ABCDEFGHIKLMNPQRSTVWXYZ]+'
for name,seq in records.items():
    if not seq or not re.fullmatch(allowed,seq):
        raise SystemExit(f'empty or invalid query sequence: {name}')
payload={'query_file':str(p),'query_file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
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
  python - "$result_json" "$expected_json" "$mode" <<'PY'
from pathlib import Path
import json,math,sys
result_path=Path(sys.argv[1]); expected_path=Path(sys.argv[2]); mode=sys.argv[3]
payload=json.loads(result_path.read_text())
expected_payload=json.loads(expected_path.read_text())
expected_lengths={x['id']:int(x['length']) for x in expected_payload['queries']}
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
        hit_counts[qid]=len(hits)
if observed != expected_lengths:
    errors.append(f'query reports mismatch: observed={observed}, expected={expected_lengths}')
control_modes={'protein_viral','nt_viral','nt_megablast'}
if mode in control_modes:
    for control in ('PNX_Duplo_A_control','PNX_Duplo_B_control'):
        if hit_counts.get(control,0) < 1:
            errors.append(f'positive control has no hit: {control}')
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
for attempt in 1 2; do
  attempts="$attempt"
  attempt_stdout="$OUT/STDOUT.attempt${attempt}.txt"
  attempt_stderr="$OUT/STDERR.attempt${attempt}.txt"
  printf '[%s] attempt=%s\n' "$(date -u +%FT%TZ)" "$attempt" > "$attempt_stderr"
  : > "$attempt_stdout"
  rm -f "$OUT/RESULTS.asn" "$OUT/RESULTS.json" "$OUT/HITS.tsv"
  if timeout 105m "${cmd[@]}" >> "$attempt_stdout" 2>> "$attempt_stderr" && \
     [[ -s "$OUT/RESULTS.asn" ]] && \
     blast_formatter -archive "$OUT/RESULTS.asn" -outfmt 15 -out "$OUT/RESULTS.json" \
       >> "$attempt_stdout" 2>> "$attempt_stderr" && \
     blast_formatter -archive "$OUT/RESULTS.asn" \
       -outfmt '6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore qcovs staxids sscinames stitle sseq' \
       -out "$OUT/HITS.tsv" >> "$attempt_stdout" 2>> "$attempt_stderr" && \
     [[ -s "$OUT/RESULTS.json" ]] && \
     validate_remote_archive "$OUT/RESULTS.json" "$OUT/EXPECTED_QUERIES.json" "$MODE" \
       >> "$attempt_stdout" 2>> "$attempt_stderr"; then
    success=1
    success_attempt="$attempt"
    break
  fi
  (( attempt < 2 )) && sleep $((attempt * 30))
done
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
fields=['qseqid','saccver','pident','length','qlen','slen','qstart','qend','sstart','send',
        'evalue','bitscore','qcovs','staxids','sscinames','stitle','sseq']
rows=[]
for values in csv.reader((out/'HITS.tsv').open(errors='replace'),delimiter='\t'):
    if not values: continue
    if len(values) != len(fields):
        raise SystemExit(f'malformed BLAST result row with {len(values)} fields')
    rows.append(dict(zip(fields,values)))
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
status={
    'generated_utc':datetime.now(timezone.utc).isoformat(), 'mode':mode, 'database':database,
    'query_file':str(query), 'query_sha256':hashlib.sha256(query.read_bytes()).hexdigest(),
    'query_count':len(expected_lengths), 'query_ids':list(expected_lengths),
    'expected_query_lengths':expected_lengths, 'json_query_lengths':json_queries,
    'command_completed_successfully':success, 'result_archive_valid':valid_json,
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
