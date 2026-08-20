#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash lantern_submission/run_frozen_v7.sh \
    --mapping /private/active_mapping.tsv \
    --manifest /private/input_manifest.tsv \
    --out /private/work \
    --threads 32 --memory-gb 256 \
    --expected-individuals 9 --expected-timepoints 4 [--dry-run]

The script never downloads or uploads restricted data. All paths must be on private
local/HPC storage. The frozen rule is not changed in response to any output.
EOF
}

mapping=""
manifest=""
out=""
threads=""
memory_gb=""
expected_individuals=""
expected_timepoints=""
dry_run=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mapping) mapping="$2"; shift 2 ;;
    --manifest) manifest="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --memory-gb) memory_gb="$2"; shift 2 ;;
    --expected-individuals) expected_individuals="$2"; shift 2 ;;
    --expected-timepoints) expected_timepoints="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in mapping manifest out threads memory_gb expected_individuals expected_timepoints; do
  if [[ -z "${!value}" ]]; then
    echo "Missing required option: $value" >&2
    usage >&2
    exit 2
  fi
done

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mapping=$(realpath "$mapping")
manifest=$(realpath "$manifest")
out=$(mkdir -p "$out" && realpath "$out")
rule="$repo_root/lantern_submission/configs/frozen_v7_rule.json"
freeze_root="$out/freeze"
plan_root="$out/plan"
submission_root="$out/submission"
mkdir -p "$freeze_root/mapping" "$freeze_root/inputs" "$plan_root" "$submission_root" "$out/provenance"

expected_samples=$(python - "$manifest" <<'PY'
import csv,sys
with open(sys.argv[1], newline='', encoding='utf-8-sig') as fh:
    rows=list(csv.DictReader(fh,delimiter='\t'))
print(','.join(row['sample_id'].strip() for row in rows if row.get('sample_id')))
PY
)

python "$repo_root/lantern_submission/scripts/load_explicit_mapping.py" \
  --mapping "$mapping" \
  --manifest "$manifest" \
  --expected-individuals "$expected_individuals" \
  --expected-timepoints "$expected_timepoints" \
  --expected-sample-ids "$expected_samples" \
  --require-consecutive-timepoints \
  --out "$freeze_root/mapping"

python "$repo_root/lantern_submission/scripts/validate_inputs.py" \
  --manifest "$manifest" \
  --mapping-freeze "$freeze_root/mapping/MAPPING_FREEZE.json" \
  --out "$freeze_root/inputs"

python "$repo_root/lantern_submission/scripts/build_execution_plan.py" \
  --mapping-freeze "$freeze_root/mapping/MAPPING_FREEZE.json" \
  --input-freeze "$freeze_root/inputs/INPUT_MANIFEST_FREEZE.json" \
  --rule "$rule" \
  --work-dir "$out/work" \
  --threads "$threads" \
  --memory-gb "$memory_gb" \
  --out "$plan_root"

python - "$repo_root" "$out" "$mapping" "$manifest" "$rule" <<'PY'
import hashlib,json,platform,sys
from pathlib import Path
repo,out,mapping,manifest,rule=map(Path,sys.argv[1:])
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()
obj={
  'status':'ACTIVE_PROTOCOL_FROZEN_BEFORE_EXECUTION',
  'repository_path':str(repo),
  'mapping_path':str(mapping),
  'mapping_sha256':sha(mapping),
  'input_manifest_path':str(manifest),
  'input_manifest_sha256':sha(manifest),
  'rule_path':str(rule),
  'rule_sha256':sha(rule),
  'platform':platform.platform(),
  'python':platform.python_version(),
  'post_result_tuning_allowed':False,
  'restricted_data_storage':'private local/HPC only',
}
(out/'ACTIVE_FREEZE_MANIFEST.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
PY

{
  echo "python=$(python --version 2>&1 || true)"
  echo "megahit=$(megahit --version 2>&1 | head -1 || true)"
  echo "spades=$(spades.py --version 2>&1 | head -1 || metaspades.py --version 2>&1 | head -1 || true)"
  echo "flye=$(flye --version 2>&1 | head -1 || true)"
  echo "minimap2=$(minimap2 --version 2>&1 | head -1 || true)"
  echo "samtools=$(samtools --version 2>&1 | head -1 || true)"
  echo "seqkit=$(seqkit version 2>&1 | head -1 || true)"
  echo "pigz=$(pigz --version 2>&1 | head -1 || true)"
} > "$out/provenance/SOFTWARE_VERSIONS.txt"

if [[ "$dry_run" == true ]]; then
  echo "DRY_RUN_COMPLETE"
  echo "Plan: $plan_root/RUN_PLAN.sh"
  exit 0
fi

(
  cd "$repo_root"
  bash "$plan_root/RUN_PLAN.sh"
)

printf 'individual_id\tfasta\n' > "$submission_root/INDIVIDUAL_FASTAS.tsv"
python - "$out" "$submission_root/INDIVIDUAL_FASTAS.tsv" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); dest=Path(sys.argv[2])
freeze=json.loads((root/'freeze/mapping/MAPPING_FREEZE.json').read_text())
with dest.open('a',encoding='utf-8') as f:
    for group in freeze['groups']:
        ind=str(group['individual_id'])
        fasta=root/f'work/individual_{ind}/submission/LANTERN_V7_INDIVIDUAL_{ind}.fasta'
        if not fasta.is_file() or fasta.stat().st_size==0:
            raise SystemExit(f'missing individual FASTA: {fasta}')
        f.write(f'{ind}\t{fasta}\n')
PY

python "$repo_root/lantern_submission/scripts/combine_submission_fastas.py" \
  --manifest "$submission_root/INDIVIDUAL_FASTAS.tsv" \
  --out-fasta "$submission_root/LANTERN_CAMI3_ASSEMBLY.fasta" \
  --out-json "$submission_root/COMBINATION_SUMMARY.json"

python "$repo_root/lantern_cami3/scripts/verify_cami_assembly.py" \
  "$submission_root/LANTERN_CAMI3_ASSEMBLY.fasta" \
  --out "$submission_root/VALIDATION.json"
sha256sum "$submission_root/LANTERN_CAMI3_ASSEMBLY.fasta" > "$submission_root/SHA256.txt"
find "$out" -type f ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > "$out/MANIFEST.sha256"
echo "FROZEN_RUN_COMPLETE_AWAITING_EXTERNAL_SUBMISSION_APPROVAL"
