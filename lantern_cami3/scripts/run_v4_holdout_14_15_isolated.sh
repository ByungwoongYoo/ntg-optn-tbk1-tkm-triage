#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-v4_holdout_14_15_isolated}"
HOLDOUT_CONFIG="${HOLDOUT_CONFIG:-lantern_cami3/config/lantern_v4_holdout_14_15_20260820.json}"
READ_ARTIFACT_ROOT="${READ_ARTIFACT_ROOT:-inputs/read_artifacts}"
BASELINE_ARTIFACT_ROOT="${BASELINE_ARTIFACT_ROOT:-inputs/baseline_artifacts}"
THREADS="${THREADS:-4}"

mkdir -p "$ROOT"/{inputs,reads,candidates,evidence,configs,variants,evaluation,pseudo_novel,abundance,decision,submission,validation,provenance,logs,failure_audit,truth_download,truth}

# Extract the exact frozen v4 configuration selected on samples 0/1.
python - "$HOLDOUT_CONFIG" "$ROOT/configs/v4_balanced.json" <<'PY'
import json,sys
from pathlib import Path
src=Path(sys.argv[1]); dst=Path(sys.argv[2])
x=json.loads(src.read_text())
assert x['holdout_pair']==[14,15], x
assert x['official_pairing_correction']['correct_pair']==[14,15], x
assert x['backbone']=='metaspades_hybrid_pair', x
assert x['frozen_config']['version']=='LANTERN-v4-balanced-20260820', x
assert x['frozen_config']['truth_blind'] is True, x
assert x['method_origin'].startswith('The v4_balanced configuration was selected on Toy samples 0/1'), x
assert x['read_input']=={'short_fraction':0.08,'long_fraction':0.04,'seed':20260819}, x
assert x['strict_success_gates']['genome_fraction_gain_minimum_percentage_points']==0.5, x
assert x['strict_success_gates']['mean_genome_recovery_gain_minimum_percentage_points']==1.0, x
assert x['strict_success_gates']['low_abundance_gain_minimum_percentage_points']==2.0, x
assert x['strict_success_gates']['longitudinal_ablation_drop_minimum_percentage_points']==0.5, x
dst.write_text(json.dumps(x['frozen_config'],indent=2,sort_keys=True)+'\n')
PY

# Isolate immutable reads created before any Toy 14/15 gold access.
for s in 14 15; do
  r1=$(find "$READ_ARTIFACT_ROOT" -type f -name "S${s}_R1.fastq.gz" -print -quit)
  test -s "$r1"
  d=$(dirname "$r1")
  for kind in R1 R2 long; do
    src="$d/S${s}_${kind}.fastq.gz"
    test -s "$src"
    ln -s "$PWD/$src" "$ROOT/reads/S${s}_${kind}.fastq.gz"
    gzip -t "$ROOT/reads/S${s}_${kind}.fastq.gz"
  done
  boundary=$(find "$d" -type f -path '*/provenance/FROZEN_INPUTS.txt' -print -quit || true)
  if test -n "$boundary"; then
    grep -q '^gold_accessed=false$' "$boundary"
    cp "$boundary" "$ROOT/provenance/READ_S${s}_FROZEN_INPUTS.txt"
  fi
done

# Isolate all strong baseline assemblies generated truth-blind in the upstream run.
methods=(megahit_short_pair metaspades_short_pair metaspades_hybrid_pair flye_long_pair flye_long_s14 flye_long_s15)
printf 'method\tavailability\tsource\n' > "$ROOT/provenance/BASELINE_AVAILABILITY.tsv"
for method in "${methods[@]}"; do
  f=$(find "$BASELINE_ARTIFACT_ROOT" -type f -path "*${method}*/final_assembly.fasta" -print -quit || true)
  if test -n "$f" && test -s "$f"; then
    cp "$f" "$ROOT/inputs/${method}.fasta"
    printf '%s\tavailable\t%s\n' "$method" "$f" >> "$ROOT/provenance/BASELINE_AVAILABILITY.tsv"
    status=$(find "$(dirname "$f")" -type f -name METHOD_STATUS.json -print -quit || true)
    if test -n "$status"; then
      grep -q '"truth_accessed"[[:space:]]*:[[:space:]]*false' "$status"
      cp "$status" "$ROOT/provenance/${method}_METHOD_STATUS.json"
    fi
  else
    printf '%s\tmissing\t\n' "$method" >> "$ROOT/provenance/BASELINE_AVAILABILITY.tsv"
  fi
done
for required in megahit_short_pair metaspades_short_pair metaspades_hybrid_pair flye_long_pair; do
  test -s "$ROOT/inputs/${required}.fasta"
done

cat > "$ROOT/inputs/source_manifest.tsv" <<'EOF'
source_id	assembler	scope	mode
megahit_short_pair	megahit	longitudinal	short
metaspades_short_pair	metaspades	longitudinal	short
flye_long_pair	flye	longitudinal	long
flye_long_s14	flye	timepoint_support	long
flye_long_s15	flye	timepoint_support	long
EOF

git rev-parse HEAD > "$ROOT/provenance/GIT_HEAD.txt"
uname -a > "$ROOT/provenance/UNAME.txt"
python --version > "$ROOT/provenance/PYTHON_VERSION.txt"
minimap2 --version > "$ROOT/provenance/MINIMAP2_VERSION.txt"
samtools --version > "$ROOT/provenance/SAMTOOLS_VERSION.txt"
mash --version > "$ROOT/provenance/MASH_VERSION.txt"
seqkit version > "$ROOT/provenance/SEQKIT_VERSION.txt" 2>&1 || true
micromamba list > "$ROOT/provenance/MICROMAMBA_LIST.txt" 2>&1 || true
cp "$HOLDOUT_CONFIG" "$ROOT/provenance/"
sha256sum "$ROOT"/inputs/*.fasta "$ROOT"/reads/*.fastq.gz "$ROOT/configs/v4_balanced.json" > "$ROOT/provenance/IMMUTABLE_INPUT_SHA256SUMS.txt"
seqkit stats -T "$ROOT"/inputs/*.fasta "$ROOT"/reads/*.fastq.gz > "$ROOT/provenance/INPUT_STATS.tsv"
printf 'holdout_pair=14,15\ntruth_accessed=false\nsource_read_run=32332720846\nsource_baseline_run=32332720846\n' > "$ROOT/provenance/ISOLATED_EXECUTION_BOUNDARY.txt"

# Build non-backbone candidate pool.
assembly_args=(
  --assembly "megahit_short_pair=$ROOT/inputs/megahit_short_pair.fasta"
  --assembly "metaspades_short_pair=$ROOT/inputs/metaspades_short_pair.fasta"
  --assembly "flye_long_pair=$ROOT/inputs/flye_long_pair.fasta"
)
test -s "$ROOT/inputs/flye_long_s14.fasta" && assembly_args+=(--assembly "flye_long_s14=$ROOT/inputs/flye_long_s14.fasta")
test -s "$ROOT/inputs/flye_long_s15.fasta" && assembly_args+=(--assembly "flye_long_s15=$ROOT/inputs/flye_long_s15.fasta")
python lantern_cami3/scripts/prepare_candidates.py "${assembly_args[@]}" \
  --source-manifest "$ROOT/inputs/source_manifest.tsv" --min-length 1000 --out "$ROOT/candidates/prepared" \
  > "$ROOT/logs/PREPARE_CANDIDATES.log" 2>&1
minimap2 -x asm5 -c -N 20 -t "$THREADS" "$ROOT/candidates/prepared/combined_candidates.fasta" \
  "$ROOT/candidates/prepared/combined_candidates.fasta" > "$ROOT/candidates/all_vs_all.paf" \
  2> "$ROOT/logs/CANDIDATE_ALL_VS_ALL.log"
python lantern_cami3/scripts/cluster_candidates.py \
  --fasta "$ROOT/candidates/prepared/combined_candidates.fasta" \
  --metadata "$ROOT/candidates/prepared/candidates.tsv" --paf "$ROOT/candidates/all_vs_all.paf" \
  --min-identity 0.97 --min-shorter-coverage 0.85 --out "$ROOT/candidates/clusters" \
  > "$ROOT/logs/CLUSTER_CANDIDATES.log" 2>&1
minimap2 -x asm5 -c -N 100 -t "$THREADS" "$ROOT/inputs/metaspades_hybrid_pair.fasta" \
  "$ROOT/candidates/clusters/cluster_representatives.fasta" > "$ROOT/candidates/candidate_to_backbone.paf" \
  2> "$ROOT/logs/CANDIDATE_TO_BACKBONE.log"
test -s "$ROOT/candidates/candidate_to_backbone.paf"
rm -f "$ROOT/candidates/all_vs_all.paf" "$ROOT/candidates/prepared/combined_candidates.fasta"

# Build per-timepoint short- and long-read evidence.
reps="$ROOT/candidates/clusters/cluster_representatives.fasta"
for s in 14 15; do
  minimap2 -t "$THREADS" -ax sr "$reps" "$ROOT/reads/S${s}_R1.fastq.gz" "$ROOT/reads/S${s}_R2.fastq.gz" \
    2> "$ROOT/logs/MAP_S${s}_SHORT.log" | samtools sort -@ 2 -o "$ROOT/evidence/S${s}_short.bam" -
  samtools index "$ROOT/evidence/S${s}_short.bam"
  samtools coverage "$ROOT/evidence/S${s}_short.bam" > "$ROOT/evidence/S${s}_short_coverage.tsv"
  minimap2 -t "$THREADS" -ax map-ont "$reps" "$ROOT/reads/S${s}_long.fastq.gz" \
    2> "$ROOT/logs/MAP_S${s}_LONG.log" | samtools sort -@ 2 -o "$ROOT/evidence/S${s}_long.bam" -
  samtools index "$ROOT/evidence/S${s}_long.bam"
  samtools coverage "$ROOT/evidence/S${s}_long.bam" > "$ROOT/evidence/S${s}_long_coverage.tsv"
  minimap2 -t "$THREADS" -x map-ont -c -N 10 "$reps" "$ROOT/reads/S${s}_long.fastq.gz" \
    > "$ROOT/evidence/S${s}_long.paf" 2> "$ROOT/logs/MAP_S${s}_LONG_PAF.log"
done
python lantern_cami3/scripts/merge_mapping_evidence.py \
  --short-coverage "S14=$ROOT/evidence/S14_short_coverage.tsv" \
  --short-coverage "S15=$ROOT/evidence/S15_short_coverage.tsv" \
  --long-coverage "S14=$ROOT/evidence/S14_long_coverage.tsv" \
  --long-coverage "S15=$ROOT/evidence/S15_long_coverage.tsv" \
  --long-paf "S14=$ROOT/evidence/S14_long.paf" \
  --long-paf "S15=$ROOT/evidence/S15_long.paf" \
  --end-margin 250 --min-mapq 20 --min-align 500 --out "$ROOT/evidence/merged" \
  > "$ROOT/logs/MERGE_EVIDENCE.log" 2>&1
rm -f "$ROOT/evidence"/*.bam "$ROOT/evidence"/*.bai "$ROOT/evidence"/*.paf
rm -rf "$ROOT/reads"

# Construct full method and all predeclared ablations before gold access.
for ablation in full no_longitudinal no_long no_consensus; do
  out="$ROOT/variants/v4_balanced_${ablation}"
  python lantern_cami3/scripts/baseline_preserving_augment.py \
    --backbone "$ROOT/inputs/metaspades_hybrid_pair.fasta" \
    --candidates "$ROOT/candidates/clusters/cluster_representatives.fasta" \
    --metadata "$ROOT/candidates/clusters/representative_metadata.tsv" \
    --evidence "$ROOT/evidence/merged/mapping_evidence.tsv" \
    --candidate-to-backbone-paf "$ROOT/candidates/candidate_to_backbone.paf" \
    --config "$ROOT/configs/v4_balanced.json" --ablation "$ablation" --out "$out" \
    > "$ROOT/logs/BUILD_${ablation}.log" 2>&1
  python lantern_cami3/scripts/verify_cami_assembly.py "$out/LANTERN_BACKBONE_AUGMENTED.fasta" \
    --out "$out/VALIDATION.json"
done

# Immutable truth-blind freeze gate.
python - "$ROOT" "$HOLDOUT_CONFIG" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
root=Path(sys.argv[1]); cfg=Path(sys.argv[2]); files=[]
include=[
 cfg, root/'configs/v4_balanced.json', root/'inputs/metaspades_hybrid_pair.fasta',
 root/'candidates/clusters/representative_metadata.tsv', root/'candidates/clusters/cluster_members.tsv',
 root/'evidence/merged/mapping_evidence.tsv', Path('lantern_cami3/scripts/baseline_preserving_augment.py'),
 Path('lantern_cami3/scripts/evaluate_frozen_holdout.py'), Path('lantern_cami3/scripts/run_v4_holdout_14_15_isolated.sh')
]
include += sorted(root.glob('variants/*/LANTERN_BACKBONE_AUGMENTED.fasta'))
for p in include:
    if not p.is_file(): raise SystemExit(f'missing freeze input: {p}')
    files.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
obj={
 'freeze_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()), 'truth_accessed':False,
 'holdout_pair':[14,15], 'n_files':len(files), 'files':files,
 'boundary':'All corrected Toy pair 14/15 inputs, evidence, v4 output assemblies, ablations, code and thresholds were frozen before either GSA was downloaded.'
}
(root/'provenance/V4_HOLDOUT_TRUTH_BLIND_FREEZE.json').write_text(json.dumps(obj,indent=2)+'\n')
PY
sha256sum "$ROOT/provenance/V4_HOLDOUT_TRUTH_BLIND_FREEZE.json" > "$ROOT/provenance/V4_HOLDOUT_TRUTH_BLIND_FREEZE_SHA256.txt"
rm -f "$ROOT/candidates/candidate_to_backbone.paf" "$ROOT/candidates/clusters/cluster_representatives.fasta"

# Only now open public Toy gold.
for s in 14 15; do
  curl -fsSL --retry 8 --retry-all-errors \
    "https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_${s}_gsa.tar.gz" \
    -o "$ROOT/truth_download/sample_${s}_gsa.tar.gz"
done
sha256sum "$ROOT/truth_download"/*.tar.gz > "$ROOT/provenance/GSA_SHA256SUMS.txt"
python lantern_cami3/scripts/extract_gsa_truth.py \
  --gsa-tar "S14=$ROOT/truth_download/sample_14_gsa.tar.gz" \
  --gsa-tar "S15=$ROOT/truth_download/sample_15_gsa.tar.gz" --out "$ROOT/truth" \
  > "$ROOT/logs/EXTRACT_TRUTH.log" 2>&1
minimap2 -I 4G -d "$ROOT/truth/toy_truth.mmi" "$ROOT/truth/combined_truth.fasta" \
  > "$ROOT/logs/TRUTH_INDEX.stdout" 2> "$ROOT/logs/TRUTH_INDEX.log"

cat > "$ROOT/provenance/METHOD_MANIFEST.tsv" <<'EOF'
method	role	mode	scope	ablation	eligible_primary_baseline
megahit_short_pair	baseline	short	longitudinal	none	true
metaspades_short_pair	baseline	short	longitudinal	none	true
metaspades_hybrid_pair	baseline	hybrid	longitudinal	none	true
flye_long_pair	baseline	long	longitudinal	none	true
flye_long_s14	baseline	long	single	none	true
flye_long_s15	baseline	long	single	none	true
v4_balanced	lantern	hybrid	longitudinal	full	false
v4_balanced_no_longitudinal	ablation	hybrid	longitudinal	no_longitudinal	false
v4_balanced_no_long	ablation	hybrid	longitudinal	no_long	false
v4_balanced_no_consensus	ablation	hybrid	longitudinal	no_consensus	false
EOF
printf 'method\tfasta\n' > "$ROOT/provenance/EVALUATED_METHODS.tsv"
for method in "${methods[@]}"; do
  test -s "$ROOT/inputs/${method}.fasta" || continue
  printf '%s\t%s\n' "$method" "$ROOT/inputs/${method}.fasta" >> "$ROOT/provenance/EVALUATED_METHODS.tsv"
done
for ablation in full no_longitudinal no_long no_consensus; do
  method=v4_balanced
  test "$ablation" = full || method="v4_balanced_${ablation}"
  printf '%s\t%s\n' "$method" "$ROOT/variants/v4_balanced_${ablation}/LANTERN_BACKBONE_AUGMENTED.fasta" >> "$ROOT/provenance/EVALUATED_METHODS.tsv"
done

tail -n +2 "$ROOT/provenance/EVALUATED_METHODS.tsv" | while IFS=$'\t' read -r method fasta; do
  dir="$ROOT/evaluation/$method"; mkdir -p "$dir"
  minimap2 -x asm5 -c -N 100 -t "$THREADS" "$ROOT/truth/toy_truth.mmi" "$fasta" \
    > "$dir/assembly_to_truth.paf" 2> "$dir/minimap2.log"
  python lantern_cami3/scripts/evaluate_gold_coverage.py \
    --paf "$dir/assembly_to_truth.paf" --truth-mapping "$ROOT/truth/truth_mapping.tsv" \
    --assembly-fasta "$fasta" --min-identity 0.90 --min-alignment 500 --out "$dir" \
    > "$dir/evaluate_console.log" 2>&1
  rm -f "$dir/assembly_to_truth.paf"
done

# Pseudo-novel evaluation-only tiers.
python - "$ROOT" <<'PY'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1]); rows=list(csv.DictReader(open(root/'truth/genome_reference_manifest.tsv'),delimiter='\t')); paths=[]
for r in rows:
    if int(r['truth_bp'])<50000: continue
    p=root/'truth/references'/Path(r['reference_path']).name
    if p.is_file(): paths.append(str(p))
(root/'pseudo_novel/mash_reference_files.txt').write_text('\n'.join(paths)+'\n')
PY
mash sketch -k 21 -s 1000 -l "$ROOT/pseudo_novel/mash_reference_files.txt" \
  -o "$ROOT/pseudo_novel/truth_genomes" > "$ROOT/logs/MASH_SKETCH.log" 2>&1
mash dist "$ROOT/pseudo_novel/truth_genomes.msh" "$ROOT/pseudo_novel/truth_genomes.msh" \
  > "$ROOT/pseudo_novel/truth_genome_distances.tsv"
recovery_args=()
while IFS=$'\t' read -r method fasta; do
  test "$method" = method && continue
  recovery_args+=(--recovery "$method=$ROOT/evaluation/$method/per_genome_recovery.tsv")
done < "$ROOT/provenance/EVALUATED_METHODS.tsv"
python lantern_cami3/scripts/build_pseudo_novel_audit.py \
  --manifest "$ROOT/truth/genome_reference_manifest.tsv" \
  --mash "$ROOT/pseudo_novel/truth_genome_distances.tsv" "${recovery_args[@]}" \
  --min-truth-bp 50000 --out "$ROOT/pseudo_novel" > "$ROOT/logs/PSEUDO_NOVEL.log" 2>&1

# Read-abundance stratification for the predefined low-abundance gate.
for s in 14 15; do
  python lantern_cami3/scripts/extract_read_abundance.py \
    --url "https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_${s}_reads.tar.gz" \
    --sample-id "sample${s}" --paired --out "$ROOT/abundance/sample${s}_abundance.tsv" \
    > "$ROOT/logs/ABUNDANCE_S${s}.log" 2>&1
done
python lantern_cami3/scripts/summarize_abundance_recovery.py \
  --abundance "$ROOT/abundance/sample14_abundance.tsv" \
  --abundance "$ROOT/abundance/sample15_abundance.tsv" "${recovery_args[@]}" \
  --low-percent 0.01 --high-percent 0.1 --out "$ROOT/abundance/results" \
  > "$ROOT/logs/ABUNDANCE_SUMMARY.log" 2>&1

# Apply the frozen strict decision gate.
method_args=()
while IFS=$'\t' read -r method fasta; do
  test "$method" = method && continue
  method_args+=(--method "$method=$ROOT/evaluation/$method")
done < "$ROOT/provenance/EVALUATED_METHODS.tsv"
python lantern_cami3/scripts/evaluate_frozen_holdout.py \
  "${method_args[@]}" --manifest "$ROOT/provenance/METHOD_MANIFEST.tsv" \
  --holdout-config "$HOLDOUT_CONFIG" --full-method v4_balanced \
  --no-longitudinal-method v4_balanced_no_longitudinal --no-long-method v4_balanced_no_long \
  --no-consensus-method v4_balanced_no_consensus \
  --abundance-summary "$ROOT/abundance/results/ABUNDANCE_STRATIFIED_SUMMARY.tsv" \
  --pseudo-summary "$ROOT/pseudo_novel/PSEUDO_NOVEL_RECOVERY_SUMMARY.tsv" --out "$ROOT/decision" \
  > "$ROOT/logs/HOLDOUT_DECISION.log" 2>&1

# CAMI-compatible candidate output for audit; this is not submitted automatically.
cp "$ROOT/variants/v4_balanced_full/LANTERN_BACKBONE_AUGMENTED.fasta" "$ROOT/submission/LANTERN_V4_HOLDOUT_14_15.fasta"
python lantern_cami3/scripts/verify_cami_assembly.py "$ROOT/submission/LANTERN_V4_HOLDOUT_14_15.fasta" \
  --out "$ROOT/validation/SUBMISSION_FASTA_VALIDATION.json"
sha256sum "$ROOT/submission/LANTERN_V4_HOLDOUT_14_15.fasta" > "$ROOT/validation/SUBMISSION_FASTA_SHA256.txt"

# Compact but retain all decision-critical evidence and submission FASTA.
rm -rf "$ROOT/truth_download" "$ROOT/truth/references"
rm -f "$ROOT/truth/combined_truth.fasta" "$ROOT/truth/toy_truth.mmi"
rm -f "$ROOT/pseudo_novel/truth_genomes.msh" "$ROOT/pseudo_novel/truth_genome_distances.tsv" "$ROOT/pseudo_novel/mash_reference_files.txt"
rm -f "$ROOT/inputs"/*.fasta
rm -f "$ROOT/candidates/prepared/candidates.tsv" "$ROOT/evidence"/*_coverage.tsv
for d in "$ROOT"/variants/v4_balanced_no_*; do rm -f "$d/LANTERN_BACKBONE_AUGMENTED.fasta"; done
find "$ROOT/evaluation" -type f -name per_truth_contig.tsv -delete 2>/dev/null || true
find "$ROOT" -type f ! -name SHA256SUMS_ALL.txt -print0 | sort -z | xargs -0 sha256sum > "$ROOT/SHA256SUMS_ALL.txt"
find "$ROOT" -type f -printf '%P\t%s\n' | sort > "$ROOT/CONTENTS.tsv"
du -ah "$ROOT" | sort -h > "$ROOT/FILE_SIZES.txt"

cat "$ROOT/decision/HOLDOUT_DECISION.md"
