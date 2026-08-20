#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-lantern_short_benchmark_v2}"
FRACTION="${LANTERN_SHORT_FRACTION:-0.03}"
SEED="${LANTERN_SEED:-20260819}"
THREADS="${LANTERN_THREADS:-4}"
MEM_GB="${LANTERN_MEMORY_GB:-12}"
SCRIPTS="$PWD/lantern_cami3/scripts"
CONFIG="$PWD/lantern_cami3/config/lantern_frozen_v1.json"
DECISION="$PWD/lantern_cami3/config/decision_gate_v1.json"
export PYTHONPATH="$SCRIPTS"
mkdir -p "$ROOT"/{reads,assemblies,logs,provenance,truth_blind,evaluation,pseudo_novel,decision,failure_audit,results}

SHORT0='https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_0_reads.tar.gz'
SHORT1='https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_1_reads.tar.gz'
GSA0='https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_0_gsa.tar.gz'
GSA1='https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_1_gsa.tar.gz'

record_versions() {
  uname -a > "$ROOT/provenance/UNAME.txt"
  python --version > "$ROOT/provenance/PYTHON_VERSION.txt" 2>&1
  megahit --version > "$ROOT/provenance/MEGAHIT_VERSION.txt" 2>&1 || true
  metaspades.py --version > "$ROOT/provenance/METASPADES_VERSION.txt" 2>&1 || true
  minimap2 --version > "$ROOT/provenance/MINIMAP2_VERSION.txt" 2>&1 || true
  samtools --version > "$ROOT/provenance/SAMTOOLS_VERSION.txt" 2>&1 || true
  mash --version > "$ROOT/provenance/MASH_VERSION.txt" 2>&1 || true
  micromamba list > "$ROOT/provenance/MICROMAMBA_LIST.txt" 2>&1 || true
  git rev-parse HEAD > "$ROOT/provenance/GIT_HEAD.txt"
  cp "$CONFIG" "$DECISION" "$ROOT/provenance/"
  printf 'short_fraction=%s\nseed=%s\nthreads=%s\nmemory_gb=%s\n' "$FRACTION" "$SEED" "$THREADS" "$MEM_GB" > "$ROOT/provenance/RUN_PARAMETERS.txt"
  printf '%s\n%s\n%s\n%s\n' "$SHORT0" "$SHORT1" "$GSA0" "$GSA1" > "$ROOT/provenance/SOURCE_URLS.txt"
  df -h > "$ROOT/provenance/DISK_BEFORE.txt"
}
record_versions

python "$SCRIPTS/stream_downsample_cami_tar.py" --url "$SHORT0" --sample-id sample0 --mode short --fraction "$FRACTION" --seed "$SEED" --out-dir "$ROOT/reads" --min-kept 10000 2>&1 | tee "$ROOT/logs/DOWNSAMPLE_SHORT0.log"
python "$SCRIPTS/stream_downsample_cami_tar.py" --url "$SHORT1" --sample-id sample1 --mode short --fraction "$FRACTION" --seed "$SEED" --out-dir "$ROOT/reads" --min-kept 10000 2>&1 | tee "$ROOT/logs/DOWNSAMPLE_SHORT1.log"
sha256sum "$ROOT"/reads/*.fastq.gz "$ROOT"/reads/*_downsample.json > "$ROOT/provenance/DOWNSAMPLED_READ_SHA256.txt"

run_assembler() {
  local name="$1"; shift
  local log="$ROOT/logs/${name}.log" status="$ROOT/failure_audit/${name}.status"
  set +e
  /usr/bin/time -v timeout --signal=TERM --kill-after=120s 150m "$@" >"$log" 2>&1
  local rc=$?
  set -e
  printf 'name=%s\nexit_code=%s\ncommand=' "$name" "$rc" > "$status"
  printf '%q ' "$@" >> "$status"; printf '\n' >> "$status"
  return 0
}

run_assembler megahit_single0 megahit -1 "$ROOT/reads/sample0_R1.fastq.gz" -2 "$ROOT/reads/sample0_R2.fastq.gz" -o "$ROOT/assemblies/megahit_single0" --min-contig-len 1000 -t "$THREADS" --memory 0.8
run_assembler megahit_single1 megahit -1 "$ROOT/reads/sample1_R1.fastq.gz" -2 "$ROOT/reads/sample1_R2.fastq.gz" -o "$ROOT/assemblies/megahit_single1" --min-contig-len 1000 -t "$THREADS" --memory 0.8
run_assembler megahit_longitudinal megahit -1 "$ROOT/reads/sample0_R1.fastq.gz,$ROOT/reads/sample1_R1.fastq.gz" -2 "$ROOT/reads/sample0_R2.fastq.gz,$ROOT/reads/sample1_R2.fastq.gz" -o "$ROOT/assemblies/megahit_longitudinal" --min-contig-len 1000 -t "$THREADS" --memory 0.8

# metaSPAdes accepts one paired-end library in meta mode. Concatenated gzip members
# preserve all reads while presenting the two longitudinal time points as one library.
cat "$ROOT/reads/sample0_R1.fastq.gz" "$ROOT/reads/sample1_R1.fastq.gz" > "$ROOT/reads/longitudinal_R1.fastq.gz"
cat "$ROOT/reads/sample0_R2.fastq.gz" "$ROOT/reads/sample1_R2.fastq.gz" > "$ROOT/reads/longitudinal_R2.fastq.gz"
gzip -t "$ROOT/reads/longitudinal_R1.fastq.gz"
gzip -t "$ROOT/reads/longitudinal_R2.fastq.gz"
run_assembler metaspades_longitudinal metaspades.py --only-assembler -t "$THREADS" -m "$MEM_GB" -1 "$ROOT/reads/longitudinal_R1.fastq.gz" -2 "$ROOT/reads/longitudinal_R2.fastq.gz" -o "$ROOT/assemblies/metaspades_longitudinal"

cat > "$ROOT/assemblies/source_manifest.tsv" <<'EOF_MANIFEST'
source_id	assembler	scope	mode
EOF_MANIFEST
ASSEMBLY_ARGS=()
n_assemblies=0
add_assembly() {
  local id="$1" assembler="$2" scope="$3" mode="$4" path="$5"
  if [[ -s "$path" ]]; then
    printf '%s\t%s\t%s\t%s\n' "$id" "$assembler" "$scope" "$mode" >> "$ROOT/assemblies/source_manifest.tsv"
    ASSEMBLY_ARGS+=(--assembly "$id=$path")
    sha256sum "$path" >> "$ROOT/provenance/BASELINE_FASTA_SHA256.txt"
    n_assemblies=$((n_assemblies+1))
  fi
  return 0
}
add_assembly megahit_s0 megahit single short "$ROOT/assemblies/megahit_single0/final.contigs.fa"
add_assembly megahit_s1 megahit single short "$ROOT/assemblies/megahit_single1/final.contigs.fa"
add_assembly megahit_pair megahit longitudinal short "$ROOT/assemblies/megahit_longitudinal/final.contigs.fa"
add_assembly metaspades_pair metaspades longitudinal short "$ROOT/assemblies/metaspades_longitudinal/contigs.fasta"
[[ "$n_assemblies" -ge 3 ]] || { echo 'Too few successful baseline assemblies' >&2; exit 2; }

mkdir -p "$ROOT/truth_blind/prepared" "$ROOT/truth_blind/clusters" "$ROOT/truth_blind/evidence"
python "$SCRIPTS/prepare_candidates.py" "${ASSEMBLY_ARGS[@]}" --source-manifest "$ROOT/assemblies/source_manifest.tsv" --min-length 1000 --out "$ROOT/truth_blind/prepared" 2>&1 | tee "$ROOT/logs/PREPARE_CANDIDATES.log"
minimap2 -x asm5 -c -N 100 -t "$THREADS" "$ROOT/truth_blind/prepared/combined_candidates.fasta" "$ROOT/truth_blind/prepared/combined_candidates.fasta" > "$ROOT/truth_blind/all_vs_all.paf" 2> "$ROOT/logs/ALL_VS_ALL_MINIMAP2.log"
python "$SCRIPTS/cluster_candidates.py" --fasta "$ROOT/truth_blind/prepared/combined_candidates.fasta" --metadata "$ROOT/truth_blind/prepared/candidates.tsv" --paf "$ROOT/truth_blind/all_vs_all.paf" --min-identity 0.97 --min-shorter-coverage 0.85 --out "$ROOT/truth_blind/clusters" 2>&1 | tee "$ROOT/logs/CLUSTER.log"

for sample in sample0 sample1; do
  minimap2 -t "$THREADS" -ax sr "$ROOT/truth_blind/clusters/cluster_representatives.fasta" "$ROOT/reads/${sample}_R1.fastq.gz" "$ROOT/reads/${sample}_R2.fastq.gz" 2> "$ROOT/logs/MAP_${sample}_SHORT.log" | samtools sort -@ 2 -o "$ROOT/truth_blind/evidence/${sample}_short.bam" -
  samtools index "$ROOT/truth_blind/evidence/${sample}_short.bam"
  samtools coverage "$ROOT/truth_blind/evidence/${sample}_short.bam" > "$ROOT/truth_blind/evidence/${sample}_short_coverage.tsv"
done
python "$SCRIPTS/merge_mapping_evidence.py" --short-coverage "sample0=$ROOT/truth_blind/evidence/sample0_short_coverage.tsv" --short-coverage "sample1=$ROOT/truth_blind/evidence/sample1_short_coverage.tsv" --out "$ROOT/truth_blind/evidence/merged" 2>&1 | tee "$ROOT/logs/MERGE_EVIDENCE.log"

for ablation in full no_longitudinal no_consensus; do
  python "$SCRIPTS/lantern_select.py" --fasta "$ROOT/truth_blind/clusters/cluster_representatives.fasta" --metadata "$ROOT/truth_blind/clusters/representative_metadata.tsv" --evidence "$ROOT/truth_blind/evidence/merged/mapping_evidence.tsv" --config "$CONFIG" --ablation "$ablation" --out "$ROOT/truth_blind/LANTERN_${ablation}" 2>&1 | tee "$ROOT/logs/LANTERN_${ablation}.log"
  python "$SCRIPTS/verify_cami_assembly.py" "$ROOT/truth_blind/LANTERN_${ablation}/LANTERN_ASSEMBLY.fasta" --out "$ROOT/truth_blind/LANTERN_${ablation}/VALIDATION.json"
done

python - "$ROOT" "$CONFIG" "$DECISION" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
root=Path(sys.argv[1]); files=[]
for base in [root/'assemblies',root/'truth_blind',Path(sys.argv[2]),Path(sys.argv[3])]:
    paths=[base] if base.is_file() else sorted(base.rglob('*'))
    for p in paths:
        if p.is_file() and not p.name.endswith('.bam') and not p.name.endswith('.bai'):
            files.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
out={'freeze_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'git_head':(root/'provenance/GIT_HEAD.txt').read_text().strip(),'truth_accessed':False,'n_files':len(files),'files':files,'boundary':'All assembly, clustering, mapping-evidence selection, thresholds, ablations and hashes were frozen before public Toy gold-standard files were downloaded.'}
(root/'truth_blind/TRUTH_BLIND_FREEZE.json').write_text(json.dumps(out,indent=2)+'\n')
PY
sha256sum "$ROOT/truth_blind/TRUTH_BLIND_FREEZE.json" > "$ROOT/provenance/TRUTH_BLIND_FREEZE_SHA256.txt"

mkdir -p "$ROOT/truth_download" "$ROOT/truth"
curl -fsSL --retry 8 --retry-all-errors "$GSA0" -o "$ROOT/truth_download/sample0_gsa.tar.gz"
curl -fsSL --retry 8 --retry-all-errors "$GSA1" -o "$ROOT/truth_download/sample1_gsa.tar.gz"
sha256sum "$ROOT/truth_download"/*.tar.gz > "$ROOT/provenance/GSA_INPUT_SHA256.txt"
python "$SCRIPTS/extract_gsa_truth.py" --gsa-tar "sample0=$ROOT/truth_download/sample0_gsa.tar.gz" --gsa-tar "sample1=$ROOT/truth_download/sample1_gsa.tar.gz" --out "$ROOT/truth" 2>&1 | tee "$ROOT/logs/EXTRACT_TRUTH.log"

cat > "$ROOT/METHOD_MANIFEST.tsv" <<'EOF_METHODS'
method	role	mode	scope	ablation	eligible_primary_baseline
megahit_s0	baseline	short	single	none	true
megahit_s1	baseline	short	single	none	true
megahit_pair	baseline	short	longitudinal	none	true
metaspades_pair	baseline	short	longitudinal	none	true
lantern_full	lantern	short	longitudinal	full	false
lantern_no_longitudinal	ablation	short	longitudinal	no_longitudinal	false
lantern_no_consensus	ablation	short	longitudinal	no_consensus	false
EOF_METHODS

METHODS=()
add_method() {
  local name="$1" path="$2"
  if [[ -s "$path" ]]; then METHODS+=("$name=$path"); fi
  return 0
}
add_method megahit_s0 "$ROOT/assemblies/megahit_single0/final.contigs.fa"
add_method megahit_s1 "$ROOT/assemblies/megahit_single1/final.contigs.fa"
add_method megahit_pair "$ROOT/assemblies/megahit_longitudinal/final.contigs.fa"
add_method metaspades_pair "$ROOT/assemblies/metaspades_longitudinal/contigs.fasta"
add_method lantern_full "$ROOT/truth_blind/LANTERN_full/LANTERN_ASSEMBLY.fasta"
add_method lantern_no_longitudinal "$ROOT/truth_blind/LANTERN_no_longitudinal/LANTERN_ASSEMBLY.fasta"
add_method lantern_no_consensus "$ROOT/truth_blind/LANTERN_no_consensus/LANTERN_ASSEMBLY.fasta"
[[ ${#METHODS[@]} -ge 6 ]]

minimap2 -I 4G -d "$ROOT/truth/toy_truth.mmi" "$ROOT/truth/combined_truth.fasta" > "$ROOT/logs/TRUTH_INDEX.stdout" 2> "$ROOT/logs/TRUTH_INDEX.log"
RECOVERY_ARGS=()
METHOD_ARGS=()
for spec in "${METHODS[@]}"; do
  name="${spec%%=*}"; fasta="${spec#*=}"; dir="$ROOT/evaluation/$name"; mkdir -p "$dir"
  sha256sum "$fasta" > "$dir/ASSEMBLY_SHA256.txt"
  minimap2 -x asm5 -c -N 100 -t "$THREADS" "$ROOT/truth/toy_truth.mmi" "$fasta" > "$dir/assembly_to_truth.paf" 2> "$dir/minimap2.log"
  python "$SCRIPTS/evaluate_gold_coverage.py" --paf "$dir/assembly_to_truth.paf" --truth-mapping "$ROOT/truth/truth_mapping.tsv" --assembly-fasta "$fasta" --min-identity 0.90 --min-alignment 500 --out "$dir" > "$dir/evaluate_console.log" 2>&1
  RECOVERY_ARGS+=(--recovery "$name=$dir/per_genome_recovery.tsv")
  METHOD_ARGS+=(--method "$name=$dir")
  rm -f "$dir/assembly_to_truth.paf"
done

python - "$ROOT" <<'PY'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1]); rows=list(csv.DictReader(open(root/'truth/genome_reference_manifest.tsv'),delimiter='\t')); paths=[]
for r in rows:
    if int(r['truth_bp'])<50000:continue
    p=root/'truth/references'/Path(r['reference_path']).name
    if not p.is_file():raise SystemExit(f'missing truth reference {p}')
    paths.append(str(p))
(root/'pseudo_novel/mash_reference_files.txt').write_text('\n'.join(paths)+'\n')
PY
if [[ $(wc -l < "$ROOT/pseudo_novel/mash_reference_files.txt") -ge 2 ]]; then
  mash sketch -k 21 -s 1000 -l "$ROOT/pseudo_novel/mash_reference_files.txt" -o "$ROOT/pseudo_novel/truth_genomes" > "$ROOT/logs/MASH_SKETCH.log" 2>&1
  mash dist "$ROOT/pseudo_novel/truth_genomes.msh" "$ROOT/pseudo_novel/truth_genomes.msh" > "$ROOT/pseudo_novel/truth_genome_distances.tsv"
  python "$SCRIPTS/build_pseudo_novel_audit.py" --manifest "$ROOT/truth/genome_reference_manifest.tsv" --mash "$ROOT/pseudo_novel/truth_genome_distances.tsv" "${RECOVERY_ARGS[@]}" --min-truth-bp 50000 --out "$ROOT/pseudo_novel" 2>&1 | tee "$ROOT/logs/PSEUDO_NOVEL.log"
  python - "$ROOT" <<'PY'
import csv,sys
from pathlib import Path
root=Path(sys.argv[1]);rows=list(csv.DictReader(open(root/'pseudo_novel/PSEUDO_NOVEL_TARGETS.tsv'),delimiter='\t'))
tiers=['exact_target_withheld','species_neighbour_withheld','genus_neighbour_withheld','deep_neighbour_withheld']
with open(root/'pseudo_novel/PSEUDO_NOVEL_SPLITS.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['genome_id','tier'],delimiter='\t');w.writeheader()
    for r in rows:
        for tier in tiers:
            if r[tier]=='true':w.writerow({'genome_id':r['genome_id'],'tier':tier})
PY
else
  printf 'genome_id\ttier\n' > "$ROOT/pseudo_novel/PSEUDO_NOVEL_SPLITS.tsv"
  printf '{"status":"INSUFFICIENT_TRUTH_GENOMES_FOR_MASH"}\n' > "$ROOT/pseudo_novel/PSEUDO_NOVEL_AUDIT.json"
fi

python "$SCRIPTS/summarize_benchmark.py" "${METHOD_ARGS[@]}" --manifest "$ROOT/METHOD_MANIFEST.tsv" --decision-config "$DECISION" --full-method lantern_full --no-longitudinal-method lantern_no_longitudinal --pseudo-splits "$ROOT/pseudo_novel/PSEUDO_NOVEL_SPLITS.tsv" --seed "$SEED" --bootstrap 5000 --out "$ROOT/decision" 2>&1 | tee "$ROOT/logs/FINAL_DECISION.log"

# Preserve compact reproducible evidence, not public raw reads or bulky intermediates.
mkdir -p "$ROOT/results/baseline_fastas"
for spec in "${METHODS[@]}"; do
  name="${spec%%=*}"; fasta="${spec#*=}"
  if [[ "$name" == megahit_* || "$name" == metaspades_* ]]; then cp "$fasta" "$ROOT/results/baseline_fastas/${name}.fasta"; fi
done
rm -rf "$ROOT/reads" "$ROOT/truth_download" "$ROOT/truth/references" "$ROOT/assemblies/megahit_single0/intermediate_contigs" "$ROOT/assemblies/megahit_single1/intermediate_contigs" "$ROOT/assemblies/megahit_longitudinal/intermediate_contigs" "$ROOT/assemblies/metaspades_longitudinal"
rm -f "$ROOT/truth/combined_truth.fasta" "$ROOT/truth/toy_truth.mmi" "$ROOT/pseudo_novel/truth_genomes.msh" "$ROOT/pseudo_novel/truth_genome_distances.tsv"
rm -f "$ROOT/truth_blind/evidence"/*.bam "$ROOT/truth_blind/evidence"/*.bai "$ROOT/truth_blind/all_vs_all.paf"
df -h > "$ROOT/provenance/DISK_AFTER.txt"
find "$ROOT" -type f ! -name SHA256SUMS_ALL.txt -print0 | sort -z | xargs -0 sha256sum > "$ROOT/SHA256SUMS_ALL.txt"
du -ah "$ROOT" | sort -h > "$ROOT/FILE_SIZES.txt"
cat "$ROOT/decision/FINAL_DECISION.md"
