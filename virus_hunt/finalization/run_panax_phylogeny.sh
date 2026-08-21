#!/usr/bin/env bash
# Two-analysis sensitivity audit of homologous PF00680 cores.
set -Eeuo pipefail

CORES="${1:?PF00680 core FASTA required}"
MANIFEST="${2:?curated reference manifest required}"
CURRENT_PANEL_MANIFEST="${3:?current nr panel manifest required}"
OUT="${4:?output directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"
if find "$OUT" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: $OUT" >&2
  exit 2
fi

for command in python mafft trimal; do
  command -v "$command" >/dev/null || { echo "missing required command: $command" >&2; exit 2; }
done
if command -v iqtree2 >/dev/null; then
  IQTREE=iqtree2
elif command -v iqtree >/dev/null; then
  IQTREE=iqtree
else
  echo "missing IQ-TREE 2 executable" >&2
  exit 2
fi

date -u +%FT%TZ > "$OUT/STARTED_UTC.txt"
cp "$CORES" "$OUT/RDRP_CORES.faa"
cp "$MANIFEST" "$OUT/CURATED_REFERENCE_MANIFEST.tsv"
cp "$CURRENT_PANEL_MANIFEST" "$OUT/CURRENT_NR_REFERENCE_CONTRACT.tsv"
{
  printf 'mafft\t'; mafft --version 2>&1 | sed -n '1p'
  printf 'trimal\t'; trimal --version 2>&1 | sed -n '1p'
  printf 'trimal_source_repository\t%s\n' "${TRIMAL_SOURCE_REPOSITORY:-unrecorded}"
  printf 'trimal_source_commit\t%s\n' "${TRIMAL_SOURCE_COMMIT:-unrecorded}"
  printf 'iqtree\t'; "$IQTREE" -version 2>&1 | sed -n '1p'
  printf 'python\t'; python --version 2>&1
} > "$OUT/TOOL_VERSIONS.txt"

printf 'mafft --localpair --maxiterate 1000 --reorder RDRP_CORES.faa\n' > "$OUT/COMMANDS.txt"
mafft --localpair --maxiterate 1000 --reorder "$OUT/RDRP_CORES.faa" \
  > "$OUT/RDRP_CORES.untrimmed.aln.faa" 2> "$OUT/MAFFT.stderr.txt"
printf 'trimal -in RDRP_CORES.untrimmed.aln.faa -out RDRP_CORES.trimmed.aln.faa -automated1\n' >> "$OUT/COMMANDS.txt"
trimal -in "$OUT/RDRP_CORES.untrimmed.aln.faa" -out "$OUT/RDRP_CORES.trimmed.aln.faa" -automated1 \
  > "$OUT/TRIMAL.stdout.txt" 2> "$OUT/TRIMAL.stderr.txt"

for label in untrimmed trimmed; do
  alignment="$OUT/RDRP_CORES.${label}.aln.faa"
  prefix="$OUT/${label}"
  printf '%s -s %s -m MFP -bb 1000 -alrt 1000 -seed 20260821 -nt 2 -pre %s\n' \
    "$IQTREE" "$(basename "$alignment")" "$label" >> "$OUT/COMMANDS.txt"
  "$IQTREE" -s "$alignment" -m MFP -bb 1000 -alrt 1000 -seed 20260821 -nt 2 \
    -pre "$prefix" -redo > "$OUT/${label}.console.txt" 2>&1
  [[ -s "$prefix.treefile" && -s "$prefix.iqtree" && -s "$prefix.log" ]]
done

date -u +%FT%TZ > "$OUT/FINISHED_UTC.txt"
python "$SCRIPT_DIR/audit_panax_phylogeny.py" \
  --cores "$OUT/RDRP_CORES.faa" \
  --manifest "$OUT/CURATED_REFERENCE_MANIFEST.tsv" \
  --current-panel-manifest "$OUT/CURRENT_NR_REFERENCE_CONTRACT.tsv" \
  --untrimmed-alignment "$OUT/RDRP_CORES.untrimmed.aln.faa" \
  --trimmed-alignment "$OUT/RDRP_CORES.trimmed.aln.faa" \
  --untrimmed-tree "$OUT/untrimmed.treefile" --trimmed-tree "$OUT/trimmed.treefile" \
  --out "$OUT"
