#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EVIDENCE_ROOT="${1:-$ROOT}"
BUILD="$ROOT/verification/drat_build"
rm -rf "$BUILD"
mkdir -p "$BUILD"

git clone https://github.com/marijnheule/drat-trim.git "$BUILD/drat-trim"
git -C "$BUILD/drat-trim" checkout 2e3b2dc0ecf938addbd779d42877b6ed69d9a985
make -C "$BUILD/drat-trim" -j"$(nproc)"

: > "$BUILD/DRAT_REPLAY.log"
for g in 4 5 6 7; do
  d=$(find "$EVIDENCE_ROOT/evidence/g4_g7_drat/g${g}" -mindepth 1 -maxdepth 1 -type d -print -quit)
  cnf="$d/gap${g}.cnf"
  proof="$d/gap${g}.drat"
  "$BUILD/drat-trim/drat-trim" "$cnf" "$proof" > "$BUILD/g${g}.log" 2>&1
  grep -q '^s VERIFIED' "$BUILD/g${g}.log"
  cat "$BUILD/g${g}.log" >> "$BUILD/DRAT_REPLAY.log"
  sha256sum "$cnf" "$proof" >> "$BUILD/DRAT_REPLAY.log"
done
echo VERIFIED_ALL_G4_G7_DRAT | tee -a "$BUILD/DRAT_REPLAY.log"
