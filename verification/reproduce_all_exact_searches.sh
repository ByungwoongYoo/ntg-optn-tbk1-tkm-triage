#!/usr/bin/env bash
# Full replay of every final exact-search branch. Expensive by design.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/verification/full_replay_output}"
THREADS="${THREADS:-$(nproc)}"
LIMIT="${LIMIT:-20000}"
mkdir -p "$OUT/bin" "$OUT/g1_v10" "$OUT/g1_v8" "$OUT/g23"

g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  "$ROOT/source/b2_gap1_v_split_v10.cpp" -o "$OUT/bin/v10"
g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  "$ROOT/source/b2_gap1_u_split_v8.cpp" -o "$OUT/bin/v8"
g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  "$ROOT/source/b2_gap_g_u_split_v9.cpp" -o "$OUT/bin/v9"

for t in 2 3 4; do
  for u in $(seq $((t+1)) $((91-t))); do
    "$OUT/bin/v10" "$t" "$u" "$LIMIT" "$THREADS" \
      "$OUT/g1_v10/t${t}_u${u}.json" "$OUT/g1_v10/t${t}_u${u}.tsv"
  done
done

for t in $(seq 5 45); do
  "$OUT/bin/v8" "$t" "$LIMIT" "$THREADS" \
    "$OUT/g1_v8/t${t}.json" "$OUT/g1_v8/t${t}.tsv"
done

for t in $(seq 4 40); do
  "$OUT/bin/v9" 2 "$t" "$LIMIT" "$THREADS" \
    "$OUT/g23/g2_t${t}.json" "$OUT/g23/g2_t${t}.tsv"
done
for t in $(seq 6 35); do
  "$OUT/bin/v9" 3 "$t" "$LIMIT" "$THREADS" \
    "$OUT/g23/g3_t${t}.json" "$OUT/g23/g3_t${t}.tsv"
done

python3 - "$OUT" <<'PY'
import json, sys
from pathlib import Path
root=Path(sys.argv[1])
files=list(root.rglob("*.json"))
assert len(files)==255+41+67, len(files)
for p in files:
    d=json.loads(p.read_text())
    assert d["completed_exhaustively"] is True, p
    assert d["timed_out"] is False, p
    assert d["witness_found"] is False, p
print(f"FULL_EXACT_REPLAY_PASS files={len(files)}")
PY
