#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EVIDENCE_ROOT="${1:-$ROOT}"
cd "$ROOT"

echo "[1/3] Verify all raw files, branch coverage, hashes, witness, and stored proof logs"
python3 verification/verify_raw_evidence.py "$EVIDENCE_ROOT"

echo "[2/3] Compile all three exact-search engines"
mkdir -p verification/build
g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  source/b2_gap1_v_split_v10.cpp \
  -o verification/build/b2_v10
g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  source/b2_gap1_u_split_v8.cpp \
  -o verification/build/b2_v8
g++ -std=c++20 -O3 -march=native -pthread -DNDEBUG \
  source/b2_gap_g_u_split_v9.cpp \
  -o verification/build/b2_v9

echo "[3/3] Run small exact replay controls"
verification/build/b2_v10 4 87 120 4 \
  verification/build/v10_t4_u87.json verification/build/v10_t4_u87.tsv
verification/build/b2_v8 45 120 4 \
  verification/build/v8_t45.json verification/build/v8_t45.tsv
verification/build/b2_v9 3 35 120 4 \
  verification/build/v9_g3_t35.json verification/build/v9_g3_t35.tsv

python3 - <<'PY'
import json
from pathlib import Path
p=Path("verification/build")
for name in ["v10_t4_u87.json","v8_t45.json","v9_g3_t35.json"]:
    d=json.loads((p/name).read_text())
    assert d["completed_exhaustively"] is True
    assert d["timed_out"] is False
    assert d["witness_found"] is False
print("CLEAN_ROOM_QUICK_REPLAY_PASS")
PY
