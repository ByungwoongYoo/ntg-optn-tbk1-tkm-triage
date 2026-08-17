#!/usr/bin/env python3
"""Independent small-object checks for the B2[2] Z/100Z resolution record.

This verifier checks the explicit 13-set and the elementary mod-10 support
classification used in the unit-difference lemma. The exhaustive upper-bound
search is reproduced by the C++ sources and GitHub Actions run IDs recorded in
B2_Z100_FINAL_RESULTS_20260817.json.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = json.loads((ROOT / "B2_Z100_FINAL_RESULTS_20260817.json").read_text())

assert RESULTS["answer"] == "NO"
assert RESULTS["maximum_size"] == 13
assert RESULTS["upper_bound_argument"]["unit_difference_lemma"] is True
assert RESULTS["upper_bound_argument"]["third_branches_exhausted"] == 44
assert RESULTS["upper_bound_argument"]["all_no_witness"] is True
assert RESULTS["summary_artifacts"]["g1_t2_t4"]["all_hard_exhausted_no_witness"] is True

A = RESULTS["lower_bound_witness"]
counts = {d: 0 for d in range(1, 100)}
for a in A:
    for b in A:
        if a != b:
            counts[(a - b) % 100] += 1
assert len(A) == 13 and len(set(A)) == 13
assert max(counts.values()) == 2

# Exhaustively audit the residue-support classification in the unit-difference lemma.
unit_mod_10 = {1, 3, 7, 9}
checked = 0
for mask in range(1, 1 << 10):
    support = [r for r in range(10) if (mask >> r) & 1]
    if all((a - b) % 10 not in unit_mod_10 for a in support for b in support):
        same_parity = len({r % 2 for r in support}) == 1
        mixed_pair = (
            not same_parity
            and len(support) == 2
            and (support[0] - support[1]) % 10 == 5
        )
        assert same_parity or mixed_pair, support
        checked += 1
assert checked == 67
assert 14 * 13 == 182
assert 49 * 2 == 98 < 182
assert 19 * 2 == 38 < 182

print("VERIFIED: valid 13-set; unit-difference lemma support audit; recorded 44-branch exact upper-bound result gives maximum 13.")
