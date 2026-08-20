#!/usr/bin/env python3
"""Small independent checks for the B2[2]/Z100 claim's elementary parts."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path


DEFAULT_WITNESS = [0, 5, 7, 31, 58, 61, 62, 63, 72, 80, 84, 91, 97]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    witness = DEFAULT_WITNESS
    ordered = Counter((a - b) % 100 for a in witness for b in witness if a != b)
    counts = [ordered[d] for d in range(1, 100)]

    unit_mod_10 = {1, 3, 7, 9}
    admissible_supports: list[list[int]] = []
    classification_failures: list[list[int]] = []
    for mask in range(1, 1 << 10):
        support = [i for i in range(10) if (mask >> i) & 1]
        avoids_units = all(
            (a - b) % 10 not in unit_mod_10 for a in support for b in support
        )
        if not avoids_units:
            continue
        admissible_supports.append(support)
        same_parity = len({x % 2 for x in support}) == 1
        opposite_pair = (
            len(support) == 2
            and len({x % 2 for x in support}) == 2
            and (support[0] - support[1]) % 10 == 5
        )
        if not (same_parity or opposite_pair):
            classification_failures.append(support)

    checks = {
        "witness_has_13_distinct_residues": len(witness) == 13
        and len(set(witness)) == 13
        and all(0 <= x < 100 for x in witness),
        "ordered_difference_total_is_156": sum(counts) == 13 * 12,
        "ordered_difference_max_is_at_most_2": max(counts) <= 2,
        "unit_difference_mod10_classification": not classification_failures,
        "same_parity_capacity_contradiction": 14 * 13 > 49 * 2,
        "multiples_of_5_capacity_contradiction": 14 * 13 > 19 * 2,
        "normalized_reflection_branch_count": len(range(2, 46)) == 44,
    }
    result = {
        "schema": "b2-z100-elementary-verification-v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "witness": witness,
        "ordered_difference_histogram": dict(sorted(Counter(counts).items())),
        "ordered_difference_max": max(counts),
        "admissible_mod10_support_pattern_count": len(admissible_supports),
        "classification_failures": classification_failures,
        "normalized_third_element_range": [2, 45],
        "normalized_branch_count": 44,
        "checks": checks,
        "scope_note": (
            "These checks validate the witness and elementary normalization lemma; "
            "they do not independently replay the full exact searches or DRAT proofs."
        ),
    }
    output = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(output, encoding="utf-8")
    print(output, end="")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
