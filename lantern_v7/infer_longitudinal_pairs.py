#!/usr/bin/env python3
"""Infer longitudinal Toy sample pairs from read-only Mash distances.

The pairing is determined from truth-blind read k-mer sketches.  No GSA, read-to-genome
mapping, taxonomic label, assembly score, or performance result is accepted.  For an
even number of samples the program solves the exact minimum-weight perfect matching by
dynamic programming and freezes one reserve pair using an input-only rule.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

SAMPLE_RE = re.compile(r"(?:^|[^0-9])S?(\d+)(?:[^0-9]|$)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--mash-dist", required=True)
    p.add_argument("--expected-samples", type=int, default=20)
    p.add_argument("--exclude-reserve-sample", action="append", type=int, default=[])
    p.add_argument("--out", required=True)
    return p.parse_args()


def sample_id(text: str) -> int:
    name = Path(text).name
    matches = SAMPLE_RE.findall(name)
    if not matches:
        raise ValueError(f"cannot parse sample ID from {text!r}")
    return int(matches[-1])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = Path(a.mash_dist)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    distances: dict[tuple[int, int], float] = {}
    shared: dict[tuple[int, int], str] = {}
    samples: set[int] = set()
    with path.open("rt", encoding="utf-8", errors="strict") as fh:
        for line_number, raw in enumerate(fh, 1):
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            left = sample_id(cols[0])
            right = sample_id(cols[1])
            try:
                distance = float(cols[2])
            except ValueError as exc:
                raise ValueError(f"invalid Mash distance at line {line_number}") from exc
            if not math.isfinite(distance) or distance < 0:
                raise ValueError(f"invalid distance {distance} at line {line_number}")
            samples.update((left, right))
            if left == right:
                continue
            key = tuple(sorted((left, right)))
            if key in distances and abs(distances[key] - distance) > 1e-12:
                raise ValueError(f"asymmetric/conflicting Mash distance for {key}")
            distances[key] = distance
            if len(cols) >= 5:
                shared[key] = cols[4]

    ordered = sorted(samples)
    if len(ordered) != a.expected_samples:
        raise SystemExit(f"expected {a.expected_samples} samples, found {ordered}")
    if len(ordered) % 2:
        raise SystemExit("sample count must be even")
    for i, left in enumerate(ordered):
        for right in ordered[i + 1 :]:
            if (left, right) not in distances:
                raise SystemExit(f"missing pairwise distance {(left, right)}")

    index = {sample: i for i, sample in enumerate(ordered)}
    full_mask = (1 << len(ordered)) - 1

    @lru_cache(maxsize=None)
    def solve(mask: int) -> tuple[float, tuple[tuple[int, int], ...]]:
        if mask == 0:
            return 0.0, ()
        first_index = (mask & -mask).bit_length() - 1
        first = ordered[first_index]
        remaining = mask ^ (1 << first_index)
        best_cost = math.inf
        best_pairs: tuple[tuple[int, int], ...] | None = None
        candidates = [i for i in range(first_index + 1, len(ordered)) if remaining & (1 << i)]
        for other_index in candidates:
            other = ordered[other_index]
            next_mask = remaining ^ (1 << other_index)
            sub_cost, sub_pairs = solve(next_mask)
            pair = tuple(sorted((first, other)))
            cost = distances[pair] + sub_cost
            pairs = tuple(sorted((pair,) + sub_pairs))
            if cost < best_cost - 1e-15 or (
                abs(cost - best_cost) <= 1e-15 and (best_pairs is None or pairs < best_pairs)
            ):
                best_cost = cost
                best_pairs = pairs
        if best_pairs is None:
            raise RuntimeError("no perfect matching")
        return best_cost, best_pairs

    total_cost, inferred_pairs = solve(full_mask)

    nearest_rows: list[dict[str, Any]] = []
    for sample in ordered:
        neighbors = sorted(
            (distances[tuple(sorted((sample, other)))], other)
            for other in ordered if other != sample
        )
        nearest = neighbors[0]
        second = neighbors[1]
        nearest_rows.append(
            {
                "sample": sample,
                "nearest_sample": nearest[1],
                "nearest_distance": nearest[0],
                "second_nearest_sample": second[1],
                "second_nearest_distance": second[0],
                "distance_gap": second[0] - nearest[0],
                "nearest_is_reciprocal": str(
                    min(
                        (distances[tuple(sorted((nearest[1], other)))], other)
                        for other in ordered if other != nearest[1]
                    )[1] == sample
                ).lower(),
            }
        )

    pair_rows: list[dict[str, Any]] = []
    for left, right in inferred_pairs:
        left_nearest = next(row for row in nearest_rows if row["sample"] == left)
        right_nearest = next(row for row in nearest_rows if row["sample"] == right)
        pair_rows.append(
            {
                "sample_a": left,
                "sample_b": right,
                "mash_distance": distances[(left, right)],
                "shared_hashes": shared.get((left, right), ""),
                "reciprocal_nearest_neighbors": str(
                    left_nearest["nearest_sample"] == right
                    and right_nearest["nearest_sample"] == left
                ).lower(),
                "minimum_distance_gap": min(
                    float(left_nearest["distance_gap"]),
                    float(right_nearest["distance_gap"]),
                ),
            }
        )

    excluded = set(a.exclude_reserve_sample)
    eligible = [
        row for row in pair_rows
        if int(row["sample_a"]) not in excluded and int(row["sample_b"]) not in excluded
    ]
    if not eligible:
        raise SystemExit("no eligible reserve pair after exclusions")
    eligible_sorted = sorted(eligible, key=lambda row: (float(row["mash_distance"]), int(row["sample_a"]), int(row["sample_b"])))
    reserve = eligible_sorted[(len(eligible_sorted) - 1) // 2]

    with (out / "INFERRED_LONGITUDINAL_PAIRS.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(pair_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)
    with (out / "SAMPLE_NEAREST_NEIGHBORS.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(nearest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(nearest_rows)

    result = {
        "status": "TRUTH_BLIND_PAIRING_FROZEN",
        "samples": ordered,
        "n_pairs": len(pair_rows),
        "minimum_weight_matching_total_distance": total_cost,
        "pairs": [[int(row["sample_a"]), int(row["sample_b"])] for row in pair_rows],
        "reciprocal_nearest_pair_count": sum(row["reciprocal_nearest_neighbors"] == "true" for row in pair_rows),
        "reserve_excluded_samples": sorted(excluded),
        "reserved_untouched_pair": [int(reserve["sample_a"]), int(reserve["sample_b"])],
        "reserve_rule": "median Mash distance among matching pairs that avoid prespecified previously exposed samples; ties lexicographic",
        "input_sha256": sha256(path),
        "truth_accessed": False,
        "forbidden_inputs": ["GSA", "read-to-genome mapping", "taxonomy", "assembly performance", "CAMI gold score"],
        "boundary": (
            "Sample pairing is inferred solely from truth-blind short-read k-mer sketches. "
            "It is an operational pairing for method development and is not claimed as an official CAMI assignment."
        ),
    }
    (out / "TRUTH_BLIND_PAIRING_FREEZE.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
