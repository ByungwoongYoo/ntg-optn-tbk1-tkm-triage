#!/usr/bin/env python3
"""Attach proper-fragment depth metrics to a Panax fragment audit table.

The input depth files must have been generated from BAMs containing only
qualifying paired fragments selected by ``audit_panax_fragments.py``.  This
script validates every reference coordinate before calculating breadth,
depth, and zero-run summaries.  It defines no biological detection threshold.
"""
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                lengths[name] = length
            name = line[1:].split()[0]
            if not name or name in lengths:
                raise ValueError(f"invalid or duplicate FASTA identifier: {name!r}")
            length = 0
        elif name is None:
            raise ValueError("sequence before FASTA header")
        else:
            length += len(line)
    if name is not None:
        lengths[name] = length
    if not lengths or any(length < 1 for length in lengths.values()):
        raise ValueError("empty reference FASTA record")
    return lengths


def read_depth(path: Path, reference: str, length: int) -> list[int]:
    if not path.is_file():
        raise ValueError(f"missing depth file: {path}")
    depths: list[int] = []
    with path.open() as handle:
        for expected_position, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"malformed depth row in {path}: {raw.rstrip()!r}")
            observed_reference, position, depth = fields
            if observed_reference != reference or int(position) != expected_position:
                raise ValueError(
                    f"non-contiguous depth coordinates in {path}: expected "
                    f"{reference}:{expected_position}, observed {observed_reference}:{position}"
                )
            value = int(depth)
            if value < 0:
                raise ValueError(f"negative depth in {path}: {value}")
            depths.append(value)
    if len(depths) != length:
        raise ValueError(
            f"depth/reference length mismatch for {reference}: "
            f"depth={len(depths)}, reference={length}, file={path}"
        )
    return depths


def longest_zero_run(values: list[int]) -> int:
    longest = current = 0
    for value in values:
        if value == 0:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def summarize(depths: list[int], prefix: str) -> dict[str, str | int]:
    covered = [index for index, depth in enumerate(depths) if depth > 0]
    internal = depths[covered[0] : covered[-1] + 1] if covered else []
    total = len(depths)
    return {
        f"{prefix}_breadth_1x": f"{sum(value >= 1 for value in depths) / total:.6f}",
        f"{prefix}_breadth_5x": f"{sum(value >= 5 for value in depths) / total:.6f}",
        f"{prefix}_breadth_10x": f"{sum(value >= 10 for value in depths) / total:.6f}",
        f"{prefix}_mean_depth": f"{sum(depths) / total:.6f}",
        f"{prefix}_median_depth": f"{statistics.median(depths):.3f}",
        f"{prefix}_max_zero_run": longest_zero_run(depths),
        f"{prefix}_max_internal_zero_run": longest_zero_run(internal),
    }


def augment(
    rows: list[dict[str, str]], lengths: dict[str, int], depth_dir: Path
) -> list[dict[str, str | int]]:
    if not rows:
        raise ValueError("empty fragment metrics table")
    references = [row.get("reference", "") for row in rows]
    if len(references) != len(set(references)) or set(references) != set(lengths):
        raise ValueError(
            f"fragment/reference set mismatch: metrics={references}, FASTA={sorted(lengths)}"
        )
    output: list[dict[str, str | int]] = []
    for row in rows:
        reference = row["reference"]
        combined: dict[str, str | int] = dict(row)
        for label, prefix in (
            ("preduplicate", "proper_preduplicate"),
            ("nonduplicate", "proper_nonduplicate"),
        ):
            depths = read_depth(
                depth_dir / f"{reference}.{label}.depth.tsv",
                reference,
                lengths[reference],
            )
            combined.update(summarize(depths, prefix))
        output.append(combined)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--depth-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    with args.metrics.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = augment(rows, read_fasta_lengths(args.reference), args.depth_dir)
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite output: {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
