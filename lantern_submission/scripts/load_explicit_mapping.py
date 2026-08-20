#!/usr/bin/env python3
"""Load and freeze an explicit sample-to-individual mapping.

This is the only production mapping path for the LANTERN CAMI III submission
package. It never infers relationships from sample order, sequence similarity,
Mash distances, taxonomy, assemblies, or benchmark results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

REQUIRED_COLUMNS = ("individual_id", "sample_id", "timepoint")
NATURAL_TOKEN = re.compile(r"(\d+)")


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(token) if token.isdigit() else token.casefold() for token in NATURAL_TOKEN.split(value))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_expected_samples(raw: str | None) -> set[str] | None:
    if raw is None:
        return None
    values = {part.strip() for part in raw.split(",") if part.strip()}
    if not values:
        raise ValueError("--expected-sample-ids was supplied but contained no IDs")
    return values


def read_manifest_samples(path: Path) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "sample_id" not in reader.fieldnames:
            raise ValueError("manifest must contain a sample_id column")
        samples: list[str] = []
        for line_number, row in enumerate(reader, 2):
            sample = (row.get("sample_id") or "").strip()
            if not sample:
                raise ValueError(f"empty sample_id in manifest at line {line_number}")
            samples.append(sample)
    if len(samples) != len(set(samples)):
        raise ValueError("duplicate sample_id in manifest")
    return set(samples)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--expected-individuals", type=int, required=True)
    parser.add_argument("--expected-timepoints", type=int, required=True)
    parser.add_argument("--expected-sample-ids", default=None,
                        help="Comma-separated exact sample ID set")
    parser.add_argument("--manifest", type=Path, default=None,
                        help="Optional TSV whose sample_id set must match the mapping")
    parser.add_argument("--require-consecutive-timepoints", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def write_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(REQUIRED_COLUMNS), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.expected_individuals <= 0 or args.expected_timepoints <= 0:
        raise SystemExit("expected counts must be positive")
    mapping = args.mapping
    if not mapping.is_file() or mapping.stat().st_size == 0:
        raise FileNotFoundError(mapping)

    with mapping.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("mapping has no header")
        missing_columns = sorted(set(REQUIRED_COLUMNS) - set(reader.fieldnames))
        if missing_columns:
            raise ValueError(f"mapping missing required columns: {missing_columns}")
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, 2):
            row = {column: (raw.get(column) or "").strip() for column in REQUIRED_COLUMNS}
            empty = [column for column, value in row.items() if not value]
            if empty:
                raise ValueError(f"empty fields at line {line_number}: {empty}")
            rows.append(row)

    if not rows:
        raise ValueError("mapping contains no data rows")

    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        duplicates = sorted({sample for sample in sample_ids if sample_ids.count(sample) > 1}, key=natural_key)
        raise ValueError(f"duplicate sample_id values: {duplicates}")

    individual_timepoints = [(row["individual_id"], row["timepoint"]) for row in rows]
    if len(individual_timepoints) != len(set(individual_timepoints)):
        duplicates = sorted({pair for pair in individual_timepoints if individual_timepoints.count(pair) > 1})
        raise ValueError(f"duplicate individual_id/timepoint values: {duplicates}")

    by_individual: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_individual[row["individual_id"]].append(row)

    if len(by_individual) != args.expected_individuals:
        raise ValueError(
            f"expected {args.expected_individuals} individuals, found {len(by_individual)}: "
            f"{sorted(by_individual, key=natural_key)}"
        )

    for individual, members in by_individual.items():
        if len(members) != args.expected_timepoints:
            raise ValueError(
                f"individual {individual!r} has {len(members)} timepoints; "
                f"expected {args.expected_timepoints}"
            )
        if args.require_consecutive_timepoints:
            try:
                numeric = sorted(int(member["timepoint"]) for member in members)
            except ValueError as exc:
                raise ValueError("consecutive-timepoint checking requires integer timepoint labels") from exc
            expected = list(range(1, args.expected_timepoints + 1))
            if numeric != expected:
                raise ValueError(
                    f"individual {individual!r} has timepoints {numeric}; expected {expected}"
                )

    mapping_samples = set(sample_ids)
    expected_samples = parse_expected_samples(args.expected_sample_ids)
    if expected_samples is not None and mapping_samples != expected_samples:
        raise ValueError(
            "mapping sample set mismatch: "
            f"missing={sorted(expected_samples - mapping_samples, key=natural_key)}, "
            f"extra={sorted(mapping_samples - expected_samples, key=natural_key)}"
        )
    if args.manifest is not None:
        manifest_samples = read_manifest_samples(args.manifest)
        if mapping_samples != manifest_samples:
            raise ValueError(
                "mapping/manifest sample set mismatch: "
                f"missing_from_manifest={sorted(mapping_samples - manifest_samples, key=natural_key)}, "
                f"extra_in_manifest={sorted(manifest_samples - mapping_samples, key=natural_key)}"
            )

    normalized_rows = sorted(
        rows,
        key=lambda row: (
            natural_key(row["individual_id"]),
            natural_key(row["timepoint"]),
            natural_key(row["sample_id"]),
        ),
    )
    args.out.mkdir(parents=True, exist_ok=True)
    normalized_path = args.out / "NORMALIZED_MAPPING.tsv"
    write_tsv(normalized_path, normalized_rows)

    groups = []
    for individual in sorted(by_individual, key=natural_key):
        members = sorted(
            by_individual[individual],
            key=lambda row: (natural_key(row["timepoint"]), natural_key(row["sample_id"])),
        )
        groups.append(
            {
                "individual_id": individual,
                "samples": [member["sample_id"] for member in members],
                "timepoints": [member["timepoint"] for member in members],
                "members": members,
            }
        )

    freeze = {
        "status": "EXPLICIT_MAPPING_FROZEN",
        "mapping_method": "explicit_metadata_only",
        "required_columns": list(REQUIRED_COLUMNS),
        "input_mapping_filename": mapping.name,
        "input_mapping_sha256": sha256_file(mapping),
        "normalized_mapping_sha256": sha256_file(normalized_path),
        "expected_individuals": args.expected_individuals,
        "expected_timepoints_per_individual": args.expected_timepoints,
        "n_individuals": len(groups),
        "n_samples": len(normalized_rows),
        "sample_ids": sorted(mapping_samples, key=natural_key),
        "groups": groups,
        "sequential_pairing_used": False,
        "similarity_pairing_used": False,
        "group_membership_source": "mapping TSV only",
        "manifest_sample_set_checked": args.manifest is not None,
        "claim_boundary": (
            "This file freezes explicit metadata group membership only. It contains no "
            "assembly score, CAMI gold result, inferred relationship, or performance claim."
        ),
    }
    freeze_path = args.out / "MAPPING_FREEZE.json"
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.out / "MAPPING_SHA256.txt").write_text(
        f"{freeze['input_mapping_sha256']}  {mapping.name}\n"
        f"{freeze['normalized_mapping_sha256']}  {normalized_path.name}\n"
        f"{sha256_file(freeze_path)}  {freeze_path.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(freeze, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
