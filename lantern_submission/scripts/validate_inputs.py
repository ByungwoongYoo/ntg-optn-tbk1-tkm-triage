#!/usr/bin/env python3
"""Validate a private active-CAMI read manifest against an explicit mapping freeze."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

REQUIRED = (
    "sample_id",
    "short_r1",
    "short_r2",
    "long_reads",
    "short_r1_sha256",
    "short_r2_sha256",
    "long_reads_sha256",
)
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_gzip(path: Path) -> None:
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mapping-freeze", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--skip-gzip-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.manifest.is_file() or args.manifest.stat().st_size == 0:
        raise FileNotFoundError(args.manifest)
    if not args.mapping_freeze.is_file():
        raise FileNotFoundError(args.mapping_freeze)

    mapping = json.loads(args.mapping_freeze.read_text(encoding="utf-8"))
    if mapping.get("status") != "EXPLICIT_MAPPING_FROZEN":
        raise ValueError("mapping freeze status is not EXPLICIT_MAPPING_FROZEN")
    if mapping.get("mapping_method") != "explicit_metadata_only":
        raise ValueError("mapping freeze is not explicit-metadata-only")
    expected_samples = set(str(value) for value in mapping.get("sample_ids", []))
    if not expected_samples:
        raise ValueError("mapping freeze contains no sample IDs")

    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("manifest has no header")
        missing = sorted(set(REQUIRED) - set(reader.fieldnames))
        if missing:
            raise ValueError(f"manifest missing columns: {missing}")
        rows = []
        for line_number, raw in enumerate(reader, 2):
            row = {key: (raw.get(key) or "").strip() for key in REQUIRED}
            empty = [key for key, value in row.items() if not value]
            if empty:
                raise ValueError(f"empty manifest fields at line {line_number}: {empty}")
            rows.append(row)
    if not rows:
        raise ValueError("manifest contains no rows")

    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("duplicate sample_id in manifest")
    actual_samples = set(sample_ids)
    if actual_samples != expected_samples:
        raise ValueError(
            "manifest/mapping sample mismatch: "
            f"missing={sorted(expected_samples - actual_samples)}, "
            f"extra={sorted(actual_samples - expected_samples)}"
        )

    validated = []
    for row in sorted(rows, key=lambda item: int(item["sample_id"]) if item["sample_id"].isdigit() else item["sample_id"]):
        record = {"sample_id": row["sample_id"], "files": {}}
        for field, checksum_field in (
            ("short_r1", "short_r1_sha256"),
            ("short_r2", "short_r2_sha256"),
            ("long_reads", "long_reads_sha256"),
        ):
            path = Path(row[field]).expanduser()
            expected = row[checksum_field].lower()
            if not SHA256_RE.fullmatch(expected):
                raise ValueError(f"invalid SHA-256 for sample {row['sample_id']} field {field}")
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
            actual = sha256_file(path)
            if actual != expected:
                raise ValueError(
                    f"SHA-256 mismatch for sample {row['sample_id']} field {field}: "
                    f"expected {expected}, got {actual}"
                )
            if not args.skip_gzip_test:
                check_gzip(path)
            record["files"][field] = {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": actual,
                "gzip_integrity": "SKIPPED" if args.skip_gzip_test else "PASS",
            }
        validated.append(record)

    args.out.mkdir(parents=True, exist_ok=True)
    freeze = {
        "status": "ACTIVE_INPUT_MANIFEST_FROZEN",
        "mapping_freeze_sha256": sha256_file(args.mapping_freeze),
        "input_manifest_filename": args.manifest.name,
        "input_manifest_sha256": sha256_file(args.manifest),
        "n_samples": len(validated),
        "sample_ids": sorted(actual_samples, key=lambda value: int(value) if value.isdigit() else value),
        "records": validated,
        "all_checksums_match": True,
        "all_files_nonempty": True,
        "gzip_integrity_checked": not args.skip_gzip_test,
        "storage_boundary": "This freeze must remain on private local/HPC storage for restricted active data.",
    }
    output = args.out / "INPUT_MANIFEST_FREEZE.json"
    output.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.out / "MANIFEST.sha256").write_text(
        f"{freeze['mapping_freeze_sha256']}  {args.mapping_freeze.name}\n"
        f"{freeze['input_manifest_sha256']}  {args.manifest.name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(freeze, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
