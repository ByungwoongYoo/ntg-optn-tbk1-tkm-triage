#!/usr/bin/env python3
"""Combine validated per-individual FASTAs with globally unique prefixed IDs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ALLOWED = re.compile(r"^[ATCGN]+$")
SAFE_ID = re.compile(r"[^A-Za-z0-9_.|:-]+")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True,
                   help="TSV with individual_id and fasta columns")
    p.add_argument("--out-fasta", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    return p.parse_args()


def read_fasta(path: Path):
    name = None
    parts: list[str] = []
    with path.open("rt", encoding="ascii", errors="strict") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier at {path}:{line_number}")
                parts = []
            else:
                if name is None:
                    raise ValueError(f"sequence before first header at {path}:{line_number}")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def write_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for index in range(0, len(sequence), width):
        handle.write(sequence[index:index + width] + "\n")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    a = parse_args()
    if not a.manifest.is_file():
        raise FileNotFoundError(a.manifest)
    with a.manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"individual_id", "fasta"}.issubset(reader.fieldnames):
            raise ValueError("manifest must contain individual_id and fasta columns")
        rows = [
            {"individual_id": (row.get("individual_id") or "").strip(),
             "fasta": (row.get("fasta") or "").strip()}
            for row in reader
        ]
    if not rows or any(not row["individual_id"] or not row["fasta"] for row in rows):
        raise ValueError("empty or incomplete FASTA manifest")
    if len({row["individual_id"] for row in rows}) != len(rows):
        raise ValueError("duplicate individual_id in FASTA manifest")

    rows.sort(key=lambda row: int(row["individual_id"]) if row["individual_id"].isdigit() else row["individual_id"])
    a.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    seen_ids: set[str] = set()
    sequence_hashes: set[str] = set()
    details = []
    total_records = total_bp = 0
    with a.out_fasta.open("wt", encoding="ascii", newline="\n") as output:
        for row in rows:
            path = Path(row["fasta"])
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
            n_records = n_bp = 0
            for original_id, sequence in read_fasta(path):
                if not sequence or not ALLOWED.fullmatch(sequence):
                    raise ValueError(f"invalid sequence in {path}: {original_id}")
                clean = SAFE_ID.sub("_", original_id)
                output_id = f"I{row['individual_id']}|{clean}"
                if output_id in seen_ids:
                    raise ValueError(f"duplicate output identifier: {output_id}")
                seen_ids.add(output_id)
                digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
                if digest in sequence_hashes:
                    continue
                sequence_hashes.add(digest)
                write_record(output, output_id, sequence)
                n_records += 1
                n_bp += len(sequence)
                total_records += 1
                total_bp += len(sequence)
            details.append({
                "individual_id": row["individual_id"],
                "input_fasta": str(path),
                "input_sha256": sha256_file(path),
                "records_emitted": n_records,
                "bases_emitted": n_bp,
            })
    if total_records == 0:
        raise ValueError("combined FASTA is empty")
    summary = {
        "status": "PASS",
        "output_fasta": str(a.out_fasta),
        "output_sha256": sha256_file(a.out_fasta),
        "records": total_records,
        "total_bp": total_bp,
        "duplicate_output_ids": 0,
        "alphabet": "ATCGN",
        "individuals": details,
    }
    a.out_json.parent.mkdir(parents=True, exist_ok=True)
    a.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
