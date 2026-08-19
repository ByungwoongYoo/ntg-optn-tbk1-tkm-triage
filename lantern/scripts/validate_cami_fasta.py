#!/usr/bin/env python3
"""CAMI FASTA integrity and summary validator.

Two validation modes are intentionally separated:

- ``submission``: CAMI assembly-output syntax, requiring uppercase A/T/C/G/N.
- ``reference``: official source/gold-standard inputs, allowing upper/lowercase DNA
  IUPAC symbols while reporting lowercase and non-ATCGN ambiguity explicitly.

Passing validates syntax and file integrity only. It does not imply biological accuracy,
assembly quality, novelty, or CAMI acceptance.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import TextIO

HEADER_RE = re.compile(r"^>[A-Za-z0-9\s\[\]_:;,\.\|\-]+$")
SUBMISSION_SEQ_RE = re.compile(r"^[ATCGN]+$")
REFERENCE_SEQ_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$")
CANONICAL_UPPER = set("ATCGN")
CANONICAL_BOTH_CASES = set("ATCGNatcgn")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("fasta", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--mode",
        choices=("submission", "reference"),
        default="submission",
        help="submission=strict CAMI output; reference=official reference/GSA input",
    )
    return p.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return path.open("rt", encoding="ascii", errors="strict")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    target = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    raise AssertionError("unreachable")


def main() -> None:
    args = parse_args()
    path = args.fasta
    if not path.is_file():
        raise FileNotFoundError(path)

    sequence_re = SUBMISSION_SEQ_RE if args.mode == "submission" else REFERENCE_SEQ_RE
    seen: set[str] = set()
    lengths: list[int] = []
    errors: list[str] = []
    current_id: str | None = None
    current_length = 0
    sequence_lines = 0
    lowercase_bases = 0
    iupac_ambiguity_bases = 0
    canonical_n_bases = 0

    with open_text(path) as fh:
        for line_number, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if current_id is not None:
                    lengths.append(current_length)
                if not HEADER_RE.fullmatch(line):
                    errors.append(f"line {line_number}: invalid FASTA header characters")
                token = line[1:].split(None, 1)[0] if len(line) > 1 else ""
                if not token:
                    errors.append(f"line {line_number}: empty sequence identifier")
                elif token in seen:
                    errors.append(f"line {line_number}: duplicate sequence identifier {token}")
                else:
                    seen.add(token)
                current_id = token or None
                current_length = 0
                continue

            if current_id is None:
                errors.append(f"line {line_number}: sequence before first header")
                continue
            if not line:
                errors.append(f"line {line_number}: blank sequence line")
                continue
            if not sequence_re.fullmatch(line):
                requirement = (
                    "uppercase ATCGN"
                    if args.mode == "submission"
                    else "upper/lowercase DNA IUPAC symbols"
                )
                errors.append(f"line {line_number}: sequence violates {requirement}")

            lowercase_bases += sum(ch.islower() for ch in line)
            iupac_ambiguity_bases += sum(ch not in CANONICAL_BOTH_CASES for ch in line)
            canonical_n_bases += sum(ch in {"N", "n"} for ch in line)
            current_length += len(line)
            sequence_lines += 1

    if current_id is not None:
        lengths.append(current_length)
    if not lengths:
        errors.append("no FASTA records")
    if any(length == 0 for length in lengths):
        errors.append("one or more zero-length records")

    summary = {
        "path": str(path),
        "mode": args.mode,
        "bytes": path.stat().st_size,
        "sha256": file_sha256(path),
        "records": len(lengths),
        "sequence_lines": sequence_lines,
        "total_bases": sum(lengths),
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "mean_length": (sum(lengths) / len(lengths)) if lengths else 0,
        "n50": n50(lengths),
        "lowercase_bases": lowercase_bases,
        "iupac_ambiguity_bases_excluding_n": iupac_ambiguity_bases,
        "n_bases": canonical_n_bases,
        "submission_compatible": (
            not errors and lowercase_bases == 0 and iupac_ambiguity_bases == 0
        ),
        "errors": errors,
        "valid": not errors,
        "claim_boundary": (
            "Syntax/integrity validation only; not biological assembly evaluation. "
            "Reference-mode validity is not CAMI submission-output compatibility."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
