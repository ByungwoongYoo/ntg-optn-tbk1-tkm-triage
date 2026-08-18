#!/usr/bin/env python3
"""Stream-filter FASTA records by length and nucleotide alphabet.

Reads stdin and writes stdout. A compact stderr JSON summary makes every scan
auditable without retaining the full Logan assembly.
"""
from __future__ import annotations

import argparse
import json
import sys


def emit(header: str | None, seq_parts: list[str], min_len: int, max_len: int) -> tuple[int, int]:
    if header is None:
        return 0, 0
    seq = "".join(seq_parts).replace(" ", "").upper()
    length = len(seq)
    if min_len <= length <= max_len and seq and sum(c in "ACGTUN" for c in seq) / length >= 0.95:
        sys.stdout.write(header + "\n")
        for i in range(0, length, 80):
            sys.stdout.write(seq[i : i + 80] + "\n")
        return 1, length
    return 0, 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-len", type=int, default=600)
    parser.add_argument("--max-len", type=int, default=100000)
    args = parser.parse_args()

    header: str | None = None
    seq_parts: list[str] = []
    seen = kept = kept_bases = 0
    for raw in sys.stdin:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seen += 1
                n, b = emit(header, seq_parts, args.min_len, args.max_len)
                kept += n
                kept_bases += b
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        seen += 1
        n, b = emit(header, seq_parts, args.min_len, args.max_len)
        kept += n
        kept_bases += b

    print(
        json.dumps(
            {
                "records_seen": seen,
                "records_kept": kept,
                "bases_kept": kept_bases,
                "min_len": args.min_len,
                "max_len": args.max_len,
            },
            sort_keys=True,
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
