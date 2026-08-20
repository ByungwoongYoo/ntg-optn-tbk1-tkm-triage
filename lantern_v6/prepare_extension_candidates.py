#!/usr/bin/env python3
"""Prepare exact-deduplicated auxiliary contigs for LANTERN-v6 extension rescue.

This script is truth blind. It accepts one or more SOURCE=FASTA assemblies, filters by
length/N fraction, collapses exact sequence duplicates across sources, and emits stable
sequence-hash identifiers plus source-count metadata. It does not use taxonomy, gold
assemblies, or performance results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

DNA = re.compile(r"^[ACGTN]+$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--source", action="append", required=True, help="SOURCE=FASTA")
    p.add_argument("--min-length", type=int, default=500)
    p.add_argument("--max-n-fraction", type=float, default=0.05)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name = None
    parts: list[str] = []
    with path.open("rt", encoding="ascii", errors="strict") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier at line {line_no}")
                parts = []
            else:
                if name is None:
                    raise ValueError(f"sequence before header at line {line_no}")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def write_record(fh, name: str, sequence: str, width: int = 80) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(sequence), width):
        fh.write(sequence[i : i + width] + "\n")


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    by_hash: dict[str, dict[str, object]] = {}
    source_stats: list[dict[str, object]] = []

    for spec in a.source:
        source, raw_path = spec.split("=", 1)
        path = Path(raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        seen = kept = bases = 0
        for original_id, seq in read_fasta(path):
            seen += 1
            if not seq or not DNA.fullmatch(seq):
                raise ValueError(f"invalid DNA sequence: {source}:{original_id}")
            if len(seq) < a.min_length or seq.count("N") / len(seq) > a.max_n_fraction:
                continue
            kept += 1
            bases += len(seq)
            digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
            row = by_hash.setdefault(
                digest,
                {
                    "sequence": seq,
                    "members": [],
                    "sources": set(),
                    "original_ids": [],
                },
            )
            row["members"].append((source, original_id))
            row["sources"].add(source)
            row["original_ids"].append(original_id)
        source_stats.append(
            {
                "source": source,
                "path": str(path),
                "records_seen": seen,
                "records_kept": kept,
                "kept_bp": bases,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )

    ordered = sorted(by_hash.items(), key=lambda item: (-len(item[1]["sequence"]), item[0]))
    metadata: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    with (out / "extension_candidates.fasta").open("wt", encoding="ascii") as fasta:
        for digest, row in ordered:
            cid = "LANTERN_V6_C" + digest[:16].upper()
            seq = row["sequence"]
            write_record(fasta, cid, seq)
            sources = sorted(row["sources"])
            metadata.append(
                {
                    "candidate_id": cid,
                    "representative_id": cid,
                    "length": len(seq),
                    "n_fraction": seq.count("N") / len(seq),
                    "source_count": len(sources),
                    "assembler_count": len(sources),
                    "sources": ",".join(sources),
                    "sequence_sha256": digest,
                    "exact_member_count": len(row["members"]),
                }
            )
            for source, original_id in sorted(row["members"]):
                member_rows.append(
                    {
                        "candidate_id": cid,
                        "source": source,
                        "original_id": original_id,
                        "sequence_sha256": digest,
                    }
                )

    for name, rows, fields in [
        (
            "extension_candidate_metadata.tsv",
            metadata,
            [
                "candidate_id",
                "representative_id",
                "length",
                "n_fraction",
                "source_count",
                "assembler_count",
                "sources",
                "sequence_sha256",
                "exact_member_count",
            ],
        ),
        (
            "extension_candidate_members.tsv",
            member_rows,
            ["candidate_id", "source", "original_id", "sequence_sha256"],
        ),
        (
            "source_stats.tsv",
            source_stats,
            ["source", "path", "records_seen", "records_kept", "kept_bp", "sha256"],
        ),
    ]:
        with (out / name).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "status": "PASS",
        "minimum_length": a.min_length,
        "maximum_n_fraction": a.max_n_fraction,
        "sources": [row["source"] for row in source_stats],
        "input_records": sum(int(row["records_seen"]) for row in source_stats),
        "filtered_records": sum(int(row["records_kept"]) for row in source_stats),
        "exact_deduplicated_candidates": len(metadata),
        "candidate_bp": sum(int(row["length"]) for row in metadata),
        "truth_input_accepted": False,
        "boundary": "Exact sequence deduplication and basic sequence QC only; no gold, taxonomy, or performance information was used.",
    }
    (out / "PREPARATION_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
