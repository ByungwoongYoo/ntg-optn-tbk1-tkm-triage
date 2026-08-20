#!/usr/bin/env python3
"""Construct a truth-blind baseline-preserving extension union.

The backbone is copied without replacing or trimming any record. Auxiliary candidates
are aligned to the backbone before any gold-standard access. Depending on the frozen
mode, the program appends either sufficiently non-contained whole contigs or only
terminal overhangs that extend beyond high-identity backbone alignments. Exact sequence
duplicates are rejected. No taxonomy or gold input is accepted.
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
    p.add_argument("--backbone", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--candidate-to-backbone-paf", required=True)
    p.add_argument("--mode", choices=("whole", "terminal", "hybrid"), required=True)
    p.add_argument("--min-candidate-length", type=int, default=500)
    p.add_argument("--min-terminal-bp", type=int, default=500)
    p.add_argument("--alignment-identity", type=float, default=0.97)
    p.add_argument("--alignment-min-bp", type=int, default=500)
    p.add_argument("--whole-max-aligned-fraction", type=float, default=0.50)
    p.add_argument("--containment-fraction", type=float, default=0.95)
    p.add_argument("--min-source-count", type=int, default=1)
    p.add_argument("--single-target-terminal", action="store_true")
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


def write_record(fh, name: str, seq: str, width: int = 80) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for left, right in sorted(intervals)[1:]:
        if left <= end:
            end = max(end, right)
        else:
            total += end - start
            start, end = left, right
    return total + end - start


def merged_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    result: list[tuple[int, int]] = []
    start, end = sorted(intervals)[0]
    for left, right in sorted(intervals)[1:]:
        if left <= end:
            end = max(end, right)
        else:
            result.append((start, end))
            start, end = left, right
    result.append((start, end))
    return result


def sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    backbone = list(read_fasta(Path(a.backbone)))
    candidates = dict(read_fasta(Path(a.candidates)))
    if not backbone or not candidates:
        raise SystemExit("empty backbone or candidate FASTA")
    if len({name for name, _ in backbone}) != len(backbone):
        raise SystemExit("duplicate backbone identifiers")
    for name, seq in backbone + list(candidates.items()):
        if not seq or not DNA.fullmatch(seq):
            raise SystemExit(f"invalid DNA sequence: {name}")

    metadata: dict[str, dict[str, str]] = {}
    with Path(a.metadata).open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = row.get("candidate_id") or row.get("representative_id")
            if not key:
                raise SystemExit("metadata lacks candidate identifier")
            metadata[key] = row
    missing = sorted(set(candidates) - set(metadata))
    if missing:
        raise SystemExit(f"candidate metadata missing: {missing[:10]}")

    alignments: dict[str, list[dict[str, object]]] = defaultdict(list)
    with Path(a.candidate_to_backbone_paf).open("rt", encoding="utf-8") as fh:
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            q, ql, qs, qe, strand, target, tl, ts, te, matches, aln_len, mapq = cols[:12]
            if q not in candidates:
                continue
            aln = int(aln_len)
            identity = int(matches) / aln if aln else 0.0
            if identity < a.alignment_identity or aln < a.alignment_min_bp:
                continue
            alignments[q].append(
                {
                    "ql": int(ql),
                    "qs": int(qs),
                    "qe": int(qe),
                    "strand": strand,
                    "target": target,
                    "tl": int(tl),
                    "ts": int(ts),
                    "te": int(te),
                    "identity": identity,
                    "alignment_length": aln,
                    "mapq": int(mapq),
                }
            )

    output_segments: list[tuple[str, str, str, str]] = []
    decisions: list[dict[str, object]] = []
    seen_hashes = {sha(seq) for _, seq in backbone}

    for cid, seq in sorted(candidates.items(), key=lambda item: (-len(item[1]), item[0])):
        meta = metadata[cid]
        source_count = int(meta.get("source_count", meta.get("assembler_count", "0")) or 0)
        alns = alignments.get(cid, [])
        intervals = merged_intervals([(int(row["qs"]), int(row["qe"])) for row in alns])
        aligned_bp = union_length(intervals)
        aligned_fraction = aligned_bp / len(seq) if seq else 0.0
        targets = sorted({str(row["target"]) for row in alns})
        exact_duplicate = sha(seq) in seen_hashes
        basic = (
            len(seq) >= a.min_candidate_length
            and seq.count("N") / len(seq) <= 0.05
            and source_count >= a.min_source_count
            and not exact_duplicate
        )
        segments: list[tuple[str, str]] = []
        reason = "selected"

        if not basic:
            if exact_duplicate:
                reason = "exact_duplicate"
            elif source_count < a.min_source_count:
                reason = "source_count_below_threshold"
            else:
                reason = "basic_qc_failed"
        elif a.mode == "whole":
            if aligned_fraction >= a.containment_fraction:
                reason = "contained_in_backbone"
            else:
                segments = [("whole", seq)]
        elif a.mode == "terminal":
            if not alns:
                segments = [("unaligned_whole", seq)]
            elif a.single_target_terminal and len(targets) != 1:
                reason = "multi_target_alignment"
            else:
                left = intervals[0][0]
                right = len(seq) - intervals[-1][1]
                if left >= a.min_terminal_bp:
                    segments.append(("left_terminal", seq[:left]))
                if right >= a.min_terminal_bp:
                    segments.append(("right_terminal", seq[len(seq) - right :]))
                if not segments:
                    reason = "no_terminal_extension"
        else:  # hybrid
            if not alns or aligned_fraction < a.whole_max_aligned_fraction:
                segments = [("whole_low_overlap", seq)]
            elif aligned_fraction >= a.containment_fraction:
                reason = "contained_in_backbone"
            elif a.single_target_terminal and len(targets) != 1:
                reason = "multi_target_alignment"
            else:
                left = intervals[0][0]
                right = len(seq) - intervals[-1][1]
                if left >= a.min_terminal_bp:
                    segments.append(("left_terminal", seq[:left]))
                if right >= a.min_terminal_bp:
                    segments.append(("right_terminal", seq[len(seq) - right :]))
                if not segments:
                    reason = "no_eligible_extension"

        accepted = 0
        accepted_bp = 0
        for segment_type, segment in segments:
            digest = sha(segment)
            if digest in seen_hashes:
                continue
            seen_hashes.add(digest)
            output_id = "LANTERN_V6_EXT_" + digest[:16].upper()
            output_segments.append((output_id, segment, cid, segment_type))
            accepted += 1
            accepted_bp += len(segment)
        if segments and accepted == 0:
            reason = "segment_duplicate"

        decisions.append(
            {
                "candidate_id": cid,
                "candidate_length": len(seq),
                "source_count": source_count,
                "sources": meta.get("sources", ""),
                "high_identity_alignment_count": len(alns),
                "aligned_query_bp": aligned_bp,
                "aligned_query_fraction": aligned_fraction,
                "backbone_target_count": len(targets),
                "backbone_targets": ",".join(targets[:20]),
                "selected_segment_count": accepted,
                "selected_segment_bp": accepted_bp,
                "selected": str(accepted > 0).lower(),
                "reason": reason if accepted == 0 else "selected",
                "mode": a.mode,
            }
        )

    output_path = out / "LANTERN_V6_EXTENSION_UNION.fasta"
    with output_path.open("wt", encoding="ascii") as fh:
        for name, seq in backbone:
            write_record(fh, name, seq)
        for output_id, seq, _source, _kind in output_segments:
            write_record(fh, output_id, seq)

    with (out / "EXTENSION_DECISIONS.tsv").open("w", newline="", encoding="utf-8") as fh:
        fields = list(decisions[0]) if decisions else ["candidate_id"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(decisions)
    with (out / "EXTENSION_SEGMENTS.tsv").open("w", newline="", encoding="utf-8") as fh:
        fields = ["output_id", "source_candidate_id", "segment_type", "length", "sha256"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for output_id, seq, source, kind in output_segments:
            writer.writerow(
                {
                    "output_id": output_id,
                    "source_candidate_id": source,
                    "segment_type": kind,
                    "length": len(seq),
                    "sha256": sha(seq),
                }
            )

    summary = {
        "status": "PASS",
        "mode": a.mode,
        "backbone_records": len(backbone),
        "backbone_bp": sum(len(seq) for _, seq in backbone),
        "candidate_records": len(candidates),
        "selected_segments": len(output_segments),
        "selected_segment_bp": sum(len(seq) for _, seq, _, _ in output_segments),
        "final_records": len(backbone) + len(output_segments),
        "final_bp": sum(len(seq) for _, seq in backbone)
        + sum(len(seq) for _, seq, _, _ in output_segments),
        "backbone_preserved_record_for_record": True,
        "parameters": {
            "min_candidate_length": a.min_candidate_length,
            "min_terminal_bp": a.min_terminal_bp,
            "alignment_identity": a.alignment_identity,
            "alignment_min_bp": a.alignment_min_bp,
            "whole_max_aligned_fraction": a.whole_max_aligned_fraction,
            "containment_fraction": a.containment_fraction,
            "min_source_count": a.min_source_count,
            "single_target_terminal": a.single_target_terminal,
        },
        "truth_input_accepted": False,
        "boundary": "Backbone-preserving construction from assembly-to-assembly alignments only; no gold, taxonomy, or performance input was accepted.",
    }
    (out / "EXTENSION_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out / "ASSEMBLY_SHA256.txt").write_text(
        f"{hashlib.sha256(output_path.read_bytes()).hexdigest()}  {output_path.name}\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
