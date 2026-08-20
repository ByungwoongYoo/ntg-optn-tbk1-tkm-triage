#!/usr/bin/env python3
"""Build a baseline-preserving additive longitudinal-rescue assembly.

The baseline assembly is copied without substitution. Candidate contigs are added only
from clusters that do not contain the selected baseline source. This prevents the
representative-selection step from replacing a baseline contig with another assembler's
representative. No gold-standard, taxonomic reference, or performance result is read.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from fasta_utils import read_fasta, write_record


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--baseline-fasta", required=True)
    p.add_argument("--representative-fasta", required=True)
    p.add_argument("--representative-metadata", required=True)
    p.add_argument("--selection", required=True)
    p.add_argument("--baseline-source-id", default="megahit_pair")
    p.add_argument(
        "--mode",
        choices=("rescue_all", "selected_nonbaseline", "rescue_high", "rescue_strict", "rescue_length"),
        required=True,
    )
    p.add_argument("--min-score", type=float, default=0.0)
    p.add_argument("--min-breadth", type=float, default=0.0)
    p.add_argument("--min-length", type=int, default=1000)
    p.add_argument("--min-timepoints", type=int, default=1)
    p.add_argument("--out", required=True)
    return p.parse_args()


def sha(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def load_tsv(path: str, key: str) -> dict[str, dict[str, str]]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    out = {row[key]: row for row in rows}
    if len(out) != len(rows):
        raise ValueError(f"duplicate {key} in {path}")
    return out


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    baseline = list(read_fasta(a.baseline_fasta))
    representatives = dict(read_fasta(a.representative_fasta))
    metadata = load_tsv(a.representative_metadata, "representative_id")
    selection = load_tsv(a.selection, "contig_id")
    if set(representatives) != set(metadata) or set(representatives) != set(selection):
        raise ValueError("representative FASTA, metadata, and selection IDs differ")
    if not baseline:
        raise ValueError("empty baseline FASTA")

    seen_hashes = {sha(seq) for _, seq in baseline}
    additions: list[tuple[str, str, dict[str, str]]] = []
    audit_rows: list[dict[str, object]] = []

    for cid, seq in representatives.items():
        m = metadata[cid]
        s = selection[cid]
        source_ids = set(filter(None, m.get("source_ids", "").split(",")))
        selected = s.get("selected", "").lower() == "true"
        rescue = s.get("rescue", "").lower() == "true"
        score = float(s.get("score") or 0)
        breadth = float(s.get("max_short_breadth") or 0)
        timepoints = int(float(s.get("timepoints_supported") or 0))
        length = len(seq)

        eligible_kind = selected if a.mode == "selected_nonbaseline" else rescue
        passes = (
            eligible_kind
            and a.baseline_source_id not in source_ids
            and score >= a.min_score
            and breadth >= a.min_breadth
            and length >= a.min_length
            and timepoints >= a.min_timepoints
        )
        reason = "eligible"
        if not eligible_kind:
            reason = "not_selected_kind"
        elif a.baseline_source_id in source_ids:
            reason = "baseline_cluster_present"
        elif score < a.min_score:
            reason = "score_below_threshold"
        elif breadth < a.min_breadth:
            reason = "breadth_below_threshold"
        elif length < a.min_length:
            reason = "length_below_threshold"
        elif timepoints < a.min_timepoints:
            reason = "timepoints_below_threshold"
        elif sha(seq) in seen_hashes:
            reason = "exact_sequence_duplicate"
            passes = False

        audit_rows.append(
            {
                "contig_id": cid,
                "accepted": str(passes).lower(),
                "reason": reason,
                "selected": str(selected).lower(),
                "rescue": str(rescue).lower(),
                "score": score,
                "max_short_breadth": breadth,
                "timepoints_supported": timepoints,
                "length": length,
                "source_ids": ",".join(sorted(source_ids)),
                "sequence_sha256": sha(seq),
            }
        )
        if passes:
            seen_hashes.add(sha(seq))
            additions.append((cid, seq, s))

    assembly = out / "ADDITIVE_ASSEMBLY.fasta"
    manifest_rows: list[dict[str, object]] = []
    with assembly.open("w") as f:
        for i, (original_id, seq) in enumerate(baseline, 1):
            sid = f"BASE_{i:08d}"
            write_record(f, sid, seq)
            manifest_rows.append(
                {
                    "output_id": sid,
                    "origin": "baseline",
                    "original_id": original_id,
                    "length": len(seq),
                    "sequence_sha256": sha(seq),
                }
            )
        for i, (cid, seq, _) in enumerate(additions, 1):
            sid = f"RESCUE_{i:08d}"
            write_record(f, sid, seq)
            manifest_rows.append(
                {
                    "output_id": sid,
                    "origin": "additive_rescue",
                    "original_id": cid,
                    "length": len(seq),
                    "sequence_sha256": sha(seq),
                }
            )

    for path, rows in [
        (out / "ADDITIVE_MANIFEST.tsv", manifest_rows),
        (out / "ADDITIVE_CANDIDATE_AUDIT.tsv", audit_rows),
    ]:
        fields = list(rows[0]) if rows else ["output_id"]
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    summary = {
        "mode": a.mode,
        "baseline_source_id": a.baseline_source_id,
        "thresholds": {
            "min_score": a.min_score,
            "min_breadth": a.min_breadth,
            "min_length": a.min_length,
            "min_timepoints": a.min_timepoints,
        },
        "baseline_contigs": len(baseline),
        "baseline_bp": sum(len(seq) for _, seq in baseline),
        "added_contigs": len(additions),
        "added_bp": sum(len(seq) for _, seq, _ in additions),
        "output_contigs": len(manifest_rows),
        "output_bp": sum(int(row["length"]) for row in manifest_rows),
        "truth_blind": True,
        "boundary": "Baseline-preserving additive construction uses only assembly provenance and frozen mapping evidence; no gold or taxonomy is read.",
    }
    (out / "ADDITIVE_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
