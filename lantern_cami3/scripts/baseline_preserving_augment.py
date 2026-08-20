#!/usr/bin/env python3
"""Truth-blind baseline-preserving LANTERN augmentation.

The backbone assembly is copied byte-for-byte at the sequence-record level. Candidate
contigs from auxiliary assemblers are appended only when they are non-redundant to the
backbone and pass pre-specified evidence gates. The program never replaces a backbone
contig with a longer representative from another assembler.

This is designed to correct the failure mode observed in LANTERN-v1 where clustering
could replace a strong backbone contig and slightly reduce gold-standard genome fraction.
No truth/gold-standard input is accepted by this program.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterator

DNA_RE = re.compile(r"^[ACGTN]+$")


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    parts: list[str] = []
    with open(path, "rt", encoding="ascii") as fh:
        for line_number, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ValueError(f"empty FASTA identifier at line {line_number}")
                parts = []
            else:
                if name is None:
                    raise ValueError(f"sequence before header at line {line_number}")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def write_record(fh, name: str, sequence: str, width: int = 80) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(sequence), width):
        fh.write(sequence[i : i + width] + "\n")


def sha256_text(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def safe_id(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.|:-]+", "_", text)
    return value or "unnamed"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--backbone", required=True, help="Strong baseline assembly FASTA")
    p.add_argument("--candidates", required=True, help="Auxiliary candidate FASTA")
    p.add_argument("--metadata", required=True, help="Candidate metadata TSV")
    p.add_argument("--evidence", required=True, help="Per-candidate/per-sample evidence TSV")
    p.add_argument("--candidate-to-backbone-paf", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)
    p.add_argument(
        "--ablation",
        choices=("full", "no_longitudinal", "no_long", "no_consensus"),
        default="full",
    )
    return p.parse_args()


def load_metadata(path: str | Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = row.get("representative_id") or row.get("contig_id") or row.get("candidate_id")
            if not key:
                raise ValueError("metadata lacks representative_id/contig_id/candidate_id")
            if key in result:
                raise ValueError(f"duplicate metadata id: {key}")
            result[key] = row
    return result


def load_evidence(path: str | Path) -> tuple[dict[str, dict[str, dict[str, float]]], list[str]]:
    evidence: dict[str, dict[str, dict[str, float]]] = defaultdict(dict)
    samples: set[str] = set()
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "contig_id" not in reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ValueError("evidence TSV requires contig_id and sample_id")
        for row in reader:
            cid = row["contig_id"]
            sample = row["sample_id"]
            samples.add(sample)
            values: dict[str, float] = {}
            for key, value in row.items():
                if key in {"contig_id", "sample_id"}:
                    continue
                try:
                    values[key] = float(value or 0)
                except ValueError:
                    values[key] = 0.0
            evidence[cid][sample] = values
    if not samples:
        raise ValueError("no evidence samples")
    return evidence, sorted(samples)


def parse_redundancy(
    paf_path: str | Path,
    candidate_lengths: dict[str, int],
    backbone_ids: set[str],
    min_identity: float,
    min_shorter_coverage: float,
) -> tuple[set[str], dict[str, dict[str, float | str]]]:
    redundant: set[str] = set()
    best: dict[str, dict[str, float | str]] = {}
    with open(paf_path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            q, ql, _qs, _qe, _strand, t, tl, _ts, _te, matches, aln_len, _mapq = cols[:12]
            if q not in candidate_lengths or t not in backbone_ids:
                continue
            ql_i = int(ql)
            tl_i = int(tl)
            aln_i = int(aln_len)
            matches_i = int(matches)
            identity = matches_i / aln_i if aln_i else 0.0
            coverage = aln_i / min(ql_i, tl_i) if min(ql_i, tl_i) else 0.0
            score = identity * coverage
            previous = best.get(q)
            if previous is None or score > float(previous["score"]):
                best[q] = {
                    "target": t,
                    "identity": identity,
                    "shorter_coverage": coverage,
                    "score": score,
                    "alignment_length": aln_i,
                }
            if identity >= min_identity and coverage >= min_shorter_coverage:
                redundant.add(q)
    return redundant, best


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))

    backbone_records = list(read_fasta(args.backbone))
    candidate_records = list(read_fasta(args.candidates))
    if not backbone_records:
        raise ValueError("empty backbone assembly")
    if not candidate_records:
        raise ValueError("empty candidate assembly")

    backbone_ids = {name for name, _ in backbone_records}
    if len(backbone_ids) != len(backbone_records):
        raise ValueError("duplicate backbone IDs")
    backbone_hashes = {sha256_text(seq) for _, seq in backbone_records}
    candidate_sequences = dict(candidate_records)
    if len(candidate_sequences) != len(candidate_records):
        raise ValueError("duplicate candidate IDs")
    for name, seq in backbone_records + candidate_records:
        if not seq or not DNA_RE.fullmatch(seq):
            raise ValueError(f"invalid DNA sequence for {name}")

    metadata = load_metadata(args.metadata)
    missing_meta = sorted(set(candidate_sequences) - set(metadata))
    if missing_meta:
        raise ValueError(f"candidate metadata missing for {missing_meta[:10]}")
    evidence, samples = load_evidence(args.evidence)

    redundant, best_alignment = parse_redundancy(
        args.candidate_to_backbone_paf,
        {key: len(seq) for key, seq in candidate_sequences.items()},
        backbone_ids,
        float(config.get("cluster_identity", 0.97)),
        float(config.get("cluster_shorter_coverage", 0.85)),
    )

    weights = dict(config["weights"])
    if args.ablation == "no_longitudinal":
        weights["temporal_recurrence"] = 0.0
    if args.ablation == "no_long":
        weights["long_breadth"] = 0.0
        weights["long_spanning"] = 0.0
    if args.ablation == "no_consensus":
        weights["source_consensus"] = 0.0
    weight_sum = sum(float(v) for v in weights.values())
    if weight_sum <= 0:
        raise ValueError("zero weight sum")
    weights = {key: float(value) / weight_sum for key, value in weights.items()}

    rows: list[dict[str, object]] = []
    selected_ids: list[str] = []
    selected_hashes: set[str] = set(backbone_hashes)

    for cid, seq in candidate_records:
        meta = metadata[cid]
        by_sample = evidence.get(cid, {})
        short_breadths: list[float] = []
        long_breadths: list[float] = []
        supports: list[bool] = []
        long_spanning = 0
        for sample in samples:
            ev = by_sample.get(sample, {})
            sb = float(ev.get("short_breadth", 0.0))
            sd = float(ev.get("short_depth", 0.0))
            lb = 0.0 if args.ablation == "no_long" else float(ev.get("long_breadth", 0.0))
            ld = 0.0 if args.ablation == "no_long" else float(ev.get("long_depth", 0.0))
            sp = 0 if args.ablation == "no_long" else int(ev.get("long_spanning_reads", 0.0))
            short_breadths.append(sb)
            long_breadths.append(lb)
            long_spanning += sp
            supports.append(
                max(sb, lb) >= float(config["support_breadth"])
                and max(sd, ld) >= float(config["support_depth"])
            )

        timepoints_supported = sum(supports)
        assembler_count = int(meta.get("assembler_count", 0) or 0)
        scopes = set(filter(None, (meta.get("scopes", "") or "").split(",")))
        source_consensus = min(1.0, assembler_count / 3.0)
        features = {
            "length": min(1.0, math.log1p(len(seq)) / math.log1p(100000)),
            "source_consensus": source_consensus,
            "temporal_recurrence": timepoints_supported / max(1, len(samples)),
            "short_breadth": max(short_breadths or [0.0]),
            "long_breadth": max(long_breadths or [0.0]),
            "long_spanning": min(1.0, math.log1p(long_spanning) / math.log1p(6)),
        }
        score = sum(weights[key] * features[key] for key in weights)

        consensus_gate = (
            args.ablation != "no_consensus"
            and assembler_count >= int(config["minimum_assembler_sources_for_consensus"])
        )
        temporal_gate = (
            args.ablation != "no_longitudinal"
            and timepoints_supported >= int(config["minimum_rescue_timepoints"])
        )
        long_gate = (
            args.ablation != "no_long"
            and long_spanning >= int(config["minimum_unique_long_spans"])
        )
        single_supported = "single" in scopes and timepoints_supported >= 1
        basic_gate = (
            len(seq) >= int(config["minimum_contig_length"])
            and seq.count("N") / len(seq) <= float(config["maximum_n_fraction"])
            and max(short_breadths + long_breadths + [0.0]) >= float(config["support_breadth"])
        )
        sequence_hash = sha256_text(seq)
        exact_duplicate = sequence_hash in selected_hashes
        redundancy_gate = cid not in redundant and not exact_duplicate
        evidence_gate = consensus_gate or temporal_gate or long_gate or single_supported
        selected = (
            basic_gate
            and redundancy_gate
            and score >= float(config["selection_score_minimum"])
            and evidence_gate
        )

        reason = "selected"
        if not basic_gate:
            reason = "failed_basic_gate"
        elif exact_duplicate:
            reason = "exact_backbone_or_prior_duplicate"
        elif cid in redundant:
            reason = "alignment_redundant_to_backbone"
        elif score < float(config["selection_score_minimum"]):
            reason = "score_below_threshold"
        elif not evidence_gate:
            reason = "no_independent_evidence_gate"

        if selected:
            selected_ids.append(cid)
            selected_hashes.add(sequence_hash)

        best = best_alignment.get(cid, {})
        rows.append(
            {
                "contig_id": cid,
                "selected": str(selected).lower(),
                "reason": reason,
                "score": score,
                "length": len(seq),
                "n_fraction": seq.count("N") / len(seq),
                "timepoints_supported": timepoints_supported,
                "n_timepoints": len(samples),
                "assembler_count": assembler_count,
                "max_short_breadth": max(short_breadths or [0.0]),
                "max_long_breadth": max(long_breadths or [0.0]),
                "long_spanning_reads": long_spanning,
                "consensus_gate": str(consensus_gate).lower(),
                "temporal_gate": str(temporal_gate).lower(),
                "long_gate": str(long_gate).lower(),
                "exact_duplicate": str(exact_duplicate).lower(),
                "redundant_to_backbone": str(cid in redundant).lower(),
                "best_backbone_target": best.get("target", ""),
                "best_backbone_identity": best.get("identity", ""),
                "best_backbone_shorter_coverage": best.get("shorter_coverage", ""),
                "scopes": meta.get("scopes", ""),
                "assemblers": meta.get("assemblers", ""),
                "ablation": args.ablation,
                **{f"feature_{key}": value for key, value in features.items()},
            }
        )

    output_fasta = out / "LANTERN_BACKBONE_AUGMENTED.fasta"
    with output_fasta.open("wt", encoding="ascii") as fh:
        for name, seq in backbone_records:
            write_record(fh, name, seq)
        for cid in selected_ids:
            name = "LANTERN_AUG_" + sha256_text(candidate_sequences[cid])[:16].upper()
            write_record(fh, safe_id(name), candidate_sequences[cid])

    rows.sort(key=lambda row: (row["selected"] != "true", -float(row["score"]), -int(row["length"]), str(row["contig_id"])))
    with (out / "AUGMENTATION_DECISIONS.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "status": "PASS",
        "method": "truth-blind baseline-preserving augmentation",
        "ablation": args.ablation,
        "backbone_records": len(backbone_records),
        "backbone_bp": sum(len(seq) for _, seq in backbone_records),
        "candidate_records": len(candidate_records),
        "selected_augmentation_records": len(selected_ids),
        "selected_augmentation_bp": sum(len(candidate_sequences[cid]) for cid in selected_ids),
        "final_records": len(backbone_records) + len(selected_ids),
        "final_bp": sum(len(seq) for _, seq in backbone_records) + sum(len(candidate_sequences[cid]) for cid in selected_ids),
        "backbone_sequences_preserved_exactly": True,
        "alignment_redundant_candidates": len(redundant),
        "exact_duplicate_candidates_rejected": sum(row["reason"] == "exact_backbone_or_prior_duplicate" for row in rows),
        "config_version": config.get("version"),
        "truth_input_accepted": False,
        "claim_boundary": "Truth-blind construction only. Gold-standard performance must be evaluated separately after freeze.",
    }
    (out / "AUGMENTATION_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
