#!/usr/bin/env python3
"""Construct a backbone-preserving LANTERN rescue assembly.

The method keeps every contig from a prespecified participant-level co-assembly
backbone and adds only nonredundant candidate contigs that meet an explicit evidence
configuration. It never replaces or shortens backbone contigs. This is intended to
prevent the v1 ensemble from losing correctly assembled backbone sequence while still
allowing longitudinal, multi-source, or long-read-supported rescue.

This script is truth blind: it accepts no gold-standard file and performs no biological
performance evaluation. A configuration may be developed on a public training pair,
but must be frozen before application to held-out pairs or the CAMI hidden challenge.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable

from fasta_utils import read_fasta, write_record


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--backbone", required=True)
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument("--selection", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--candidate-to-backbone-paf", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def as_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, "") or default)
    except ValueError:
        return default


def as_int(row: dict[str, str], key: str, default: int = 0) -> int:
    try:
        return int(float(row.get(key, "") or default))
    except ValueError:
        return default


def sha(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def load_tsv(path: str | Path, key: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            item = row[key]
            if item in rows:
                raise ValueError(f"Duplicate {key}: {item}")
            rows[item] = row
    return rows


def redundant_candidates(
    paf: str | Path,
    query_lengths: dict[str, int],
    minimum_identity: float,
    minimum_query_coverage: float,
) -> tuple[set[str], dict[str, dict[str, float | str]]]:
    redundant: set[str] = set()
    best: dict[str, dict[str, float | str]] = {}
    with open(paf) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q, qlen, qs, qe, strand, target, tlen, ts, te, matches, aln_len, mapq = fields[:12]
            if q not in query_lengths:
                continue
            qlen_i = int(qlen)
            aln_i = int(aln_len)
            matches_i = int(matches)
            identity = matches_i / aln_i if aln_i else 0.0
            query_coverage = aln_i / qlen_i if qlen_i else 0.0
            score = identity * query_coverage
            if q not in best or score > float(best[q]["score"]):
                best[q] = {
                    "target": target,
                    "identity": identity,
                    "query_coverage": query_coverage,
                    "alignment_length": aln_i,
                    "score": score,
                }
            if identity >= minimum_identity and query_coverage >= minimum_query_coverage:
                redundant.add(q)
    return redundant, best


def contains_token(csv_text: str, token: str) -> bool:
    return token in {item.strip() for item in str(csv_text).split(",") if item.strip()}


def candidate_gate(
    cid: str,
    seq: str,
    selection: dict[str, str],
    metadata: dict[str, str],
    config: dict,
    redundant: set[str],
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if cid in redundant:
        reasons.append("redundant_to_backbone")
    if config.get("require_v1_selected", True) and not truthy(selection.get("selected")):
        reasons.append("not_v1_selected")
    if len(seq) < int(config["minimum_contig_length"]):
        reasons.append("too_short")
    if seq.count("N") / max(1, len(seq)) > float(config.get("maximum_n_fraction", 0.01)):
        reasons.append("too_many_N")
    if as_float(selection, "score") < float(config["minimum_v1_score"]):
        reasons.append("score_below_threshold")
    if as_int(selection, "timepoints_supported") < int(config["minimum_timepoints_supported"]):
        reasons.append("insufficient_temporal_support")
    if as_float(selection, "max_short_breadth") < float(config["minimum_short_breadth"]):
        reasons.append("short_breadth_below_threshold")
    min_long_breadth = float(config.get("minimum_long_breadth", 0.0))
    min_long_spans = int(config.get("minimum_long_spanning_reads", 0))
    long_gate_enabled = min_long_breadth > 0 or min_long_spans > 0
    if long_gate_enabled:
        long_ok = (
            as_float(selection, "max_long_breadth") >= min_long_breadth
            or as_int(selection, "long_spanning_reads") >= min_long_spans
        )
        if not long_ok:
            reasons.append("long_support_below_threshold")
    if int(config.get("minimum_source_count", 1)) > as_int(selection, "source_count"):
        reasons.append("source_count_below_threshold")
    if int(config.get("minimum_assembler_count", 1)) > as_int(selection, "assembler_count"):
        reasons.append("assembler_count_below_threshold")
    if config.get("require_single_scope", False) and not contains_token(metadata.get("scopes", ""), "single"):
        reasons.append("no_single_scope_source")
    for token in config.get("exclude_source_ids", []):
        if contains_token(metadata.get("source_ids", ""), token):
            reasons.append(f"excluded_source:{token}")
    for token in config.get("require_any_source_ids", []):
        # Applied as an OR set below.
        pass
    required = list(config.get("require_any_source_ids", []))
    if required and not any(contains_token(metadata.get("source_ids", ""), token) for token in required):
        reasons.append("missing_required_source")
    return not reasons, reasons


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(args.config).read_text())
    backbone_records = list(read_fasta(args.backbone))
    candidate_records = list(read_fasta(args.candidate_fasta))
    backbone = dict(backbone_records)
    candidates = dict(candidate_records)
    if len(backbone) != len(backbone_records):
        raise ValueError("Duplicate backbone identifiers")
    if len(candidates) != len(candidate_records):
        raise ValueError("Duplicate candidate identifiers")
    selection = load_tsv(args.selection, "contig_id")
    metadata = load_tsv(args.metadata, "representative_id")
    if set(candidates) - set(selection):
        raise ValueError("Candidate FASTA contains IDs absent from selection table")
    if set(candidates) - set(metadata):
        raise ValueError("Candidate FASTA contains IDs absent from metadata table")

    redundant, best_alignments = redundant_candidates(
        args.candidate_to_backbone_paf,
        {cid: len(seq) for cid, seq in candidate_records},
        float(config["redundancy_minimum_identity"]),
        float(config["redundancy_minimum_query_coverage"]),
    )

    existing_hashes = {sha(seq) for seq in backbone.values()}
    audit_rows: list[dict] = []
    rescues: list[tuple[str, str]] = []
    for cid, seq in candidate_records:
        keep, reasons = candidate_gate(cid, seq, selection[cid], metadata[cid], config, redundant)
        seq_hash = sha(seq)
        if seq_hash in existing_hashes:
            keep = False
            reasons.append("exact_sequence_duplicate")
        alignment = best_alignments.get(cid, {})
        row = {
            "candidate_id": cid,
            "selected_as_rescue": str(keep).lower(),
            "length": len(seq),
            "sequence_sha256": seq_hash,
            "reasons_excluded": ",".join(sorted(set(reasons))),
            "v1_selected": selection[cid].get("selected", ""),
            "v1_rescue": selection[cid].get("rescue", ""),
            "v1_score": selection[cid].get("score", ""),
            "timepoints_supported": selection[cid].get("timepoints_supported", ""),
            "max_short_breadth": selection[cid].get("max_short_breadth", ""),
            "max_long_breadth": selection[cid].get("max_long_breadth", ""),
            "long_spanning_reads": selection[cid].get("long_spanning_reads", ""),
            "source_count": selection[cid].get("source_count", ""),
            "assembler_count": selection[cid].get("assembler_count", ""),
            "source_ids": metadata[cid].get("source_ids", ""),
            "assemblers": metadata[cid].get("assemblers", ""),
            "scopes": metadata[cid].get("scopes", ""),
            "best_backbone_target": alignment.get("target", ""),
            "best_backbone_identity": alignment.get("identity", ""),
            "best_backbone_query_coverage": alignment.get("query_coverage", ""),
        }
        audit_rows.append(row)
        if keep:
            rescues.append((cid, seq))
            existing_hashes.add(seq_hash)

    max_rescues = int(config.get("maximum_rescue_contigs", 0))
    max_rescue_bp = int(config.get("maximum_rescue_bp", 0))
    if max_rescues or max_rescue_bp:
        rescues.sort(
            key=lambda item: (
                -as_float(selection[item[0]], "score"),
                -as_int(selection[item[0]], "timepoints_supported"),
                -as_float(selection[item[0]], "max_short_breadth"),
                -len(item[1]),
                item[0],
            )
        )
        accepted: list[tuple[str, str]] = []
        total_bp = 0
        accepted_ids: set[str] = set()
        for item in rescues:
            if max_rescues and len(accepted) >= max_rescues:
                break
            if max_rescue_bp and total_bp + len(item[1]) > max_rescue_bp:
                continue
            accepted.append(item)
            accepted_ids.add(item[0])
            total_bp += len(item[1])
        for row in audit_rows:
            if row["selected_as_rescue"] == "true" and row["candidate_id"] not in accepted_ids:
                row["selected_as_rescue"] = "false"
                row["reasons_excluded"] = ",".join(filter(None, [row["reasons_excluded"], "rescue_cap"]))
        rescues = accepted

    assembly_path = out / "LANTERN_BACKBONE_RESCUE.fasta"
    with assembly_path.open("w") as fh:
        for name, seq in backbone_records:
            write_record(fh, f"BACKBONE_{name}", seq)
        for index, (name, seq) in enumerate(rescues, 1):
            write_record(fh, f"RESCUE_{index:06d}_{name}", seq)

    audit_rows.sort(
        key=lambda row: (
            row["selected_as_rescue"] != "true",
            -float(row["v1_score"] or 0),
            -int(row["length"]),
            row["candidate_id"],
        )
    )
    with (out / "RESCUE_CANDIDATE_AUDIT.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(audit_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    summary = {
        "method": config.get("method", "LANTERN-backbone-rescue"),
        "config_version": config.get("version"),
        "backbone_contigs": len(backbone_records),
        "backbone_bp": sum(len(seq) for _, seq in backbone_records),
        "candidate_contigs": len(candidate_records),
        "redundant_candidates": len(redundant),
        "rescue_contigs": len(rescues),
        "rescue_bp": sum(len(seq) for _, seq in rescues),
        "assembly_contigs": len(backbone_records) + len(rescues),
        "assembly_bp": sum(len(seq) for _, seq in backbone_records) + sum(len(seq) for _, seq in rescues),
        "assembly_sha256": hashlib.sha256(assembly_path.read_bytes()).hexdigest(),
        "config": config,
        "truth_blind": True,
        "claim_boundary": "Backbone-preserving sequence construction only; no gold-standard or hidden challenge result is used by this script.",
    }
    (out / "LANTERN_BACKBONE_RESCUE_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
