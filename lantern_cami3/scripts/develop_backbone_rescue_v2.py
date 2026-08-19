#!/usr/bin/env python3
"""Develop a backbone-preserving LANTERN v2 config on a public Toy training pair.

This script is deliberately a *development* tool: it may inspect the public Toy gold
standard to select one configuration. The resulting config is not evidence of external
performance. It must be frozen and then applied without retuning to different Toy
individuals and, subject to data-access authorization, the CAMI hidden challenge.

For efficiency, backbone and candidate assemblies are aligned to the gold standard once.
Each grid configuration is then scored exactly from the subset of candidate alignments
that it would add to the immutable backbone.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import random
from collections import defaultdict
from pathlib import Path

from fasta_utils import read_fasta
from lantern_backbone_rescue import (
    as_float,
    as_int,
    candidate_gate,
    load_tsv,
    redundant_candidates,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--backbone-fasta", required=True)
    p.add_argument("--candidate-fasta", required=True)
    p.add_argument("--selection", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--candidate-to-backbone-paf", required=True)
    p.add_argument("--backbone-to-truth-paf", required=True)
    p.add_argument("--candidate-to-truth-paf", required=True)
    p.add_argument("--truth-mapping", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=20260820)
    p.add_argument("--bootstrap", type=int, default=2000)
    return p.parse_args()


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    start, end = intervals[0]
    total = 0
    for a, b in intervals[1:]:
        if a <= end:
            end = max(end, b)
        else:
            total += end - start
            start, end = a, b
    return total + end - start


def load_truth(path: str | Path) -> tuple[dict[str, str], dict[str, int]]:
    sequence_to_genome: dict[str, str] = {}
    genome_lengths: dict[str, int] = defaultdict(int)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sequence_to_genome[row["sequence_id"]] = row["genome_id"]
            genome_lengths[row["genome_id"]] += int(row["length"])
    if not sequence_to_genome:
        raise ValueError("Empty truth mapping")
    return sequence_to_genome, dict(genome_lengths)


def accepted_paf_hits(
    path: str | Path,
    allowed_queries: set[str],
    truth: dict[str, str],
    minimum_identity: float = 0.90,
    minimum_alignment: int = 500,
) -> dict[str, list[tuple[str, int, int, int, int]]]:
    """Return q -> (target sequence, target start, target end, query start, query end)."""
    hits: dict[str, list[tuple[str, int, int, int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q, qlen, qs, qe, strand, target, tlen, ts, te, matches, aln_len, mapq = fields[:12]
            if q not in allowed_queries or target not in truth:
                continue
            aln = int(aln_len)
            identity = int(matches) / aln if aln else 0.0
            if identity < minimum_identity or aln < minimum_alignment:
                continue
            hits[q].append((target, int(ts), int(te), int(qs), int(qe)))
    return hits


def assembly_metrics(
    query_lengths: dict[str, int],
    hits: dict[str, list[tuple[str, int, int, int, int]]],
    truth: dict[str, str],
    genome_lengths: dict[str, int],
) -> tuple[dict, dict[str, float]]:
    per_truth_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    query_genome_intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    accepted_alignment_bp = 0
    for query in query_lengths:
        for target, ts, te, qs, qe in hits.get(query, []):
            per_truth_intervals[target].append((ts, te))
            query_genome_intervals[query][truth[target]].append((qs, qe))
            accepted_alignment_bp += te - ts

    genome_covered: dict[str, int] = defaultdict(int)
    for target, intervals in per_truth_intervals.items():
        genome_covered[truth[target]] += union_length(intervals)
    recovery = {
        genome: genome_covered.get(genome, 0) / length if length else 0.0
        for genome, length in genome_lengths.items()
    }
    assembly_bp = sum(query_lengths.values())
    chimeric_bp = 0
    chimeric_contigs = 0
    aligned_contigs = 0
    for query, qlen in query_lengths.items():
        by_genome = query_genome_intervals.get(query, {})
        support = sorted((union_length(v), genome) for genome, v in by_genome.items())
        if support:
            aligned_contigs += 1
        support.sort(reverse=True)
        second = support[1][0] if len(support) > 1 else 0
        threshold = max(1000, int(0.10 * qlen))
        if len(support) > 1 and second >= threshold:
            chimeric_contigs += 1
            chimeric_bp += qlen
    truth_bp = sum(genome_lengths.values())
    unique_bp = sum(genome_covered.values())
    values = list(recovery.values())
    summary = {
        "n_genomes": len(genome_lengths),
        "truth_bp_total": truth_bp,
        "unique_truth_bp_covered": unique_bp,
        "genome_fraction_percent": 100 * unique_bp / truth_bp if truth_bp else 0.0,
        "mean_genome_recovery": sum(values) / len(values) if values else 0.0,
        "median_genome_recovery": sorted(values)[len(values) // 2] if values else 0.0,
        "genomes_recovered_50": sum(value >= 0.50 for value in values),
        "genomes_recovered_90": sum(value >= 0.90 for value in values),
        "assembly_contigs_total": len(query_lengths),
        "assembly_total_bp": assembly_bp,
        "assembly_contigs_with_accepted_alignment": aligned_contigs,
        "assembly_contigs_without_accepted_alignment": len(query_lengths) - aligned_contigs,
        "accepted_alignment_bp": accepted_alignment_bp,
        "alignment_to_unique_truth_ratio": accepted_alignment_bp / unique_bp if unique_bp else None,
        "cross_binid_chimeric_contigs": chimeric_contigs,
        "cross_binid_chimeric_bp": chimeric_bp,
        "cross_binid_chimeric_bp_fraction": chimeric_bp / assembly_bp if assembly_bp else 0.0,
    }
    return summary, recovery


def paired_bootstrap(
    full: dict[str, float],
    baseline: dict[str, float],
    seed: int,
    repetitions: int,
) -> tuple[float, float, float]:
    genomes = sorted(set(full) & set(baseline))
    differences = [100 * (full[g] - baseline[g]) for g in genomes]
    observed = sum(differences) / len(differences) if differences else 0.0
    rng = random.Random(seed)
    draws: list[float] = []
    if differences:
        for _ in range(repetitions):
            sample = [differences[rng.randrange(len(differences))] for _ in differences]
            draws.append(sum(sample) / len(sample))
    draws.sort()
    if not draws:
        return observed, 0.0, 0.0
    lo = draws[max(0, int(0.025 * len(draws)) - 1)]
    hi = draws[min(len(draws) - 1, int(0.975 * len(draws)))]
    return observed, lo, hi


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    backbone_records = list(read_fasta(args.backbone_fasta))
    candidate_records = list(read_fasta(args.candidate_fasta))
    backbone = dict(backbone_records)
    candidates = dict(candidate_records)
    selection = load_tsv(args.selection, "contig_id")
    metadata = load_tsv(args.metadata, "representative_id")
    if set(candidates) - set(selection) or set(candidates) - set(metadata):
        raise ValueError("Candidate IDs are not aligned across FASTA, selection, and metadata")

    truth, genome_lengths = load_truth(args.truth_mapping)
    backbone_hits = accepted_paf_hits(args.backbone_to_truth_paf, set(backbone), truth)
    candidate_hits = accepted_paf_hits(args.candidate_to_truth_paf, set(candidates), truth)
    redundant, best_backbone = redundant_candidates(
        args.candidate_to_backbone_paf,
        {cid: len(seq) for cid, seq in candidate_records},
        0.97,
        0.85,
    )
    baseline_summary, baseline_recovery = assembly_metrics(
        {name: len(seq) for name, seq in backbone_records},
        backbone_hits,
        truth,
        genome_lengths,
    )
    (out / "BACKBONE_BASELINE_SUMMARY.json").write_text(json.dumps(baseline_summary, indent=2) + "\n")

    score_grid = [0.28, 0.34, 0.40, 0.46, 0.52, 0.60]
    breadth_grid = [0.50, 0.60, 0.70, 0.80, 0.90]
    length_grid = [1000, 1500, 2000, 3000]
    source_count_grid = [1, 2]
    require_single_grid = [False, True]
    rows: list[dict] = []
    recoveries: dict[str, dict[str, float]] = {}
    configs: dict[str, dict] = {}

    for index, (score_min, breadth_min, length_min, source_min, require_single) in enumerate(
        itertools.product(score_grid, breadth_grid, length_grid, source_count_grid, require_single_grid), 1
    ):
        config_id = f"C{index:04d}"
        config = {
            "version": f"LANTERN-backbone-rescue-development-{config_id}",
            "method": "LANTERN-backbone-rescue-v2",
            "minimum_contig_length": length_min,
            "maximum_n_fraction": 0.05,
            "minimum_v1_score": score_min,
            "minimum_timepoints_supported": 2,
            "minimum_short_breadth": breadth_min,
            "minimum_long_breadth": 0.0,
            "minimum_long_spanning_reads": 0,
            "minimum_source_count": source_min,
            "minimum_assembler_count": 1,
            "require_single_scope": require_single,
            "exclude_source_ids": ["megahit_pair"],
            "require_any_source_ids": [],
            "require_v1_selected": True,
            "redundancy_minimum_identity": 0.97,
            "redundancy_minimum_query_coverage": 0.85,
            "maximum_rescue_contigs": 0,
            "maximum_rescue_bp": 0,
            "development_pair": "CAMI III Toy samples 0 and 1",
            "development_uses_public_gold": True,
        }
        selected_ids: list[str] = []
        for cid, seq in candidate_records:
            keep, reasons = candidate_gate(cid, seq, selection[cid], metadata[cid], config, redundant)
            if keep and hashlib.sha256(seq.encode("ascii")).hexdigest() not in {
                hashlib.sha256(value.encode("ascii")).hexdigest() for value in backbone.values()
            }:
                selected_ids.append(cid)
        query_lengths = {name: len(seq) for name, seq in backbone_records}
        query_lengths.update({cid: len(candidates[cid]) for cid in selected_ids})
        hits = defaultdict(list)
        for q, values in backbone_hits.items():
            hits[q].extend(values)
        for q in selected_ids:
            hits[q].extend(candidate_hits.get(q, []))
        summary, recovery = assembly_metrics(query_lengths, hits, truth, genome_lengths)
        relative_chimera = (
            summary["cross_binid_chimeric_bp_fraction"] / baseline_summary["cross_binid_chimeric_bp_fraction"] - 1
            if baseline_summary["cross_binid_chimeric_bp_fraction"]
            else 0.0
        )
        gf_gain = summary["genome_fraction_percent"] - baseline_summary["genome_fraction_percent"]
        mean_gain = 100 * (summary["mean_genome_recovery"] - baseline_summary["mean_genome_recovery"])
        row = {
            "config_id": config_id,
            "minimum_v1_score": score_min,
            "minimum_short_breadth": breadth_min,
            "minimum_contig_length": length_min,
            "minimum_source_count": source_min,
            "require_single_scope": require_single,
            "n_rescue": len(selected_ids),
            "rescue_bp": sum(len(candidates[cid]) for cid in selected_ids),
            "genome_fraction_percent": summary["genome_fraction_percent"],
            "genome_fraction_gain_percentage_points": gf_gain,
            "mean_genome_recovery_percent": 100 * summary["mean_genome_recovery"],
            "mean_genome_recovery_gain_percentage_points": mean_gain,
            "recovered_50": summary["genomes_recovered_50"],
            "recovered_90": summary["genomes_recovered_90"],
            "chimera_bp_fraction": summary["cross_binid_chimeric_bp_fraction"],
            "relative_chimera_change": relative_chimera,
            "unaligned_contigs": summary["assembly_contigs_without_accepted_alignment"],
            "alignment_to_unique_truth_ratio": summary["alignment_to_unique_truth_ratio"],
            "assembly_total_bp": summary["assembly_total_bp"],
            "assembly_bp_ratio_to_backbone": summary["assembly_total_bp"] / baseline_summary["assembly_total_bp"],
            "training_feasible": (
                relative_chimera <= 0.10
                and summary["assembly_total_bp"] <= 1.50 * baseline_summary["assembly_total_bp"]
                and summary["genomes_recovered_90"] >= baseline_summary["genomes_recovered_90"]
            ),
        }
        # Primary ranking is gold recovery; penalties only break ties and prevent an
        # unconstrained union from winning by adding large amounts of unsupported DNA.
        row["development_objective"] = (
            gf_gain
            + 0.25 * mean_gain
            + 0.05 * (summary["genomes_recovered_50"] - baseline_summary["genomes_recovered_50"])
            - 0.25 * max(0.0, relative_chimera)
            - 0.02 * max(0.0, row["assembly_bp_ratio_to_backbone"] - 1.0)
        )
        rows.append(row)
        configs[config_id] = config
        recoveries[config_id] = recovery

    feasible = [row for row in rows if row["training_feasible"]]
    if not feasible:
        feasible = rows
    feasible.sort(
        key=lambda row: (
            -row["genome_fraction_gain_percentage_points"],
            -row["mean_genome_recovery_gain_percentage_points"],
            row["relative_chimera_change"],
            row["n_rescue"],
            row["config_id"],
        )
    )
    best = feasible[0]
    best_config = configs[best["config_id"]]
    best_config["version"] = "LANTERN-backbone-rescue-v2-training-freeze-20260820"
    best_config["selected_training_config_id"] = best["config_id"]
    best_config["training_selection_rule"] = (
        "Highest public-training genome-fraction gain among configurations satisfying the prespecified "
        "chimera, assembly-size, and recovered-90 feasibility constraints; ties use mean recovery, "
        "lower chimera, and fewer rescue contigs."
    )
    best_config["external_validation_required"] = True
    best_config["claim_boundary"] = (
        "Selected using public CAMI III Toy samples 0 and 1. This is a training result, not external validation. "
        "The config must be committed unchanged before held-out-pair evaluation."
    )
    observed, lo, hi = paired_bootstrap(
        recoveries[best["config_id"]], baseline_recovery, args.seed, args.bootstrap
    )
    best["paired_bootstrap_mean_gain_percentage_points"] = observed
    best["paired_bootstrap_95_ci_low"] = lo
    best["paired_bootstrap_95_ci_high"] = hi

    rows.sort(key=lambda row: (-row["development_objective"], row["config_id"]))
    with (out / "BACKBONE_RESCUE_GRID_RESULTS.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "LANTERN_BACKBONE_RESCUE_V2_CONFIG.json").write_text(json.dumps(best_config, indent=2) + "\n")
    decision = {
        "status": "TRAINING_CONFIG_SELECTED",
        "best_training_result": best,
        "backbone_baseline": baseline_summary,
        "n_grid_configurations": len(rows),
        "n_feasible_configurations": sum(row["training_feasible"] for row in rows),
        "external_validation_required": True,
        "claim_boundary": best_config["claim_boundary"],
    }
    (out / "V2_DEVELOPMENT_DECISION.json").write_text(json.dumps(decision, indent=2) + "\n")
    with (out / "TOP_30_GRID_RESULTS.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows[:30])
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
