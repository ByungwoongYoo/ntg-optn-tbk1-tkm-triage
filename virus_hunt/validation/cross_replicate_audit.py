#!/usr/bin/env python3
"""Recompute cross-run sequence replication for strict Panax RdRP candidates.

This script intentionally separates candidate detection from discovery claims. It
verifies that each retained lineage has complete A/B/C palm motifs and asks whether
independently assembled runs contain closely matching protein and nucleotide regions.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq


def fasta(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def local_stats(a: str, b: str, alphabet: str) -> dict[str, float | int | str]:
    aligner = Align.PairwiseAligner(mode="local")
    if alphabet == "nt":
        aligner.match_score = 2.0
        aligner.mismatch_score = -2.0
        aligner.open_gap_score = -5.0
        aligner.extend_gap_score = -1.0
        variants = [(b, "+"), (str(Seq(b).reverse_complement()), "-")]
    else:
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -5.0
        aligner.extend_gap_score = -0.5
        variants = [(b, "+")]

    best: dict[str, float | int | str] | None = None
    for target, orientation in variants:
        alignment = aligner.align(a, target)[0]
        matches = 0
        aligned = 0
        for (s1, e1), (s2, e2) in zip(alignment.aligned[0], alignment.aligned[1]):
            x = a[s1:e1]
            y = target[s2:e2]
            matches += sum(c1 == c2 for c1, c2 in zip(x, y))
            aligned += len(x)
        result: dict[str, float | int | str] = {
            "identity": matches / aligned if aligned else 0.0,
            "aligned": aligned,
            "coverage_1": aligned / len(a) if a else 0.0,
            "coverage_2": aligned / len(b) if b else 0.0,
            "orientation_2": orientation,
            "score": float(alignment.score),
        }
        if best is None or float(result["score"]) > float(best["score"]):
            best = result
    assert best is not None
    return best


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--members", type=Path, required=True)
    p.add_argument("--aa", type=Path, required=True)
    p.add_argument("--nt", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    aa = fasta(args.aa)
    nt = fasta(args.nt)
    rows: list[dict[str, str]] = list(csv.DictReader(args.members.open(), delimiter="\t"))
    by_lineage: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_lineage[row["lineage"]].append(row)

    pair_fields = [
        "lineage", "query_1", "query_2", "accession_1", "accession_2",
        "aa_identity", "aa_aligned", "aa_coverage_1", "aa_coverage_2",
        "nt_identity", "nt_aligned", "nt_coverage_1", "nt_coverage_2", "nt_orientation_2",
    ]
    pair_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for lineage, members in sorted(by_lineage.items()):
        pairs = []
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                r1, r2 = members[i], members[j]
                if r1["accession"] == r2["accession"]:
                    continue
                a = local_stats(aa[r1["query"]], aa[r2["query"]], "aa")
                n = local_stats(nt[r1["nt_contig"]], nt[r2["nt_contig"]], "nt")
                out = {
                    "lineage": lineage,
                    "query_1": r1["query"], "query_2": r2["query"],
                    "accession_1": r1["accession"], "accession_2": r2["accession"],
                    "aa_identity": f"{100*float(a['identity']):.3f}",
                    "aa_aligned": str(a["aligned"]),
                    "aa_coverage_1": f"{100*float(a['coverage_1']):.3f}",
                    "aa_coverage_2": f"{100*float(a['coverage_2']):.3f}",
                    "nt_identity": f"{100*float(n['identity']):.3f}",
                    "nt_aligned": str(n["aligned"]),
                    "nt_coverage_1": f"{100*float(n['coverage_1']):.3f}",
                    "nt_coverage_2": f"{100*float(n['coverage_2']):.3f}",
                    "nt_orientation_2": str(n["orientation_2"]),
                }
                pair_rows.append(out)
                pairs.append(out)

        runs = sorted({m["accession"] for m in members})
        complete_abc = all(
            m.get("pssm_ABC") in {"ABC", "CAB"}
            and m.get("motif_A") and m.get("motif_B") and m.get("motif_C")
            for m in members
        )
        strong_pairs = [
            x for x in pairs
            if float(x["aa_identity"]) >= 90.0
            and min(float(x["aa_coverage_1"]), float(x["aa_coverage_2"])) >= 50.0
            and float(x["nt_identity"]) >= 90.0
            and min(float(x["nt_coverage_1"]), float(x["nt_coverage_2"])) >= 50.0
        ]
        replication = "strong" if len(runs) >= 2 and complete_abc and strong_pairs else "not_established"
        summary_rows.append({
            "lineage": lineage,
            "member_count": str(len(members)),
            "distinct_runs": str(len(runs)),
            "runs": ";".join(runs),
            "complete_motif_evidence_all_members": str(complete_abc).lower(),
            "strong_cross_run_pairs": str(len(strong_pairs)),
            "cross_run_replication": replication,
            "minimum_palmdb_identity": f"{min(float(m['pident']) for m in members):.3f}",
            "maximum_palmdb_identity": f"{max(float(m['pident']) for m in members):.3f}",
        })

    with (args.out_dir / "cross_run_pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields, delimiter="\t")
        writer.writeheader(); writer.writerows(pair_rows)
    summary_fields = list(summary_rows[0]) if summary_rows else ["lineage"]
    with (args.out_dir / "lineage_replication_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader(); writer.writerows(summary_rows)
    (args.out_dir / "lineage_replication_summary.json").write_text(json.dumps(summary_rows, indent=2) + "\n")
    print(json.dumps(summary_rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
