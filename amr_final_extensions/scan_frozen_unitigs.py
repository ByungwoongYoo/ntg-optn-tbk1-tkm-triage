#!/usr/bin/env python3
"""Scan frozen unitig candidates exactly in BV-BRC genome FASTA files.

Both forward and reverse-complement sequences are queried with a single Aho-Corasick
automaton. Candidate sequences and the external cohort were frozen before this scan.
Exact short-sequence presence is an association feature, not a mechanism or functional
annotation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import ahocorasick
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-panel", required=True)
    p.add_argument("--cohort", required=True)
    p.add_argument("--assemblies-dir", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
    return sequence.upper().translate(table)[::-1]


def fasta_records(path: Path):
    header = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(args.candidate_panel, dtype=str).fillna("")
    unitigs = panel[panel["candidate_type"].eq("unitig")].copy()
    if unitigs.empty:
        raise ValueError("No frozen unitigs in candidate panel")
    if unitigs["candidate_id"].duplicated().any():
        raise ValueError("Duplicate candidate_id")
    unitigs["sequence"] = unitigs["matrix_key"].str.upper()
    if unitigs["sequence"].duplicated().any():
        raise ValueError("Duplicate frozen unitig sequences")
    invalid = unitigs[~unitigs["sequence"].str.fullmatch(r"[ACGTRYMKBDHVN]+")]
    if len(invalid):
        raise ValueError(f"Invalid unitig sequences: {invalid.candidate_id.tolist()[:10]}")

    automaton = ahocorasick.Automaton()
    pattern_map: dict[str, list[tuple[str, str]]] = {}
    for _, row in unitigs.iterrows():
        candidate = str(row["candidate_id"])
        seq = str(row["sequence"])
        rc = reverse_complement(seq)
        pattern_map.setdefault(seq, []).append((candidate, "+"))
        pattern_map.setdefault(rc, []).append((candidate, "-"))
    for pattern, values in pattern_map.items():
        automaton.add_word(pattern, (pattern, values))
    automaton.make_automaton()

    cohort = pd.read_csv(args.cohort, dtype=str).fillna("")
    cohort = cohort.drop_duplicates("genome_id")
    genome_ids = cohort["genome_id"].astype(str).tolist()
    matrix = pd.DataFrame(
        0,
        index=unitigs["candidate_id"].astype(str),
        columns=genome_ids,
        dtype=np.uint8,
    )
    evidence: list[dict] = []
    missing: list[str] = []
    for i, genome_id in enumerate(genome_ids, start=1):
        fasta = Path(args.assemblies_dir) / f"{safe_name(genome_id)}.fna"
        if not fasta.exists():
            missing.append(genome_id)
            continue
        found: set[str] = set()
        for contig, sequence in fasta_records(fasta):
            for end, (_pattern, values) in automaton.iter(sequence):
                for candidate, orientation in values:
                    if candidate in found:
                        continue
                    start = end - len(_pattern) + 1
                    found.add(candidate)
                    matrix.loc[candidate, genome_id] = 1
                    evidence.append({
                        "candidate_id": candidate,
                        "genome_id": genome_id,
                        "contig": contig,
                        "start_0based": int(start),
                        "end_0based_inclusive": int(end),
                        "orientation": orientation,
                    })
                if len(found) == len(unitigs):
                    break
            if len(found) == len(unitigs):
                break
        if i % 25 == 0 or i == len(genome_ids):
            print(f"scanned={i}/{len(genome_ids)}", flush=True)

    matrix.to_csv(out / "BVBRC_FROZEN_UNITIGS.Rtab", sep="\t")
    pd.DataFrame(evidence).to_csv(out / "BVBRC_FROZEN_UNITIG_FIRST_HITS.csv", index=False)
    unitigs[["candidate_id", "sequence"]].to_csv(out / "BVBRC_FROZEN_UNITIG_MANIFEST.csv", index=False)
    (out / "MISSING_FASTA_GENOME_IDS.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    prevalence = pd.DataFrame({
        "candidate_id": matrix.index,
        "carriers": matrix.sum(axis=1).astype(int).to_numpy(),
        "frequency": matrix.mean(axis=1).to_numpy(),
    })
    prevalence.to_csv(out / "BVBRC_FROZEN_UNITIG_PREVALENCE.csv", index=False)
    summary = {
        "n_candidates": int(len(unitigs)),
        "n_genomes": int(len(genome_ids)),
        "n_missing_fasta": int(len(missing)),
        "n_candidates_observed": int((prevalence["carriers"] > 0).sum()),
        "n_exact_candidate_genome_hits": int(matrix.to_numpy().sum()),
        "boundary": "Exact unitig occurrence is not a functional annotation, causal determinant, or novelty claim.",
    }
    (out / "BVBRC_FROZEN_UNITIG_SCAN_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
