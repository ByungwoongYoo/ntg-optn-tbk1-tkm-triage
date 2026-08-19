#!/usr/bin/env python3
"""Stream discovery-only unitig output and freeze a candidate set before validation.

The parser deliberately extracts accession identifiers and nucleotide sequences using robust
regular expressions so it can tolerate minor unitig-caller/pyseer formatting differences.
Validation samples are never read by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ACC_RE = re.compile(r"GC[AF]_\d+(?:\.\d+)?")
SEQ_RE = re.compile(r"(?<![A-Za-z])([ACGTNacgtn]{15,})(?![A-Za-z])")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sparse", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--p-screen", type=float, default=1e-4)
    p.add_argument("--max-candidates", type=int, default=5000)
    p.add_argument("--min-present", type=int, default=5)
    return p.parse_args()


def normalize_acc(v: str) -> str:
    m = ACC_RE.search(str(v))
    return m.group(0) if m else str(v)


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(a.manifest, dtype=str).drop_duplicates("assembly_ID")
    manifest = manifest[manifest["split"].eq("discovery")].copy()
    manifest["phenotype"] = manifest["phenotype"].str.upper()
    manifest = manifest[manifest["phenotype"].isin(["R", "S"])]
    manifest["assembly_ID"] = manifest["assembly_ID"].map(normalize_acc)
    label = dict(zip(manifest["assembly_ID"], manifest["phenotype"]))
    r_ids = {k for k, v in label.items() if v == "R"}
    s_ids = {k for k, v in label.items() if v == "S"}
    nr, ns = len(r_ids), len(s_ids)
    if nr < 20 or ns < 20:
        raise SystemExit(f"Insufficient discovery classes R={nr} S={ns}")

    # Only a finite set of count combinations exists. Precompute exact one-sided P values.
    p_lookup: dict[tuple[int, int], float] = {}
    for rp in range(nr + 1):
        for sp in range(ns + 1):
            if rp + sp < a.min_present:
                continue
            p_lookup[(rp, sp)] = float(
                fisher_exact([[rp, sp], [nr - rp, ns - sp]], alternative="greater").pvalue
            )

    heap: list[tuple[float, int, dict[str, object]]] = []
    n_tests = 0
    n_parse_fail = 0
    n_no_accessions = 0
    seq_seen: set[str] = set()

    with open(a.sparse, "rt", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            seq_match = SEQ_RE.search(line)
            if not seq_match:
                n_parse_fail += 1
                continue
            seq = seq_match.group(1).upper()
            if seq in seq_seen:
                continue
            seq_seen.add(seq)
            ids = {normalize_acc(x) for x in ACC_RE.findall(line)}
            if not ids:
                n_no_accessions += 1
                continue
            rp = len(ids & r_ids)
            sp = len(ids & s_ids)
            if rp + sp < a.min_present:
                continue
            n_tests += 1
            pval = p_lookup[(rp, sp)]
            # Haldane-Anscombe correction for ranking only.
            orank = ((rp + 0.5) * (ns - sp + 0.5)) / ((sp + 0.5) * (nr - rp + 0.5))
            if orank <= 1 or pval > a.p_screen:
                continue
            rec = {
                "sequence": seq,
                "sequence_length": len(seq),
                "R_present": rp,
                "S_present": sp,
                "R_absent": nr - rp,
                "S_absent": ns - sp,
                "unadjusted_or": orank,
                "one_sided_fisher_p": pval,
                "source_line": line_no,
            }
            # max-heap emulation with negative P; keep the smallest P values.
            item = (-pval, -len(seq), rec)
            if len(heap) < a.max_candidates:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)

    records = [x[2] for x in heap]
    records.sort(key=lambda z: (float(z["one_sided_fisher_p"]), -float(z["unadjusted_or"]), -int(z["sequence_length"])))
    for i, rec in enumerate(records, start=1):
        rec["candidate_id"] = f"UGWAS_{i:05d}"
        rec["bonferroni_p"] = min(1.0, float(rec["one_sided_fisher_p"]) * max(1, n_tests))
        rec["discovery_bonferroni"] = rec["bonferroni_p"] <= 0.05

    table = pd.DataFrame(records)
    table.to_csv(out / "SELECTED_DISCOVERY_UNITIGS.csv", index=False)
    with open(out / "selected_unitigs.fasta", "w") as fh:
        for rec in records:
            fh.write(f">{rec['candidate_id']}\n{rec['sequence']}\n")

    summary = {
        "n_discovery_R": nr,
        "n_discovery_S": ns,
        "n_unique_parseable_unitigs_tested": n_tests,
        "n_selected": len(records),
        "n_discovery_bonferroni": int(sum(bool(r["discovery_bonferroni"]) for r in records)),
        "n_parse_fail": n_parse_fail,
        "n_no_accessions": n_no_accessions,
        "p_screen": a.p_screen,
        "max_candidates": a.max_candidates,
        "validation_touched": False,
        "boundary": "Discovery-only unadjusted screening. Population correction and untouched validation are mandatory before any association claim.",
    }
    (out / "SCREEN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "SPARSE_FORMAT_EXAMPLES.txt").write_text("Parser used accession and nucleotide regexes; raw sparse input is not redistributed.\n")

    hashes = []
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
