#!/usr/bin/env python3
"""Summarize compact per-reference raw-read support from a mapped BAM."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
from pathlib import Path

from Bio import SeqIO


def count(cmd: list[str]) -> int:
    return int(subprocess.check_output(cmd, text=True).strip() or "0")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--bam", type=Path, required=True)
    p.add_argument("--refs", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    refs = {r.id: len(r.seq) for r in SeqIO.parse(args.refs, "fasta")}
    depths: dict[str, list[int]] = {r: [0] * n for r, n in refs.items()}
    proc = subprocess.Popen(
        ["samtools", "depth", "-aa", str(args.bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        r, pos, dep = line.rstrip().split("\t")[:3]
        if r in depths:
            i = int(pos) - 1
            if 0 <= i < len(depths[r]):
                depths[r][i] = int(dep)
    if proc.wait() != 0:
        raise RuntimeError("samtools depth failed")

    metadata_text = args.metadata.read_text(encoding="utf-8", errors="replace") if args.metadata.exists() else ""
    read_count = ""
    try:
        meta_rows = list(csv.DictReader(metadata_text.splitlines(), delimiter="\t"))
        if meta_rows:
            read_count = meta_rows[0].get("read_count", "")
    except Exception:
        pass

    fields = [
        "run", "reference", "reference_length", "input_read_count",
        "mapped_primary", "mapped_mapq20", "proper_pair_primary",
        "breadth_1x", "breadth_5x", "breadth_10x", "mean_depth", "median_depth",
        "depth_p10", "depth_p90", "depth_cv", "zero_bases",
    ]
    rows = []
    for ref, vals in depths.items():
        n = len(vals)
        mean = statistics.fmean(vals) if vals else 0.0
        med = statistics.median(vals) if vals else 0.0
        ordered = sorted(vals)
        def quantile(pct: float) -> float:
            if not ordered:
                return 0.0
            idx = min(len(ordered) - 1, max(0, int(round((len(ordered)-1)*pct))))
            return float(ordered[idx])
        sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        mapped = count(["samtools", "view", "-c", "-F", "2308", str(args.bam), ref])
        mapq20 = count(["samtools", "view", "-c", "-F", "2308", "-q", "20", str(args.bam), ref])
        proper = count(["samtools", "view", "-c", "-f", "2", "-F", "2308", str(args.bam), ref])
        rows.append({
            "run": args.run,
            "reference": ref,
            "reference_length": n,
            "input_read_count": read_count,
            "mapped_primary": mapped,
            "mapped_mapq20": mapq20,
            "proper_pair_primary": proper,
            "breadth_1x": f"{sum(x >= 1 for x in vals)/n:.6f}" if n else "0",
            "breadth_5x": f"{sum(x >= 5 for x in vals)/n:.6f}" if n else "0",
            "breadth_10x": f"{sum(x >= 10 for x in vals)/n:.6f}" if n else "0",
            "mean_depth": f"{mean:.6f}",
            "median_depth": f"{med:.6f}",
            "depth_p10": f"{quantile(0.10):.6f}",
            "depth_p90": f"{quantile(0.90):.6f}",
            "depth_cv": f"{(sd/mean if mean else 0.0):.6f}",
            "zero_bases": sum(x == 0 for x in vals),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    (args.output.with_suffix(".json")).write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
