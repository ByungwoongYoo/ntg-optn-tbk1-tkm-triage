#!/usr/bin/env python3
"""Summarize post-freeze abundance-stratified genome recovery.

Abundance files may contain different sample IDs.  The previous implementation chose
output columns from the first row only, so a later row containing (for example)
``sample15_read_percent`` failed CSV serialization after the first row exposed only
``sample14_read_percent``.  This version takes the ordered union of fields across all
rows and performs consistency checks before writing.

Abundance and genome labels are evaluation-only metadata.  They must be opened only
after truth-blind assemblies and hashes have been frozen.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--abundance", action="append", required=True)
    p.add_argument(
        "--recovery", action="append", required=True,
        help="METHOD=per_genome_recovery.tsv",
    )
    p.add_argument("--low-percent", type=float, default=0.01)
    p.add_argument("--high-percent", type=float, default=0.1)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ordered_field_union(rows: list[dict[str, Any]], fallback: list[str]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields or fallback


def write_tsv(path: Path, rows: list[dict[str, Any]], fallback: list[str]) -> None:
    fields = ordered_field_union(rows, fallback)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    a = parse_args()
    if not (0 <= a.low_percent < a.high_percent):
        raise SystemExit("require 0 <= low-percent < high-percent")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    by_genome: dict[str, dict[str, float]] = defaultdict(dict)
    abundance_inputs: list[dict[str, Any]] = []
    for raw_path in a.abundance:
        path = Path(raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        rows = read_tsv(path)
        sample_ids = sorted({row["sample_id"] for row in rows})
        abundance_inputs.append(
            {"path": str(path), "rows": len(rows), "sample_ids": sample_ids}
        )
        for row in rows:
            genome = row["genome_id"]
            sample = row["sample_id"]
            value = float(row["relative_read_percent"])
            if sample in by_genome[genome] and by_genome[genome][sample] != value:
                raise ValueError(f"conflicting abundance value for {genome}/{sample}")
            by_genome[genome][sample] = value

    profile: dict[str, dict[str, Any]] = {}
    all_sample_ids = sorted({sample for values in by_genome.values() for sample in values})
    for genome, sample_values in by_genome.items():
        values = list(sample_values.values())
        mean = sum(values) / len(values)
        minimum = min(values)
        maximum = max(values)
        if mean < a.low_percent:
            abundance_stratum = "low"
        elif mean < a.high_percent:
            abundance_stratum = "medium"
        else:
            abundance_stratum = "high"
        if maximum < a.low_percent:
            trajectory = "persistent_low"
        elif minimum < a.low_percent <= maximum:
            trajectory = "longitudinal_rescue_opportunity"
        else:
            trajectory = "consistently_detectable"
        row: dict[str, Any] = {
            "genome_id": genome,
            "n_samples": len(values),
            "mean_read_percent": mean,
            "min_read_percent": minimum,
            "max_read_percent": maximum,
            "abundance_stratum": abundance_stratum,
            "trajectory_stratum": trajectory,
        }
        for sample in all_sample_ids:
            row[f"{sample}_read_percent"] = sample_values.get(sample, "")
        profile[genome] = row

    details: list[dict[str, Any]] = []
    methods: list[str] = []
    for spec in a.recovery:
        method, raw_path = spec.split("=", 1)
        if method in methods:
            raise ValueError(f"duplicate recovery method: {method}")
        methods.append(method)
        path = Path(raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        for row in read_tsv(path):
            genome = row["genome_id"]
            if genome not in profile:
                continue
            details.append(
                {
                    "method": method,
                    **profile[genome],
                    "recovery_fraction": float(row["recovery_fraction"]),
                    "recovered_50": row["recovered_50"],
                    "recovered_90": row["recovered_90"],
                }
            )

    groups = [
        "low", "medium", "high", "persistent_low",
        "longitudinal_rescue_opportunity", "consistently_detectable",
    ]
    summaries: list[dict[str, Any]] = []
    for method in methods:
        for group in groups:
            if group in {"low", "medium", "high"}:
                rows = [
                    row for row in details
                    if row["method"] == method and row["abundance_stratum"] == group
                ]
            else:
                rows = [
                    row for row in details
                    if row["method"] == method and row["trajectory_stratum"] == group
                ]
            values = [float(row["recovery_fraction"]) for row in rows]
            summaries.append(
                {
                    "method": method,
                    "stratum": group,
                    "n_genomes": len(rows),
                    "mean_recovery": sum(values) / len(values) if values else None,
                    "median_recovery": statistics.median(values) if values else None,
                    "minimum_recovery": min(values) if values else None,
                    "recovered_50": sum(row["recovered_50"] == "true" for row in rows),
                    "recovered_90": sum(row["recovered_90"] == "true" for row in rows),
                }
            )

    write_tsv(
        out / "ABUNDANCE_RECOVERY_DETAIL.tsv", details,
        ["method", "genome_id", "recovery_fraction"],
    )
    write_tsv(
        out / "ABUNDANCE_STRATIFIED_SUMMARY.tsv", summaries,
        ["method", "stratum", "n_genomes", "mean_recovery"],
    )
    report = {
        "status": "PASS",
        "n_abundance_genomes": len(profile),
        "n_joined_rows": len(details),
        "sample_ids": all_sample_ids,
        "methods": methods,
        "abundance_inputs": abundance_inputs,
        "low_threshold_percent": a.low_percent,
        "high_threshold_percent": a.high_percent,
        "field_union_fix": True,
        "boundary": (
            "Abundance and genome labels are gold-standard evaluation metadata accessed "
            "after truth-blind output freeze. Low-abundance strata do not influence selection."
        ),
    }
    (out / "ABUNDANCE_AUDIT.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
