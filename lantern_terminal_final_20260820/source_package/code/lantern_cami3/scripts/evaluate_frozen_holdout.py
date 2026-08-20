#!/usr/bin/env python3
"""Evaluate one frozen LANTERN method on an untouched Toy holdout.

The script opens only completed gold-side evaluation tables. It cannot alter assemblies,
thresholds, or candidate selection. The strongest eligible baseline is selected by the
predeclared genome-fraction ordering, and all strict gates are read from the frozen
holdout configuration.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--method", action="append", required=True, help="METHOD=EVALUATION_DIR")
    p.add_argument("--manifest", required=True)
    p.add_argument("--holdout-config", required=True)
    p.add_argument("--full-method", required=True)
    p.add_argument("--no-longitudinal-method", required=True)
    p.add_argument("--no-long-method", default=None)
    p.add_argument("--no-consensus-method", default=None)
    p.add_argument("--abundance-summary", required=True)
    p.add_argument("--pseudo-summary", required=True)
    p.add_argument("--seed", type=int, default=20260819)
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_csv(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * p
    lo = math.floor(position)
    hi = math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def bootstrap_mean_ci(diffs: list[float], n: int, seed: int) -> tuple[float | None, float | None, float | None]:
    if not diffs:
        return None, None, None
    rng = random.Random(seed)
    m = len(diffs)
    samples = [sum(diffs[rng.randrange(m)] for _ in range(m)) / m for _ in range(n)]
    return sum(diffs) / m, percentile(samples, 0.025), percentile(samples, 0.975)


def relative_change(new: float, old: float) -> float | None:
    if old == 0:
        return 0.0 if new == 0 else None
    return (new - old) / old


def load_abundance(path: Path) -> dict[tuple[str, str], float]:
    result: dict[tuple[str, str], float] = {}
    for row in read_tsv(path):
        value = row.get("mean_recovery")
        if value not in (None, "", "None"):
            result[(row["method"], row["stratum"])] = 100.0 * float(value)
    return result


def load_pseudo(path: Path) -> dict[tuple[str, str], float]:
    result: dict[tuple[str, str], float] = {}
    for row in read_tsv(path):
        value = row.get("mean_recovery")
        if value not in (None, "", "None"):
            result[(row["method"], row["tier"])] = 100.0 * float(value)
    return result


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(a.holdout_config).read_text(encoding="utf-8"))
    strict = config["strict_success_gates"]
    manifest_rows = read_tsv(a.manifest)
    manifest = {row["method"]: row for row in manifest_rows}

    method_dirs: dict[str, Path] = {}
    for spec in a.method:
        name, path = spec.split("=", 1)
        if name in method_dirs:
            raise ValueError(f"duplicate method: {name}")
        if name not in manifest:
            raise ValueError(f"method absent from manifest: {name}")
        method_dirs[name] = Path(path)

    summaries: dict[str, dict[str, Any]] = {}
    recoveries: dict[str, dict[str, float]] = {}
    for name, directory in method_dirs.items():
        summary_path = directory / "GOLD_COVERAGE_SUMMARY.json"
        recovery_path = directory / "per_genome_recovery.tsv"
        if not summary_path.is_file() or not recovery_path.is_file():
            raise FileNotFoundError(f"missing evaluation output for {name}: {directory}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary.update(manifest[name])
        summary["method"] = name
        summaries[name] = summary
        recoveries[name] = {
            row["genome_id"]: float(row["recovery_fraction"])
            for row in read_tsv(recovery_path)
        }

    baselines = [
        summaries[name]
        for name in summaries
        if manifest[name].get("role") == "baseline"
        and manifest[name].get("eligible_primary_baseline", "true").lower() == "true"
    ]
    if not baselines:
        raise SystemExit("no eligible baseline")
    strongest = max(
        baselines,
        key=lambda row: (
            float(row["genome_fraction_percent"]),
            float(row["mean_genome_recovery"]),
            -float(row.get("cross_binid_chimeric_bp_fraction") or 0.0),
            row["method"],
        ),
    )
    full = summaries[a.full_method]
    no_longitudinal = summaries[a.no_longitudinal_method]

    common = sorted(set(recoveries[a.full_method]) & set(recoveries[strongest["method"]]))
    paired_differences_pp = [
        100.0 * (recoveries[a.full_method][gid] - recoveries[strongest["method"]][gid])
        for gid in common
    ]
    paired_mean, ci_low, ci_high = bootstrap_mean_ci(paired_differences_pp, a.bootstrap, a.seed)

    gf_gain = float(full["genome_fraction_percent"]) - float(strongest["genome_fraction_percent"])
    mean_gain = 100.0 * (
        float(full["mean_genome_recovery"]) - float(strongest["mean_genome_recovery"])
    )
    full_chimera = float(full.get("cross_binid_chimeric_bp_fraction") or 0.0)
    baseline_chimera = float(strongest.get("cross_binid_chimeric_bp_fraction") or 0.0)
    chimera_change = relative_change(full_chimera, baseline_chimera)
    longitudinal_drop = float(full["genome_fraction_percent"]) - float(
        no_longitudinal["genome_fraction_percent"]
    )

    abundance = load_abundance(Path(a.abundance_summary))
    low_full = abundance.get((a.full_method, "low"))
    low_baseline = abundance.get((strongest["method"], "low"))
    low_gain = None if low_full is None or low_baseline is None else low_full - low_baseline

    pseudo = load_pseudo(Path(a.pseudo_summary))
    tiers = [
        "exact_target_withheld",
        "species_neighbour_withheld",
        "genus_neighbour_withheld",
        "deep_neighbour_withheld",
    ]
    pseudo_gains: dict[str, float | None] = {}
    for tier in tiers:
        full_value = pseudo.get((a.full_method, tier))
        baseline_value = pseudo.get((strongest["method"], tier))
        pseudo_gains[tier] = (
            None if full_value is None or baseline_value is None else full_value - baseline_value
        )
    pseudo_all_positive = all(value is not None and value > 0 for value in pseudo_gains.values())

    gates = {
        "genome_fraction_gain": gf_gain
        >= float(strict["genome_fraction_gain_minimum_percentage_points"]),
        "mean_genome_recovery_gain": mean_gain
        >= float(strict["mean_genome_recovery_gain_minimum_percentage_points"]),
        "chimera_tradeoff": chimera_change is not None
        and chimera_change <= float(strict["relative_cross_binid_chimera_increase_maximum"]),
        "longitudinal_ablation_contribution": longitudinal_drop
        >= float(strict["longitudinal_ablation_drop_minimum_percentage_points"]),
        "low_abundance_gain": low_gain is not None
        and low_gain >= float(strict["low_abundance_gain_minimum_percentage_points"]),
        "pseudo_novel_all_tiers": pseudo_all_positive,
        "paired_bootstrap_positive": ci_low is not None and ci_low > 0,
    }
    strict_success = all(gates.values())
    status = "SUCCESS" if strict_success else (
        "PARTIAL" if gf_gain > 0 or mean_gain > 0 or pseudo_all_positive else "NEGATIVE"
    )

    decision = {
        "status": status,
        "strict_success": strict_success,
        "holdout_pair": config["holdout_pair"],
        "full_method": a.full_method,
        "strongest_baseline": strongest["method"],
        "genome_fraction_gain_percentage_points": gf_gain,
        "mean_genome_recovery_gain_percentage_points": mean_gain,
        "paired_per_genome_mean_gain_percentage_points": paired_mean,
        "paired_bootstrap_95_ci": [ci_low, ci_high],
        "n_paired_genomes": len(common),
        "full_chimera_bp_fraction": full_chimera,
        "baseline_chimera_bp_fraction": baseline_chimera,
        "relative_chimera_change": chimera_change,
        "longitudinal_ablation_drop_percentage_points": longitudinal_drop,
        "low_abundance_full_mean_recovery_percent": low_full,
        "low_abundance_baseline_mean_recovery_percent": low_baseline,
        "low_abundance_gain_percentage_points": low_gain,
        "pseudo_novel_gains_percentage_points": pseudo_gains,
        "gates": gates,
        "strict_success_gates": strict,
        "claim_boundary": config["claim_boundary"],
    }
    (out / "HOLDOUT_DECISION.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    comparison_rows: list[dict[str, Any]] = []
    for name, summary in sorted(
        summaries.items(),
        key=lambda item: (
            item[1].get("role") != "baseline",
            -float(item[1]["genome_fraction_percent"]),
            item[0],
        ),
    ):
        comparison_rows.append(
            {
                "method": name,
                "role": summary.get("role", ""),
                "mode": summary.get("mode", ""),
                "scope": summary.get("scope", ""),
                "ablation": summary.get("ablation", ""),
                "genome_fraction_percent": summary["genome_fraction_percent"],
                "mean_genome_recovery_percent": 100.0 * float(summary["mean_genome_recovery"]),
                "median_genome_recovery_percent": 100.0 * float(summary["median_genome_recovery"]),
                "genomes_recovered_50": summary["genomes_recovered_50"],
                "genomes_recovered_90": summary["genomes_recovered_90"],
                "assembly_contigs_total": summary["assembly_contigs_total"],
                "assembly_total_bp": summary["assembly_total_bp"],
                "cross_binid_chimeric_bp_fraction": summary["cross_binid_chimeric_bp_fraction"],
                "alignment_to_unique_truth_ratio": summary["alignment_to_unique_truth_ratio"],
                "strongest_baseline": str(name == strongest["method"]).lower(),
            }
        )
    write_csv(out / "METHOD_COMPARISON.csv", comparison_rows)

    paired_rows = [
        {
            "genome_id": gid,
            "full_recovery_percent": 100.0 * recoveries[a.full_method][gid],
            "baseline_recovery_percent": 100.0 * recoveries[strongest["method"]][gid],
            "difference_percentage_points": 100.0
            * (recoveries[a.full_method][gid] - recoveries[strongest["method"]][gid]),
        }
        for gid in common
    ]
    write_csv(out / "PAIRED_PER_GENOME_DIFFERENCES.csv", paired_rows)

    ablation_rows: list[dict[str, Any]] = []
    for label, method in [
        ("full", a.full_method),
        ("no_longitudinal", a.no_longitudinal_method),
        ("no_long", a.no_long_method),
        ("no_consensus", a.no_consensus_method),
    ]:
        if method and method in summaries:
            summary = summaries[method]
            ablation_rows.append(
                {
                    "ablation": label,
                    "method": method,
                    "genome_fraction_percent": summary["genome_fraction_percent"],
                    "mean_genome_recovery_percent": 100.0 * float(summary["mean_genome_recovery"]),
                    "cross_binid_chimeric_bp_fraction": summary["cross_binid_chimeric_bp_fraction"],
                }
            )
    write_csv(out / "ABLATION_RESULTS.csv", ablation_rows)

    lines = [
        "# LANTERN frozen untouched-holdout decision",
        "",
        f"- Status: **{status}**",
        f"- Strict success: **{strict_success}**",
        f"- Holdout samples: **{config['holdout_pair']}**",
        f"- Full method: `{a.full_method}`",
        f"- Strongest baseline: `{strongest['method']}`",
        f"- Genome-fraction gain: **{gf_gain:.6f} percentage points**",
        f"- Mean recovery gain: **{mean_gain:.6f} percentage points**",
        f"- Paired bootstrap mean gain (95% CI): **{paired_mean:.6f} ({ci_low:.6f}, {ci_high:.6f}) percentage points**"
        if paired_mean is not None and ci_low is not None and ci_high is not None
        else "- Paired bootstrap: not evaluable",
        f"- Longitudinal ablation drop: **{longitudinal_drop:.6f} percentage points**",
        f"- Low-abundance gain: **{low_gain} percentage points**",
        f"- Relative chimera-proxy change: **{chimera_change}**",
        "",
        "## Gates",
        "",
    ]
    for key, value in gates.items():
        lines.append(f"- {key}: **{'PASS' if value else 'FAIL'}**")
    lines.extend(["", "## Claim boundary", "", config["claim_boundary"], ""])
    (out / "HOLDOUT_DECISION.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
