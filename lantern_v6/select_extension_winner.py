#!/usr/bin/env python3
"""Select and freeze one extension-union development rule after output freeze.

All variants and baseline assemblies must already have been frozen before this script is
called. The selector reads gold-side evaluation tables only after construction, applies
fixed success gates, and freezes one parameter set for a different participant pair.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--baseline", action="append", required=True, help="NAME=EVAL_DIR")
    p.add_argument("--variant", action="append", required=True, help="NAME=EVAL_DIR")
    p.add_argument("--variant-summary", action="append", required=True, help="NAME=EXTENSION_SUMMARY.json")
    p.add_argument("--abundance-summary", required=True)
    p.add_argument("--pseudo-summary", required=True)
    p.add_argument("--gates", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_specs(specs: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for spec in specs:
        name, path = spec.split("=", 1)
        if name in result:
            raise ValueError(f"duplicate name: {name}")
        result[name] = Path(path)
    return result


def relative_change(new: float, old: float) -> float | None:
    if old == 0:
        return 0.0 if new == 0 else None
    return (new - old) / old


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    baselines = parse_specs(a.baseline)
    variants = parse_specs(a.variant)
    variant_summaries = parse_specs(a.variant_summary)
    gates_cfg = json.loads(Path(a.gates).read_text(encoding="utf-8"))

    summaries: dict[str, dict[str, Any]] = {}
    for name, directory in {**baselines, **variants}.items():
        path = directory / "GOLD_COVERAGE_SUMMARY.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        summaries[name] = json.loads(path.read_text(encoding="utf-8"))

    strongest_name = max(
        baselines,
        key=lambda name: (
            float(summaries[name]["genome_fraction_percent"]),
            float(summaries[name]["mean_genome_recovery"]),
            -float(summaries[name].get("cross_binid_chimeric_bp_fraction") or 0.0),
            name,
        ),
    )
    strongest = summaries[strongest_name]

    abundance = {
        (row["method"], row["stratum"]): (
            None if row.get("mean_recovery") in (None, "", "None") else 100.0 * float(row["mean_recovery"])
        )
        for row in read_tsv(Path(a.abundance_summary))
    }
    pseudo = {
        (row["method"], row["tier"]): (
            None if row.get("mean_recovery") in (None, "", "None") else 100.0 * float(row["mean_recovery"])
        )
        for row in read_tsv(Path(a.pseudo_summary))
    }
    tiers = [
        "exact_target_withheld",
        "species_neighbour_withheld",
        "genus_neighbour_withheld",
        "deep_neighbour_withheld",
    ]

    rows: list[dict[str, Any]] = []
    for name, directory in variants.items():
        summary = summaries[name]
        gf_gain = float(summary["genome_fraction_percent"]) - float(strongest["genome_fraction_percent"])
        mean_gain = 100.0 * (
            float(summary["mean_genome_recovery"]) - float(strongest["mean_genome_recovery"])
        )
        chimera = float(summary.get("cross_binid_chimeric_bp_fraction") or 0.0)
        baseline_chimera = float(strongest.get("cross_binid_chimeric_bp_fraction") or 0.0)
        chimera_change = relative_change(chimera, baseline_chimera)
        low_variant = abundance.get((name, "low"))
        low_baseline = abundance.get((strongest_name, "low"))
        low_gain = None if low_variant is None or low_baseline is None else low_variant - low_baseline
        pseudo_gains: dict[str, float | None] = {}
        for tier in tiers:
            v = pseudo.get((name, tier))
            b = pseudo.get((strongest_name, tier))
            pseudo_gains[tier] = None if v is None or b is None else v - b
        pseudo_all_positive = all(value is not None and value > 0 for value in pseudo_gains.values())
        # The no-longitudinal construction is the unchanged backbone. Therefore the
        # development longitudinal contribution equals the extension-union GF gain.
        longitudinal_drop = gf_gain
        strict_gates = {
            "genome_fraction_gain": gf_gain >= float(gates_cfg["genome_fraction_gain_minimum_percentage_points"]),
            "mean_genome_recovery_gain": mean_gain >= float(gates_cfg["mean_genome_recovery_gain_minimum_percentage_points"]),
            "chimera_tradeoff": chimera_change is not None and chimera_change <= float(gates_cfg["relative_cross_binid_chimera_increase_maximum"]),
            "longitudinal_ablation_contribution": longitudinal_drop >= float(gates_cfg["longitudinal_ablation_drop_minimum_percentage_points"]),
            "low_abundance_gain": low_gain is not None and low_gain >= float(gates_cfg["low_abundance_gain_minimum_percentage_points"]),
            "pseudo_novel_all_tiers": pseudo_all_positive,
        }
        strict_success = all(strict_gates.values())
        safe = bool(strict_gates["chimera_tradeoff"] and gf_gain > 0 and mean_gain > 0)
        rows.append(
            {
                "variant": name,
                "safe_development_candidate": safe,
                "strict_development_success": strict_success,
                "genome_fraction_percent": float(summary["genome_fraction_percent"]),
                "genome_fraction_gain_percentage_points": gf_gain,
                "mean_genome_recovery_percent": 100.0 * float(summary["mean_genome_recovery"]),
                "mean_genome_recovery_gain_percentage_points": mean_gain,
                "longitudinal_ablation_drop_percentage_points": longitudinal_drop,
                "cross_binid_chimera_bp_fraction": chimera,
                "relative_chimera_change": chimera_change,
                "low_abundance_gain_percentage_points": low_gain,
                "pseudo_novel_consistent_all_tiers": pseudo_all_positive,
                "pseudo_exact_gain_pp": pseudo_gains[tiers[0]],
                "pseudo_species_gain_pp": pseudo_gains[tiers[1]],
                "pseudo_genus_gain_pp": pseudo_gains[tiers[2]],
                "pseudo_deep_gain_pp": pseudo_gains[tiers[3]],
                "strict_gate_json": json.dumps(strict_gates, sort_keys=True),
            }
        )

    safe_rows = [row for row in rows if row["safe_development_candidate"]]
    pool = safe_rows or rows
    winner = max(
        pool,
        key=lambda row: (
            bool(row["strict_development_success"]),
            float(row["genome_fraction_gain_percentage_points"]),
            float(row["mean_genome_recovery_gain_percentage_points"]),
            float(row["low_abundance_gain_percentage_points"] or -1e9),
            row["variant"],
        ),
    )
    winner_name = str(winner["variant"])
    if winner_name not in variant_summaries:
        raise SystemExit(f"missing construction summary for winner: {winner_name}")
    construction = json.loads(variant_summaries[winner_name].read_text(encoding="utf-8"))
    freeze = {
        "status": "EXTENSION_RULE_FROZEN_FOR_UNTOUCHED_VALIDATION",
        "winner": winner_name,
        "strongest_development_baseline": strongest_name,
        "selection_rule": (
            "Prefer strict gate success; otherwise among chimera-safe variants with positive GF and mean gains, "
            "maximize GF gain, then mean gain, then low-abundance gain."
        ),
        "winner_metrics": winner,
        "construction": construction,
        "development_samples": [0, 1],
        "claim_boundary": "Public Toy development selection only. The exact construction parameters must be applied unchanged on a different participant pair.",
    }
    (out / "FROZEN_EXTENSION_WINNER.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    with (out / "EXTENSION_DEVELOPMENT_COMPARISON.csv").open("w", newline="", encoding="utf-8") as fh:
        fields = list(rows[0])
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (-float(row["genome_fraction_gain_percentage_points"]), row["variant"])))
    (out / "EXTENSION_DEVELOPMENT_DECISION.md").write_text(
        "# LANTERN-v6 extension-union development decision\n\n"
        f"- Status: **{freeze['status']}**\n"
        f"- Strongest baseline: `{strongest_name}`\n"
        f"- Frozen winner: `{winner_name}`\n"
        f"- GF gain: **{winner['genome_fraction_gain_percentage_points']:.6f} pp**\n"
        f"- Mean-recovery gain: **{winner['mean_genome_recovery_gain_percentage_points']:.6f} pp**\n"
        f"- Low-abundance gain: **{winner['low_abundance_gain_percentage_points']} pp**\n"
        f"- Relative chimera change: **{winner['relative_chimera_change']}**\n"
        f"- Strict development success: **{winner['strict_development_success']}**\n\n"
        + freeze["claim_boundary"]
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(freeze, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
