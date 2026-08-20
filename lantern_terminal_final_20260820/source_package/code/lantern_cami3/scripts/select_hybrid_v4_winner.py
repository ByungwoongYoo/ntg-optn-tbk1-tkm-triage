#!/usr/bin/env python3
"""Select one LANTERN-v4 development configuration under frozen safety rules.

This script is intended only for the predeclared public Toy development pair. It selects
one configuration, freezes its exact JSON and hashes, and explicitly prohibits threshold
changes before the untouched corrected Toy holdout is opened.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--grid", required=True)
    p.add_argument("--baseline-evaluation", required=True)
    p.add_argument("--full", action="append", required=True, help="VARIANT=EVALUATION_DIR")
    p.add_argument("--no-longitudinal", action="append", required=True, help="VARIANT=EVALUATION_DIR")
    p.add_argument("--abundance-summary", required=True)
    p.add_argument("--pseudo-summary", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def specs(values: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for value in values:
        name, path = value.split("=", 1)
        if name in out:
            raise ValueError(f"duplicate variant {name}")
        out[name] = Path(path)
    return out


def load_score(directory: Path) -> dict[str, Any]:
    path = directory / "GOLD_COVERAGE_SUMMARY.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def load_abundance(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    return {(r["method"], r["stratum"]): r for r in rows}


def abundance_percent(table: dict[tuple[str, str], dict[str, str]], method: str, stratum: str) -> float | None:
    row = table.get((method, stratum))
    if not row or row.get("mean_recovery") in (None, "", "None"):
        return None
    return 100.0 * float(row["mean_recovery"])


def load_pseudo(path: Path) -> dict[tuple[str, str], float]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    return {(r["method"], r["tier"]): 100.0 * float(r["mean_recovery"]) for r in rows}


def relative_change(new: float, old: float) -> float | None:
    if old == 0:
        return 0.0 if new == 0 else None
    return (new - old) / old


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    grid_path = Path(a.grid)
    grid = json.loads(grid_path.read_text())
    full_dirs = specs(a.full)
    no_long_dirs = specs(a.no_longitudinal)
    variant_names = list(grid["variants"])
    if set(full_dirs) != set(variant_names) or set(no_long_dirs) != set(variant_names):
        raise ValueError("evaluated variant names differ from frozen grid")

    baseline_name = grid["backbone"]["method"]
    baseline = load_score(Path(a.baseline_evaluation))
    abundance = load_abundance(Path(a.abundance_summary))
    pseudo = load_pseudo(Path(a.pseudo_summary))
    low_baseline = abundance_percent(abundance, baseline_name, "low")
    pseudo_tiers = [
        "exact_target_withheld",
        "species_neighbour_withheld",
        "genus_neighbour_withheld",
        "deep_neighbour_withheld",
    ]
    safety = grid["development_selection_rule"]["safety_constraints"]
    strict = grid["strict_success_gates"]

    rows: list[dict[str, Any]] = []
    for name in variant_names:
        score = load_score(full_dirs[name])
        no_long = load_score(no_long_dirs[name])
        gf_gain = float(score["genome_fraction_percent"]) - float(baseline["genome_fraction_percent"])
        mean_gain = 100.0 * (float(score["mean_genome_recovery"]) - float(baseline["mean_genome_recovery"]))
        ablation_drop = float(score["genome_fraction_percent"]) - float(no_long["genome_fraction_percent"])
        chimera_change = relative_change(
            float(score.get("cross_binid_chimeric_bp_fraction") or 0.0),
            float(baseline.get("cross_binid_chimeric_bp_fraction") or 0.0),
        )
        alignment_increase = float(score.get("alignment_to_unique_truth_ratio") or 0.0) - float(
            baseline.get("alignment_to_unique_truth_ratio") or 0.0
        )
        low_full = abundance_percent(abundance, name, "low")
        low_gain = None if low_full is None or low_baseline is None else low_full - low_baseline
        pseudo_gains = {
            tier: None
            if (name, tier) not in pseudo or (baseline_name, tier) not in pseudo
            else pseudo[(name, tier)] - pseudo[(baseline_name, tier)]
            for tier in pseudo_tiers
        }
        pseudo_consistent = all(value is not None and value > 0 for value in pseudo_gains.values())
        safe = (
            gf_gain > 0
            and chimera_change is not None
            and chimera_change <= float(safety["relative_cross_binid_chimera_increase_maximum"])
            and alignment_increase <= float(safety["alignment_to_unique_truth_ratio_increase_maximum"])
            and ablation_drop >= float(safety["longitudinal_ablation_drop_minimum_percentage_points"])
        )
        strict_gates = {
            "genome_fraction_gain": gf_gain >= float(strict["genome_fraction_gain_minimum_percentage_points"]),
            "mean_genome_recovery_gain": mean_gain >= float(strict["mean_genome_recovery_gain_minimum_percentage_points"]),
            "chimera_tradeoff": chimera_change is not None
            and chimera_change <= float(strict["relative_cross_binid_chimera_increase_maximum"]),
            "longitudinal_ablation_contribution": ablation_drop
            >= float(strict["longitudinal_ablation_drop_minimum_percentage_points"]),
            "low_abundance_gain": low_gain is not None
            and low_gain >= float(strict["low_abundance_gain_minimum_percentage_points"]),
            "pseudo_novel_all_tiers": pseudo_consistent,
        }
        rows.append(
            {
                "variant": name,
                "safe_development_candidate": safe,
                "strict_development_success": all(strict_gates.values()),
                "genome_fraction_percent": float(score["genome_fraction_percent"]),
                "genome_fraction_gain_percentage_points": gf_gain,
                "mean_genome_recovery_percent": 100.0 * float(score["mean_genome_recovery"]),
                "mean_genome_recovery_gain_percentage_points": mean_gain,
                "longitudinal_ablation_drop_percentage_points": ablation_drop,
                "cross_binid_chimera_bp_fraction": float(score.get("cross_binid_chimeric_bp_fraction") or 0.0),
                "relative_chimera_change": chimera_change,
                "alignment_to_unique_truth_ratio": float(score.get("alignment_to_unique_truth_ratio") or 0.0),
                "alignment_ratio_increase": alignment_increase,
                "low_abundance_mean_recovery_percent": low_full,
                "low_abundance_gain_percentage_points": low_gain,
                "pseudo_novel_consistent_all_tiers": pseudo_consistent,
                "pseudo_exact_gain_pp": pseudo_gains["exact_target_withheld"],
                "pseudo_species_gain_pp": pseudo_gains["species_neighbour_withheld"],
                "pseudo_genus_gain_pp": pseudo_gains["genus_neighbour_withheld"],
                "pseudo_deep_gain_pp": pseudo_gains["deep_neighbour_withheld"],
                "strict_gate_json": json.dumps(strict_gates, sort_keys=True),
            }
        )

    safe_rows = [r for r in rows if r["safe_development_candidate"]]
    winner = max(
        safe_rows,
        key=lambda r: (
            r["genome_fraction_gain_percentage_points"],
            r["mean_genome_recovery_gain_percentage_points"],
            r["low_abundance_gain_percentage_points"] if r["low_abundance_gain_percentage_points"] is not None else float("-inf"),
            -float(r["relative_chimera_change"]),
        ),
    ) if safe_rows else None

    fields = list(rows[0])
    with (out / "V4_DEVELOPMENT_COMPARISON.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (not r["safe_development_candidate"], -r["genome_fraction_gain_percentage_points"])))

    freeze: dict[str, Any] = {
        "status": "WINNER_FROZEN_FOR_UNTOUCHED_VALIDATION" if winner else "NO_SAFE_V4_DEVELOPMENT_WINNER",
        "development_pair": grid["development_pair"],
        "untouched_validation_pair": grid["untouched_validation_pair"],
        "baseline_method": baseline_name,
        "winner": winner,
        "selection_rule": grid["development_selection_rule"],
        "strict_success_gates": grid["strict_success_gates"],
        "grid_sha256": sha256(grid_path),
        "claim_boundary": "The winner, if any, is selected on the public Toy development pair only. It must be rerun unchanged on the corrected untouched Toy pair 14/16 before any performance claim.",
    }
    if winner:
        winner_name = winner["variant"]
        winner_config = grid["variants"][winner_name]
        winner_config_path = out / "FROZEN_V4_WINNER_CONFIG.json"
        winner_config_path.write_text(json.dumps(winner_config, indent=2, sort_keys=True) + "\n")
        freeze["winner_config_sha256"] = sha256(winner_config_path)
        freeze["winner_name"] = winner_name
    (out / "V4_DEVELOPMENT_FREEZE.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")

    lines = [
        "# LANTERN v4 hybrid-backbone development decision",
        "",
        f"- Status: **{freeze['status']}**",
        f"- Baseline: `{baseline_name}`",
    ]
    if winner:
        lines += [
            f"- Frozen winner: `{winner['variant']}`",
            f"- Genome-fraction gain: **{winner['genome_fraction_gain_percentage_points']:.6f} percentage points**",
            f"- Mean-recovery gain: **{winner['mean_genome_recovery_gain_percentage_points']:.6f} percentage points**",
            f"- Low-abundance gain: **{winner['low_abundance_gain_percentage_points']} percentage points**",
            f"- Longitudinal ablation drop: **{winner['longitudinal_ablation_drop_percentage_points']:.6f} percentage points**",
            f"- Relative chimera change: **{winner['relative_chimera_change']:.6f}**",
            "",
            "This result is development-only. The exact config and code must be used without retuning on corrected holdout pair 14/16.",
        ]
    else:
        lines += ["", "No configuration satisfied the frozen positive-gain and safety constraints; v4 cannot be promoted to holdout validation."]
    (out / "V4_DEVELOPMENT_DECISION.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(freeze, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
