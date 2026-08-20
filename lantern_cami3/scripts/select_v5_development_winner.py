#!/usr/bin/env python3
"""Select one LANTERN-v5 long-read method on a predeclared development pair.

All assemblies and ablations must already be frozen. The script opens gold-side tables,
selects the strongest eligible baseline, applies fixed safety constraints, and freezes at
most one backbone/configuration for untouched validation. It cannot modify thresholds or
assemblies.
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
    p.add_argument("--manifest", required=True)
    p.add_argument("--method", action="append", required=True, help="METHOD=EVALUATION_DIR")
    p.add_argument("--abundance-summary", required=True)
    p.add_argument("--pseudo-summary", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_abundance(path: Path) -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    for row in read_tsv(path):
        if row.get("mean_recovery") not in (None, "", "None"):
            out[(row["method"], row["stratum"])] = 100.0 * float(row["mean_recovery"])
    return out


def load_pseudo(path: Path) -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    for row in read_tsv(path):
        if row.get("mean_recovery") not in (None, "", "None"):
            out[(row["method"], row["tier"])] = 100.0 * float(row["mean_recovery"])
    return out


def relative_change(new: float, old: float) -> float | None:
    if old == 0:
        return 0.0 if new == 0 else None
    return (new - old) / old


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    grid_path = Path(a.grid)
    grid = json.loads(grid_path.read_text(encoding="utf-8"))
    manifest = {row["method"]: row for row in read_tsv(a.manifest)}
    directories: dict[str, Path] = {}
    for spec in a.method:
        name, path = spec.split("=", 1)
        if name in directories:
            raise ValueError(f"duplicate method {name}")
        if name not in manifest:
            raise ValueError(f"method absent from manifest: {name}")
        directories[name] = Path(path)

    summaries: dict[str, dict[str, Any]] = {}
    for name, directory in directories.items():
        path = directory / "GOLD_COVERAGE_SUMMARY.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        value = json.loads(path.read_text(encoding="utf-8"))
        value.update(manifest[name])
        value["method"] = name
        summaries[name] = value

    baselines = [
        row for name, row in summaries.items()
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

    abundance = load_abundance(Path(a.abundance_summary))
    pseudo = load_pseudo(Path(a.pseudo_summary))
    low_baseline = abundance.get((strongest["method"], "low"))
    tiers = [
        "exact_target_withheld",
        "species_neighbour_withheld",
        "genus_neighbour_withheld",
        "deep_neighbour_withheld",
    ]
    safety = grid["development_selection_rule"]["safety_constraints"]
    strict = grid["strict_success_gates"]

    rows: list[dict[str, Any]] = []
    full_methods = [
        name for name in summaries
        if manifest[name].get("role") == "lantern" and manifest[name].get("ablation") == "full"
    ]
    for name in sorted(full_methods):
        no_longitudinal_name = name + "_no_longitudinal"
        if no_longitudinal_name not in summaries:
            raise ValueError(f"missing longitudinal ablation for {name}")
        full = summaries[name]
        no_longitudinal = summaries[no_longitudinal_name]
        gf_gain = float(full["genome_fraction_percent"]) - float(strongest["genome_fraction_percent"])
        mean_gain = 100.0 * (
            float(full["mean_genome_recovery"]) - float(strongest["mean_genome_recovery"])
        )
        ablation_drop = float(full["genome_fraction_percent"]) - float(
            no_longitudinal["genome_fraction_percent"]
        )
        chimera_change = relative_change(
            float(full.get("cross_binid_chimeric_bp_fraction") or 0.0),
            float(strongest.get("cross_binid_chimeric_bp_fraction") or 0.0),
        )
        alignment_increase = float(full.get("alignment_to_unique_truth_ratio") or 0.0) - float(
            strongest.get("alignment_to_unique_truth_ratio") or 0.0
        )
        low_full = abundance.get((name, "low"))
        low_gain = None if low_full is None or low_baseline is None else low_full - low_baseline
        pseudo_gains: dict[str, float | None] = {}
        for tier in tiers:
            full_value = pseudo.get((name, tier))
            baseline_value = pseudo.get((strongest["method"], tier))
            pseudo_gains[tier] = (
                None if full_value is None or baseline_value is None else full_value - baseline_value
            )
        pseudo_all_positive = all(v is not None and v > 0 for v in pseudo_gains.values())
        safe = (
            gf_gain > 0
            and chimera_change is not None
            and chimera_change <= float(safety["relative_cross_binid_chimera_increase_maximum"])
            and alignment_increase <= float(safety["alignment_to_unique_truth_ratio_increase_maximum"])
            and ablation_drop >= float(safety["longitudinal_ablation_drop_minimum_percentage_points"])
        )
        strict_gates = {
            "genome_fraction_gain": gf_gain
            >= float(strict["genome_fraction_gain_minimum_percentage_points"]),
            "mean_genome_recovery_gain": mean_gain
            >= float(strict["mean_genome_recovery_gain_minimum_percentage_points"]),
            "chimera_tradeoff": chimera_change is not None
            and chimera_change <= float(strict["relative_cross_binid_chimera_increase_maximum"]),
            "longitudinal_ablation_contribution": ablation_drop
            >= float(strict["longitudinal_ablation_drop_minimum_percentage_points"]),
            "low_abundance_gain": low_gain is not None
            and low_gain >= float(strict["low_abundance_gain_minimum_percentage_points"]),
            "pseudo_novel_all_tiers": pseudo_all_positive,
        }
        parts = name.split("__", 1)
        if len(parts) != 2:
            raise ValueError(f"v5 full method must be BACKBONE__VARIANT: {name}")
        rows.append(
            {
                "method": name,
                "backbone": parts[0],
                "variant": parts[1],
                "safe_development_candidate": safe,
                "strict_development_success": all(strict_gates.values()),
                "strongest_baseline": strongest["method"],
                "genome_fraction_percent": float(full["genome_fraction_percent"]),
                "genome_fraction_gain_percentage_points": gf_gain,
                "mean_genome_recovery_percent": 100.0 * float(full["mean_genome_recovery"]),
                "mean_genome_recovery_gain_percentage_points": mean_gain,
                "longitudinal_ablation_drop_percentage_points": ablation_drop,
                "cross_binid_chimera_bp_fraction": float(full.get("cross_binid_chimeric_bp_fraction") or 0.0),
                "relative_chimera_change": chimera_change,
                "alignment_to_unique_truth_ratio": float(full.get("alignment_to_unique_truth_ratio") or 0.0),
                "alignment_ratio_increase": alignment_increase,
                "low_abundance_mean_recovery_percent": low_full,
                "low_abundance_gain_percentage_points": low_gain,
                "pseudo_novel_consistent_all_tiers": pseudo_all_positive,
                "pseudo_exact_gain_pp": pseudo_gains["exact_target_withheld"],
                "pseudo_species_gain_pp": pseudo_gains["species_neighbour_withheld"],
                "pseudo_genus_gain_pp": pseudo_gains["genus_neighbour_withheld"],
                "pseudo_deep_gain_pp": pseudo_gains["deep_neighbour_withheld"],
                "strict_gate_json": json.dumps(strict_gates, sort_keys=True),
            }
        )

    safe_rows = [row for row in rows if row["safe_development_candidate"]]
    winner = max(
        safe_rows,
        key=lambda row: (
            row["genome_fraction_gain_percentage_points"],
            row["mean_genome_recovery_gain_percentage_points"],
            row["low_abundance_gain_percentage_points"]
            if row["low_abundance_gain_percentage_points"] is not None
            else float("-inf"),
            -float(row["relative_chimera_change"]),
        ),
    ) if safe_rows else None

    fields = list(rows[0]) if rows else ["method"]
    with (out / "V5_DEVELOPMENT_COMPARISON.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            sorted(
                rows,
                key=lambda row: (
                    not row["safe_development_candidate"],
                    -row["genome_fraction_gain_percentage_points"],
                    row["method"],
                ),
            )
        )

    freeze: dict[str, Any] = {
        "status": "WINNER_FROZEN_FOR_UNTOUCHED_VALIDATION" if winner else "NO_SAFE_V5_DEVELOPMENT_WINNER",
        "development_pair": grid["development_pair"],
        "untouched_validation_pair": grid["untouched_validation_pair"],
        "final_replication_pair": grid["final_replication_pair"],
        "strongest_baseline": strongest["method"],
        "winner": winner,
        "selection_rule": grid["development_selection_rule"],
        "strict_success_gates": strict,
        "grid_sha256": sha256(grid_path),
        "claim_boundary": grid["claim_boundary"],
    }
    if winner:
        method = winner["method"]
        config = grid["variants"][winner["variant"]]
        config_path = out / "FROZEN_V5_WINNER_CONFIG.json"
        config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        freeze["winner_method"] = method
        freeze["winner_backbone"] = winner["backbone"]
        freeze["winner_variant"] = winner["variant"]
        freeze["winner_config_sha256"] = sha256(config_path)
    (out / "V5_DEVELOPMENT_FREEZE.json").write_text(
        json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# LANTERN v5 long-read development decision",
        "",
        f"- Status: **{freeze['status']}**",
        f"- Strongest baseline: `{strongest['method']}`",
    ]
    if winner:
        lines.extend(
            [
                f"- Frozen method: `{winner['method']}`",
                f"- Genome-fraction gain: **{winner['genome_fraction_gain_percentage_points']:.6f} percentage points**",
                f"- Mean-recovery gain: **{winner['mean_genome_recovery_gain_percentage_points']:.6f} percentage points**",
                f"- Low-abundance gain: **{winner['low_abundance_gain_percentage_points']} percentage points**",
                f"- Longitudinal ablation drop: **{winner['longitudinal_ablation_drop_percentage_points']:.6f} percentage points**",
                f"- Relative chimera change: **{winner['relative_chimera_change']:.6f}**",
                "",
                "This result is development-only. The exact method must be rerun without retuning on the untouched validation pair.",
            ]
        )
    else:
        lines.extend(["", "No method satisfied the frozen positive-gain and safety constraints."])
    (out / "V5_DEVELOPMENT_DECISION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(freeze, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
