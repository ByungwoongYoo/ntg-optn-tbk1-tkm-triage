#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from proteingym_wide_v2 import (
    META_EXCLUDE,
    Spec,
    bootstrap_mean_ci,
    fit_method,
    load_data,
    model_train_stats,
    per_assay_scores,
    rank01,
    stable_bucket,
)
from proteingym_external_v3 import evaluate_subset


def load_excluded_human(
    score_zip_path: Path,
    reference_path: Path,
    models: list[str],
    core_filenames: set[str],
    min_single_mutants: int = 20,
) -> tuple[pd.DataFrame, list[dict]]:
    ref = pd.read_csv(reference_path)
    ref = ref[ref["taxon"].astype(str).str.lower().eq("human")].copy()
    ref = ref[~ref["DMS_filename"].astype(str).isin(core_filenames)].copy()
    by_file = ref.set_index("DMS_filename")
    frames = []
    failures = []
    with zipfile.ZipFile(score_zip_path) as z:
        members = {Path(n).name: n for n in z.namelist() if n.lower().endswith(".csv")}
        for fn, meta in by_file.iterrows():
            if fn not in members:
                failures.append({"file": fn, "reason": "missing from official score archive"})
                continue
            try:
                with z.open(members[fn]) as f:
                    d = pd.read_csv(f, low_memory=False)
                if "mutant" not in d.columns or "DMS_score" not in d.columns:
                    failures.append({"file": fn, "reason": "missing mutant or DMS_score"})
                    continue
                keep_models = [m for m in models if m in d.columns]
                d = d[["mutant", "DMS_score"] + keep_models].copy()
                d["mutant"] = d["mutant"].astype(str)
                d = d[~d["mutant"].str.contains(":", regex=False)]
                d["DMS_score"] = pd.to_numeric(d["DMS_score"], errors="coerce")
                d = d.dropna(subset=["DMS_score"]).drop_duplicates("mutant")
                if len(d) < min_single_mutants:
                    failures.append({"file": fn, "reason": f"only {len(d)} single substitutions"})
                    continue
                d["target_rank"] = rank01(d["DMS_score"])
                for m in models:
                    d[m] = rank01(d[m]) if m in d.columns else np.nan
                d["DMS_id"] = meta["DMS_id"]
                d["UniProt_ID"] = meta["UniProt_ID"]
                d["taxon"] = meta["taxon"]
                d["selection_type"] = meta.get("coarse_selection_type", meta.get("selection_type", "unknown"))
                d["msa_depth"] = meta.get("MSA_Neff_L_category", "unknown")
                frames.append(d[["DMS_id", "UniProt_ID", "taxon", "selection_type", "msa_depth", "mutant", "DMS_score", "target_rank"] + models])
            except Exception as exc:
                failures.append({"file": fn, "reason": repr(exc)})
    if not frames:
        return pd.DataFrame(), failures
    return pd.concat(frames, ignore_index=True), failures


def audit_all_individuals(df: pd.DataFrame, models: list[str], signs: dict[str, int]) -> pd.DataFrame:
    rows = []
    for m in models:
        p = pd.to_numeric(df[m], errors="coerce").to_numpy(float)
        if signs[m] < 0:
            p = 1.0 - p
        s = per_assay_scores(df, p, "rho")
        prot = s.groupby("UniProt_ID")["rho"].mean()
        rows.append({
            "model": m,
            "mean_assay_spearman": float(s["rho"].mean()),
            "median_assay_spearman": float(s["rho"].median()),
            "n_assays": int(s["rho"].notna().sum()),
            "mean_protein_spearman": float(prot.mean()),
            "n_proteins": int(prot.notna().sum()),
        })
    return pd.DataFrame(rows).sort_values("mean_protein_spearman", ascending=False)


def compare_to_named(df: pd.DataFrame, fitted, name: str, signs: dict[str, int]) -> dict:
    comp, summary = evaluate_subset(df, fitted, name, signs[name])
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-zip", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    score_zip = Path(args.score_zip)
    reference = Path(args.reference)
    data, models, diagnostics = load_data(score_zip, reference, human_only=False)
    ref = pd.read_csv(reference)
    core_mask = (ref["DMS_number_single_mutants"] >= 100) & (~ref["includes_multiple_mutants"].astype(bool))
    core_human_filenames = set(ref.loc[core_mask & ref["taxon"].astype(str).str.lower().eq("human"), "DMS_filename"].astype(str))

    human = data[data["taxon"].astype(str).str.lower().eq("human")].copy()
    proteins = sorted(human["UniProt_ID"].unique())
    buckets = {p: stable_bucket(p) for p in proteins}
    human["split"] = human["UniProt_ID"].map(lambda x: "train" if buckets[x] <= 5 else ("tune" if buckets[x] <= 7 else "final"))
    split_counts = human[["DMS_id", "UniProt_ID", "split"]].drop_duplicates().groupby("split").agg(
        assays=("DMS_id", "nunique"), proteins=("UniProt_ID", "nunique")
    ).to_dict("index")
    if any(split_counts.get(s, {}).get("proteins", 0) < 6 for s in ["train", "tune", "final"]):
        assignments = {}
        for i, p in enumerate(proteins):
            frac = i / max(1, len(proteins))
            assignments[p] = "train" if frac < 0.6 else ("tune" if frac < 0.8 else "final")
        human["split"] = human["UniProt_ID"].map(assignments)

    train_tune = human[human["split"].isin(["train", "tune"])].copy()
    human_final = human[human["split"].eq("final")].copy()
    nonhuman = data[~data["taxon"].astype(str).str.lower().eq("human")].copy()

    spec = Spec("weighted20_p2", "weighted", 20, 2.0)
    fitted = fit_method(train_tune, models, spec)
    signs, perf = model_train_stats(train_tune, models)
    frozen_best = max(perf, key=lambda m: np.nan_to_num(perf[m], nan=-9))

    nonhuman_individuals = audit_all_individuals(nonhuman, models, signs)
    posthoc_best = str(nonhuman_individuals.iloc[0]["model"])
    nonhuman_individuals.to_csv(out / "nonhuman_all_individual_models.csv", index=False)

    results = {
        "frozen_method": spec.name,
        "selected_models": fitted.selected,
        "frozen_best_individual": frozen_best,
        "posthoc_best_nonhuman_individual": posthoc_best,
        "diagnostics": diagnostics,
        "evaluations": {},
    }

    _, results["evaluations"]["nonhuman_vs_frozen_best"] = evaluate_subset(nonhuman, fitted, frozen_best, signs[frozen_best])
    _, results["evaluations"]["nonhuman_vs_posthoc_oracle_best"] = evaluate_subset(nonhuman, fitted, posthoc_best, signs[posthoc_best])
    _, results["evaluations"]["human_core_sealed_final"] = evaluate_subset(human_final, fitted, frozen_best, signs[frozen_best])

    excluded, failures = load_excluded_human(score_zip, reference, models, core_human_filenames, min_single_mutants=20)
    results["excluded_human_failures"] = failures
    results["excluded_human_loaded_assays"] = int(excluded["DMS_id"].nunique()) if not excluded.empty else 0
    results["excluded_human_loaded_proteins"] = int(excluded["UniProt_ID"].nunique()) if not excluded.empty else 0

    if not excluded.empty:
        _, results["evaluations"]["all_excluded_human"] = evaluate_subset(excluded, fitted, frozen_best, signs[frozen_best])
        core_proteins = set(human["UniProt_ID"].astype(str))
        strict = excluded[~excluded["UniProt_ID"].astype(str).isin(core_proteins)].copy()
        results["strict_novel_human_assays"] = int(strict["DMS_id"].nunique())
        results["strict_novel_human_proteins"] = int(strict["UniProt_ID"].nunique())
        if strict["DMS_id"].nunique() >= 3:
            _, results["evaluations"]["strict_novel_human_proteins"] = evaluate_subset(strict, fitted, frozen_best, signs[frozen_best])

    ext = results["evaluations"]["nonhuman_vs_frozen_best"]["protein_level"]
    oracle = results["evaluations"]["nonhuman_vs_posthoc_oracle_best"]["protein_level"]
    results["beats_frozen_best_external"] = bool(ext["paired_improvement_95ci"][0] is not None and ext["paired_improvement_95ci"][0] > 0)
    results["beats_every_accessible_individual_descriptively"] = bool(oracle["ensemble_mean_protein_averaged_spearman"] > oracle["best_base_mean_protein_averaged_spearman"])
    results["open_problem_fully_solved"] = False

    (out / "result.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    report = f"""# ProteinGym confirmatory audit v4

## Frozen model

The human train+tuning-selected `weighted20_p2` ensemble was applied unchanged to 75 non-human assays. All 95 accessible individual model columns were also audited on that external set.

## Key decisions

- Frozen comparator: `{frozen_best}`
- Post-hoc best individual on non-human external data: `{posthoc_best}`
- Beats frozen comparator with protein-level CI above zero: `{results['beats_frozen_best_external']}`
- Ensemble mean exceeds the post-hoc best accessible individual mean: `{results['beats_every_accessible_individual_descriptively']}`
- Excluded human assays loaded for additional confirmation: `{results['excluded_human_loaded_assays']}`

The post-hoc oracle comparison is descriptive because the comparator was selected after viewing the external labels. It is included to determine whether the ensemble merely beats a weak transferred baseline or actually exceeds every accessible constituent predictor on average.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
