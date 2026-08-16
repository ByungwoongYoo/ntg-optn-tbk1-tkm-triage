#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from proteingym_wide_v2 import (
    SEED,
    Spec,
    bootstrap_mean_ci,
    fit_method,
    load_data,
    model_train_stats,
    per_assay_scores,
    safe_spearman,
    stable_bucket,
)


def summarize_assays(compare: pd.DataFrame) -> dict:
    d = (compare["ensemble"] - compare["best_base"]).to_numpy(float)
    d = d[np.isfinite(d)]
    try:
        p = float(wilcoxon(d, alternative="greater").pvalue) if len(d) >= 5 and np.any(d != 0) else None
    except Exception:
        p = None
    return {
        "n_assays": int(len(d)),
        "ensemble_mean_assay_spearman": float(compare["ensemble"].mean()),
        "best_base_mean_assay_spearman": float(compare["best_base"].mean()),
        "paired_mean_improvement": float(np.mean(d)),
        "paired_improvement_95ci": bootstrap_mean_ci(d),
        "fraction_assays_improved": float(np.mean(d > 0)),
        "wilcoxon_one_sided_p": p,
    }


def summarize_proteins(compare: pd.DataFrame) -> dict:
    p = compare.groupby("UniProt_ID")[["ensemble", "best_base"]].mean()
    d = (p["ensemble"] - p["best_base"]).to_numpy(float)
    d = d[np.isfinite(d)]
    try:
        pv = float(wilcoxon(d, alternative="greater").pvalue) if len(d) >= 5 and np.any(d != 0) else None
    except Exception:
        pv = None
    return {
        "n_proteins": int(len(d)),
        "ensemble_mean_protein_averaged_spearman": float(p["ensemble"].mean()),
        "best_base_mean_protein_averaged_spearman": float(p["best_base"].mean()),
        "paired_mean_improvement": float(np.mean(d)),
        "paired_improvement_95ci": bootstrap_mean_ci(d),
        "fraction_proteins_improved": float(np.mean(d > 0)),
        "wilcoxon_one_sided_p": pv,
    }


def evaluate_subset(df: pd.DataFrame, fitted, best_base: str, base_sign: int) -> tuple[pd.DataFrame, dict]:
    pred = fitted.predict(df)
    ens = per_assay_scores(df, pred, "ensemble")
    b = pd.to_numeric(df[best_base], errors="coerce").to_numpy(float)
    if base_sign < 0:
        b = 1.0 - b
    bas = per_assay_scores(df, b, "best_base")
    comp = ens.merge(bas[["DMS_id", "best_base"]], on="DMS_id", how="inner")
    return comp, {"assay_level": summarize_assays(comp), "protein_level": summarize_proteins(comp)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-zip", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    data, models, diagnostics = load_data(Path(args.score_zip), Path(args.reference), human_only=False)
    human = data[data["taxon"].astype(str).str.lower().eq("human")].copy()
    human_proteins = sorted(human["UniProt_ID"].unique())
    bucket = {p: stable_bucket(p) for p in human_proteins}
    human["split"] = human["UniProt_ID"].map(lambda x: "train" if bucket[x] <= 5 else ("tune" if bucket[x] <= 7 else "final"))
    split_counts = human[["DMS_id", "UniProt_ID", "split"]].drop_duplicates().groupby("split").agg(
        assays=("DMS_id", "nunique"), proteins=("UniProt_ID", "nunique")
    ).to_dict("index")
    # Match v2 fallback rule exactly if needed.
    if any(split_counts.get(s, {}).get("proteins", 0) < 6 for s in ["train", "tune", "final"]):
        assignments = {}
        for i, p in enumerate(human_proteins):
            frac = i / max(1, len(human_proteins))
            assignments[p] = "train" if frac < 0.6 else ("tune" if frac < 0.8 else "final")
        human["split"] = human["UniProt_ID"].map(assignments)

    human_train_tune = human[human["split"].isin(["train", "tune"])].copy()
    human_final = human[human["split"].eq("final")].copy()
    nonhuman = data[~data["taxon"].astype(str).str.lower().eq("human")].copy()

    frozen_spec = Spec("weighted20_p2", "weighted", 20, 2.0)
    fitted = fit_method(human_train_tune, models, frozen_spec)
    base_signs, base_perf = model_train_stats(human_train_tune, models)
    best_base = max(base_perf, key=lambda m: np.nan_to_num(base_perf[m], nan=-9))

    results = {
        "method_frozen_before_external_labels": frozen_spec.name,
        "training_scope": {
            "human_train_tune_assays": int(human_train_tune["DMS_id"].nunique()),
            "human_train_tune_proteins": int(human_train_tune["UniProt_ID"].nunique()),
        },
        "selected_models": fitted.selected,
        "best_base_selected_on_human_train_tune": best_base,
        "diagnostics": diagnostics,
        "subsets": {},
    }

    comp_h, sum_h = evaluate_subset(human_final, fitted, best_base, base_signs[best_base])
    comp_h["evaluation_subset"] = "human_sealed_final"
    results["subsets"]["human_sealed_final"] = sum_h

    comp_n, sum_n = evaluate_subset(nonhuman, fitted, best_base, base_signs[best_base])
    comp_n["evaluation_subset"] = "all_nonhuman_external"
    results["subsets"]["all_nonhuman_external"] = sum_n

    taxa_frames = []
    for taxon, g in nonhuman.groupby("taxon"):
        comp, summary = evaluate_subset(g, fitted, best_base, base_signs[best_base])
        comp["evaluation_subset"] = f"taxon_{taxon}"
        taxa_frames.append(comp)
        results["subsets"][f"taxon_{taxon}"] = summary

    # Combined external set excludes every human assay used to fit weights.
    external = pd.concat([human_final, nonhuman], ignore_index=True)
    comp_e, sum_e = evaluate_subset(external, fitted, best_base, base_signs[best_base])
    comp_e["evaluation_subset"] = "combined_external"
    results["subsets"]["combined_external"] = sum_e

    primary_ci = sum_n["protein_level"]["paired_improvement_95ci"]
    results["external_generalization_confirmed"] = bool(primary_ci[0] is not None and primary_ci[0] > 0)
    results["open_problem_fully_solved"] = False

    all_comp = pd.concat([comp_h, comp_n, comp_e] + taxa_frames, ignore_index=True)
    all_comp.to_csv(out / "external_assay_scores.csv", index=False)
    (out / "result.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")

    nsum = results["subsets"]["all_nonhuman_external"]
    report = f"""# ProteinGym frozen human-trained ensemble: external taxonomic validation

## Frozen rule

The exact `weighted20_p2` rule selected in the preceding human train/tuning experiment was retained. No non-human DMS label was used to choose constituent models, signs, weights, or the comparator.

## External result

- Non-human assays: {nsum['assay_level']['n_assays']}
- Non-human proteins: {nsum['protein_level']['n_proteins']}
- Ensemble mean assay Spearman: {nsum['assay_level']['ensemble_mean_assay_spearman']:.4f}
- Frozen best individual (`{best_base}`): {nsum['assay_level']['best_base_mean_assay_spearman']:.4f}
- Assay-level paired improvement: {nsum['assay_level']['paired_mean_improvement']:+.4f}
- Assay-level 95% CI: {nsum['assay_level']['paired_improvement_95ci']}
- Protein-level paired improvement: {nsum['protein_level']['paired_mean_improvement']:+.4f}
- Protein-level 95% CI: {nsum['protein_level']['paired_improvement_95ci']}

`external_generalization_confirmed = {results['external_generalization_confirmed']}`

A positive result would validate a cross-protein supervised ensemble strategy beyond humans. It would remain a benchmark advance rather than a complete solution to human missense interpretation.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
