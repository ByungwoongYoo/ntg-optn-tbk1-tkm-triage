#!/usr/bin/env python3
"""Five-fold nested exact-sequence-grouped ProteinGym clinical rank ensemble."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from proteingym_clinical_rank_v9 import (
    SPECS, bootstrap_ci, compare, load_scores, orientation_performance,
    predict, protein_auc, safe_auc, select_diverse,
)


def outer_fold(seq_hash: str) -> int:
    return int(seq_hash[:12], 16) % 5


def inner_tune(seq_hash: str) -> bool:
    return int(seq_hash[12:24], 16) % 5 == 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical-score-zip", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    data, models, diagnostics = load_scores(Path(args.clinical_score_zip))
    data["outer_fold"] = data["sequence_hash"].map(outer_fold)

    fold_results = []
    all_protein_scores = []
    all_variant_predictions = []
    all_model_outer_scores = []

    for fold in range(5):
        outer_train = data[data["outer_fold"] != fold].copy()
        outer_test = data[data["outer_fold"] == fold].copy()
        outer_train["inner_tune"] = outer_train["sequence_hash"].map(inner_tune)
        inner_train = outer_train[~outer_train["inner_tune"]].copy()
        inner_tune_df = outer_train[outer_train["inner_tune"]].copy()
        if inner_train["sequence_hash"].nunique() < 20 or inner_tune_df["sequence_hash"].nunique() < 5:
            raise RuntimeError(f"Fold {fold}: insufficient inner groups")

        signs_inner, perf_inner = orientation_performance(inner_train, models)
        candidates = []
        for spec in SPECS:
            p = predict(inner_tune_df, spec, models, signs_inner, perf_inner)
            t = protein_auc(inner_tune_df, p, spec.name)
            candidates.append({
                "method": spec.name,
                "mean_inner_tune_auc": float(t[spec.name].mean()),
                "selected_models": select_diverse(perf_inner, spec.k),
            })
        cand = pd.DataFrame(candidates).sort_values(["mean_inner_tune_auc", "method"], ascending=[False, True])
        chosen_name = str(cand.iloc[0]["method"])
        chosen = next(s for s in SPECS if s.name == chosen_name)

        signs, perf = orientation_performance(outer_train, models)
        ensemble_pred = predict(outer_test, chosen, models, signs, perf)
        ens = protein_auc(outer_test, ensemble_pred, "ensemble")

        best_selected = max(perf, key=lambda m: np.nan_to_num(perf[m], nan=-9))
        bp = outer_test[best_selected].to_numpy(float, copy=True)
        if signs[best_selected] < 0:
            bp = 1.0 - bp
        base = protein_auc(outer_test, bp, "best_selected")
        comp = ens.merge(base[["protein_file", "best_selected"]], on="protein_file")
        comp["outer_fold"] = fold
        comp["chosen_method"] = chosen_name
        comp["best_selected_model"] = best_selected
        all_protein_scores.append(comp)

        for model in models:
            q = outer_test[model].to_numpy(float, copy=True)
            if signs[model] < 0:
                q = 1.0 - q
            mt = protein_auc(outer_test, q, "auc")
            mt["model"] = model
            mt["outer_fold"] = fold
            all_model_outer_scores.append(mt)

        all_variant_predictions.append(pd.DataFrame({
            "protein_file": outer_test["protein_file"],
            "sequence_hash": outer_test["sequence_hash"],
            "mutant": outer_test["mutant"],
            "label": outer_test["label"],
            "outer_fold": fold,
            "ensemble_prediction": ensemble_pred,
        }))
        fold_results.append({
            "fold": fold,
            "outer_train_sequence_groups": int(outer_train["sequence_hash"].nunique()),
            "inner_train_sequence_groups": int(inner_train["sequence_hash"].nunique()),
            "inner_tune_sequence_groups": int(inner_tune_df["sequence_hash"].nunique()),
            "outer_test_sequence_groups": int(outer_test["sequence_hash"].nunique()),
            "outer_test_protein_files": int(outer_test["protein_file"].nunique()),
            "chosen_method": chosen_name,
            "chosen_models": select_diverse(perf, chosen.k),
            "best_individual_selected_without_outer_labels": best_selected,
            "inner_tuning_top5": cand.head(5).to_dict(orient="records"),
            "outer_ensemble_mean_auc": float(ens["ensemble"].mean()),
            "outer_selected_baseline_mean_auc": float(base["best_selected"].mean()),
        })
        print(json.dumps(fold_results[-1], indent=2), flush=True)

    comp = pd.concat(all_protein_scores, ignore_index=True)
    model_outer = pd.concat(all_model_outer_scores, ignore_index=True)
    variants = pd.concat(all_variant_predictions, ignore_index=True)

    vs_selected = compare(comp, "ensemble", "best_selected")
    model_summary = model_outer.groupby("model")["auc"].agg(["mean", "median", "count"]).sort_values("mean", ascending=False)
    posthoc_oracle = str(model_summary.index[0])
    oracle = model_outer[model_outer["model"] == posthoc_oracle][["protein_file", "auc"]].rename(columns={"auc": "oracle_best"})
    comp_oracle = comp.merge(oracle, on="protein_file", how="inner")
    vs_oracle = compare(comp_oracle, "ensemble", "oracle_best")

    result = {
        "question": "Does a nested exact-sequence-held-out rank ensemble improve ProteinGym clinical missense classification?",
        "diagnostics": diagnostics,
        "outer_folds": fold_results,
        "protein_files_evaluated": int(comp["protein_file"].nunique()),
        "sequence_groups_evaluated": int(variants["sequence_hash"].nunique()),
        "variants_evaluated": int(len(variants)),
        "vs_nested_selected_individual": vs_selected,
        "posthoc_best_individual_across_outer_predictions": posthoc_oracle,
        "vs_posthoc_oracle": vs_oracle,
        "nested_benchmark_advance_confirmed": bool(vs_selected["gain_95ci"][0] is not None and vs_selected["gain_95ci"][0] > 0),
        "beats_posthoc_oracle_confirmed": bool(vs_oracle["gain_95ci"][0] is not None and vs_oracle["gain_95ci"][0] > 0),
        "open_problem_fully_solved": False,
    }

    comp.to_csv(out / "nested_protein_auc.csv", index=False)
    model_summary.to_csv(out / "nested_individual_model_summary.csv")
    variants.to_csv(out / "nested_variant_predictions.csv.gz", index=False, compression="gzip")
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "REPORT.md").write_text(f"""# ProteinGym clinical nested benchmark v10

Five outer folds were grouped by exact protein-sequence hash. Within every outer training set, a separate sequence-grouped inner split selected the ensemble rule. No outer-fold label contributed to its own method selection, score orientation, constituent-model selection, weights, or comparator.

- Protein files: {result['protein_files_evaluated']}
- Exact sequence groups: {result['sequence_groups_evaluated']}
- Variants: {result['variants_evaluated']:,}
- Ensemble mean protein AUC: {vs_selected['ensemble_mean_auc']:.4f}
- Nested-selected individual mean AUC: {vs_selected['comparator_mean_auc']:.4f}
- Mean paired gain: {vs_selected['mean_gain']:+.4f}
- 95% CI: {vs_selected['gain_95ci']}
- Proteins improved: {100 * vs_selected['fraction_improved']:.1f}%
- One-sided Wilcoxon p: {vs_selected['wilcoxon_one_sided_p']}
- Post-hoc best individual: `{posthoc_oracle}`
- Gain versus post-hoc oracle: {vs_oracle['mean_gain']:+.4f}; 95% CI {vs_oracle['gain_95ci']}
- Nested advance confirmed: {result['nested_benchmark_advance_confirmed']}

This is a nested benchmark result on ProteinGym's current clinical annotations. It is not a prospective ClinVar test and does not solve all human missense effects.
""", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
