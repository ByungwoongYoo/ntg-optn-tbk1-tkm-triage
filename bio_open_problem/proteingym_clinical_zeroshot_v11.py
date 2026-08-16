#!/usr/bin/env python3
"""Nested ProteinGym clinical benchmark restricted to unsupervised/zero-shot protein models."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from proteingym_clinical_rank_v9 import (
    bootstrap_ci, compare, load_scores, orientation_performance, protein_auc,
)
from proteingym_clinical_nested_v10 import inner_tune, outer_fold

ALLOWED = ["PoET", "TranceptEVE_L", "GEMME", "EVE", "ESM1b"]


@dataclass(frozen=True)
class Spec:
    name: str
    kind: str
    k: int
    power: float = 1.0


SPECS = [
    Spec("mean2", "mean", 2), Spec("mean3", "mean", 3), Spec("mean5", "mean", 5),
    Spec("median3", "median", 3), Spec("median5", "median", 5),
    Spec("weighted3_p1", "weighted", 3, 1), Spec("weighted5_p1", "weighted", 5, 1),
    Spec("weighted3_p2", "weighted", 3, 2), Spec("weighted5_p2", "weighted", 5, 2),
    Spec("weighted3_p4", "weighted", 3, 4), Spec("weighted5_p4", "weighted", 5, 4),
]


def selected_models(perf: dict[str, float], k: int) -> list[str]:
    return sorted(perf, key=lambda m: (-np.nan_to_num(perf[m], nan=-9), m))[:min(k, len(perf))]


def predict(df: pd.DataFrame, spec: Spec, signs: dict[str, int], perf: dict[str, float]) -> np.ndarray:
    models = selected_models(perf, spec.k)
    x = df[models].to_numpy(dtype=float, copy=True)
    for j, m in enumerate(models):
        if signs[m] < 0:
            x[:, j] = 1.0 - x[:, j]
    if spec.kind == "mean":
        return np.nanmean(x, axis=1)
    if spec.kind == "median":
        return np.nanmedian(x, axis=1)
    w = np.array([max(perf[m] - .5, .001) ** spec.power for m in models])
    w /= w.sum()
    ok = np.isfinite(x)
    ww = ok * w.reshape(1, -1)
    den = ww.sum(axis=1)
    out = np.full(len(df), np.nan)
    valid = den > 0
    out[valid] = np.nansum(x[valid] * w.reshape(1, -1), axis=1) / den[valid]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical-score-zip", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    data, all_models, diagnostics = load_scores(Path(args.clinical_score_zip))
    models = [m for m in ALLOWED if m in all_models]
    if len(models) < 4:
        raise RuntimeError(f"Only {len(models)} allowed zero-shot columns found: {models}")
    data = data[["protein_file", "sequence_hash", "mutant", "label"] + models].copy()
    data["outer_fold"] = data["sequence_hash"].map(outer_fold)

    fold_records = []
    protein_scores = []
    individual_scores = []
    variant_predictions = []

    for fold in range(5):
        outer_train = data[data["outer_fold"] != fold].copy()
        outer_test = data[data["outer_fold"] == fold].copy()
        outer_train["inner_tune"] = outer_train["sequence_hash"].map(inner_tune)
        inner_train = outer_train[~outer_train["inner_tune"]].copy()
        inner_val = outer_train[outer_train["inner_tune"]].copy()

        signs_i, perf_i = orientation_performance(inner_train, models)
        tuning = []
        for spec in SPECS:
            p = predict(inner_val, spec, signs_i, perf_i)
            t = protein_auc(inner_val, p, spec.name)
            tuning.append({"method": spec.name, "mean_auc": float(t[spec.name].mean())})
        tune = pd.DataFrame(tuning).sort_values(["mean_auc", "method"], ascending=[False, True])
        name = str(tune.iloc[0]["method"])
        spec = next(s for s in SPECS if s.name == name)

        signs, perf = orientation_performance(outer_train, models)
        ep = predict(outer_test, spec, signs, perf)
        ens = protein_auc(outer_test, ep, "ensemble")
        best = max(perf, key=lambda m: np.nan_to_num(perf[m], nan=-9))
        bp = outer_test[best].to_numpy(float, copy=True)
        if signs[best] < 0:
            bp = 1.0 - bp
        base = protein_auc(outer_test, bp, "best_selected")
        comp = ens.merge(base[["protein_file", "best_selected"]], on="protein_file")
        comp["outer_fold"] = fold
        comp["method"] = name
        comp["best_selected_model"] = best
        protein_scores.append(comp)

        for m in models:
            q = outer_test[m].to_numpy(float, copy=True)
            if signs[m] < 0:
                q = 1.0 - q
            mt = protein_auc(outer_test, q, "auc")
            mt["model"] = m
            mt["outer_fold"] = fold
            individual_scores.append(mt)

        variant_predictions.append(pd.DataFrame({
            "protein_file": outer_test["protein_file"], "sequence_hash": outer_test["sequence_hash"],
            "mutant": outer_test["mutant"], "label": outer_test["label"],
            "outer_fold": fold, "ensemble_prediction": ep,
        }))
        fold_records.append({
            "fold": fold, "method": name,
            "constituents": selected_models(perf, spec.k),
            "best_selected": best,
            "outer_proteins": int(outer_test["protein_file"].nunique()),
            "ensemble_mean_auc": float(ens["ensemble"].mean()),
            "baseline_mean_auc": float(base["best_selected"].mean()),
            "inner_top3": tune.head(3).to_dict(orient="records"),
        })
        print(json.dumps(fold_records[-1], indent=2), flush=True)

    comp = pd.concat(protein_scores, ignore_index=True)
    indiv = pd.concat(individual_scores, ignore_index=True)
    variants = pd.concat(variant_predictions, ignore_index=True)
    vs_selected = compare(comp, "ensemble", "best_selected")
    summary = indiv.groupby("model")["auc"].agg(["mean", "median", "count"]).sort_values("mean", ascending=False)
    oracle = str(summary.index[0])
    oracle_tab = indiv[indiv["model"] == oracle][["protein_file", "auc"]].rename(columns={"auc": "oracle_best"})
    comp2 = comp.merge(oracle_tab, on="protein_file", how="inner")
    vs_oracle = compare(comp2, "ensemble", "oracle_best")

    result = {
        "question": "Can a nested ensemble of five unsupervised zero-shot protein models improve clinical missense AUC?",
        "allowed_models": models,
        "folds": fold_records,
        "proteins": int(comp["protein_file"].nunique()),
        "sequence_groups": int(variants["sequence_hash"].nunique()),
        "variants": int(len(variants)),
        "vs_nested_selected_individual": vs_selected,
        "posthoc_best_zero_shot_individual": oracle,
        "vs_posthoc_oracle": vs_oracle,
        "zero_shot_ensemble_advance_confirmed": bool(vs_selected["gain_95ci"][0] is not None and vs_selected["gain_95ci"][0] > 0),
        "beats_posthoc_oracle_confirmed": bool(vs_oracle["gain_95ci"][0] is not None and vs_oracle["gain_95ci"][0] > 0),
        "open_problem_fully_solved": False,
        "diagnostics": diagnostics,
    }
    comp.to_csv(out / "nested_zero_shot_protein_auc.csv", index=False)
    summary.to_csv(out / "zero_shot_individual_summary.csv")
    variants.to_csv(out / "nested_zero_shot_predictions.csv.gz", index=False, compression="gzip")
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "REPORT.md").write_text(f"""# ProteinGym zero-shot clinical ensemble v11

Five exact-sequence-grouped outer folds and independent inner model-selection splits were used. Only PoET, TranceptEVE_L, GEMME, EVE, and ESM1b were eligible; supervised clinical meta-predictors were excluded.

- Proteins: {result['proteins']}
- Variants: {result['variants']:,}
- Ensemble mean protein AUC: {vs_selected['ensemble_mean_auc']:.4f}
- Nested-selected zero-shot individual: {vs_selected['comparator_mean_auc']:.4f}
- Gain: {vs_selected['mean_gain']:+.4f}; 95% CI {vs_selected['gain_95ci']}
- Post-hoc best zero-shot individual: `{oracle}`
- Gain versus post-hoc oracle: {vs_oracle['mean_gain']:+.4f}; 95% CI {vs_oracle['gain_95ci']}
- Advance confirmed: {result['zero_shot_ensemble_advance_confirmed']}

This excludes clinical-label-trained component predictors, but the underlying public scores may still inherit model-development and benchmark overlap. It is not a prospective clinical test.
""", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
