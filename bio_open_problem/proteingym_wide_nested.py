#!/usr/bin/env python3
"""Leakage-controlled cross-protein ensemble of public ProteinGym predictors.

The official ProteinGym v1.3 score archive is wide: one CSV per DMS assay, with
DMS_score plus many zero-shot model columns. This script evaluates whether a
predeclared collection of public predictors can be combined to improve ranking
of single amino-acid substitutions in entirely held-out human proteins.

Primary safeguards
------------------
1. Outer folds are grouped by UniProt ID.
2. Ensemble-method and baseline-model selection occurs only inside outer
   training proteins through inner grouped cross-validation.
3. The primary score follows ProteinGym's hierarchy: assay -> UniProt x
   selection type -> selection type -> overall mean.
4. All model inputs are within-assay percentile ranks. ProteinGym's public score
   archive already uses the convention that higher predicted score is better.
5. No held-out protein label is used for orientation, weighting, fitting, or
   method selection.
"""
from __future__ import annotations

import argparse
import json
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline

SEED = 20260817
RNG = np.random.default_rng(SEED)

# Fixed before looking at human-assay outcomes. The list deliberately mixes
# evolutionary, sequence-language, structure-aware, retrieval, and inverse-
# folding predictors rather than taking many variants of one model family.
PREFERRED_MODELS = [
    "VenusREM",
    "ProSST-2048",
    "S3F_MSA",
    "PoET",
    "ESM3",
    "VespaG",
    "SaProt_650M_AF2",
    "TranceptEVE_L",
    "GEMME",
    "ProteinMPNN",
    "ESM1v_ensemble",
    "MSA_Transformer_ensemble",
    "ProtSSN_ensemble",
    "ESM2_15B",
]

CANDIDATE_METHODS = [
    "mean_all",
    "median_all",
    "performance_weighted",
    "top3_mean",
    "top5_mean",
    "ridge_0.1",
    "ridge_1",
    "ridge_10",
    "ridge_100",
    "ridge_1000",
]

META_COLUMNS = {
    "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin",
    "DMS_score_bin_manual", "fitness", "score", "target_seq", "sequence",
}


def boolish(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def safe_spearman(x: Iterable[float], y: Iterable[float]) -> float:
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if int(keep.sum()) < 20:
        return float("nan")
    if np.nanstd(a[keep]) == 0 or np.nanstd(b[keep]) == 0:
        return float("nan")
    return float(spearmanr(a[keep], b[keep]).statistic)


def archive_member_map(archive: zipfile.ZipFile) -> dict[str, str]:
    return {
        Path(name).name: name
        for name in archive.namelist()
        if name.lower().endswith(".csv") and not name.endswith("/")
    }


def assay_member(meta: pd.Series, members: dict[str, str]) -> str | None:
    candidates = [
        str(meta.get("DMS_filename", "")),
        f"{meta['DMS_id']}.csv",
    ]
    for candidate in candidates:
        if candidate in members:
            return members[candidate]
    return None


def choose_available_models(
    archive: zipfile.ZipFile,
    ref: pd.DataFrame,
    members: dict[str, str],
    min_assay_coverage: float,
) -> tuple[list[str], dict[str, int]]:
    coverage = {model: 0 for model in PREFERRED_MODELS}
    n_seen = 0
    for _, meta in ref.iterrows():
        member = assay_member(meta, members)
        if member is None:
            continue
        with archive.open(member) as handle:
            columns = set(pd.read_csv(handle, nrows=0).columns)
        n_seen += 1
        for model in PREFERRED_MODELS:
            coverage[model] += int(model in columns)
    threshold = max(1, math.ceil(min_assay_coverage * n_seen))
    selected = [m for m in PREFERRED_MODELS if coverage[m] >= threshold]
    return selected, coverage


def rank_within_assay(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.rank(method="average", pct=True)


def load_human_data(
    score_zip_path: Path,
    reference_path: Path,
    min_mutants: int,
    min_assay_coverage: float,
    min_row_model_fraction: float,
) -> tuple[pd.DataFrame, list[str], dict]:
    ref = pd.read_csv(reference_path)
    ref = ref[ref["taxon"].astype(str).str.lower().eq("human")].copy()
    ref = ref[pd.to_numeric(ref["DMS_number_single_mutants"], errors="coerce") >= min_mutants].copy()
    ref = ref[~ref["includes_multiple_mutants"].map(boolish)].copy()
    records: list[pd.DataFrame] = []
    failures: list[dict] = []

    with zipfile.ZipFile(score_zip_path) as archive:
        members = archive_member_map(archive)
        models, coverage = choose_available_models(
            archive, ref, members, min_assay_coverage=min_assay_coverage
        )
        if len(models) < 6:
            raise RuntimeError(
                f"Only {len(models)} preferred models met assay coverage: {models}; coverage={coverage}"
            )
        usecols = ["mutant", "DMS_score"] + models
        for index, (_, meta) in enumerate(ref.iterrows(), start=1):
            member = assay_member(meta, members)
            if member is None:
                failures.append({"DMS_id": meta["DMS_id"], "reason": "score file missing"})
                continue
            try:
                with archive.open(member) as handle:
                    df = pd.read_csv(handle, usecols=lambda c: c in usecols, low_memory=False)
                if "mutant" not in df or "DMS_score" not in df:
                    raise ValueError("required columns missing")
                df["mutant"] = df["mutant"].astype(str)
                df = df[~df["mutant"].str.contains(":", regex=False)].copy()
                df["DMS_score"] = pd.to_numeric(df["DMS_score"], errors="coerce")
                for model in models:
                    if model not in df:
                        df[model] = np.nan
                    df[model] = pd.to_numeric(df[model], errors="coerce")
                df = df.dropna(subset=["DMS_score"]).drop_duplicates("mutant")
                minimum_models = max(3, math.ceil(min_row_model_fraction * len(models)))
                df = df[df[models].notna().sum(axis=1) >= minimum_models].copy()
                if len(df) < min_mutants:
                    raise ValueError(f"only {len(df)} rows after coverage filtering")
                df["DMS_id"] = str(meta["DMS_id"])
                df["UniProt_ID"] = str(meta["UniProt_ID"])
                df["selection_type"] = str(meta.get("coarse_selection_type", "Unknown"))
                df["msa_depth"] = str(meta.get("MSA_Neff_L_category", "Unknown"))
                df["target_rank"] = rank_within_assay(df["DMS_score"])
                for model in models:
                    df[f"rank__{model}"] = rank_within_assay(df[model])
                keep = [
                    "DMS_id", "UniProt_ID", "selection_type", "msa_depth",
                    "mutant", "DMS_score", "target_rank",
                ] + [f"rank__{m}" for m in models]
                records.append(df[keep])
                print(f"loaded {index}/{len(ref)}: {meta['DMS_id']} n={len(df)}", flush=True)
            except Exception as exc:
                failures.append({"DMS_id": str(meta["DMS_id"]), "reason": repr(exc)})

    if len(records) < 20:
        raise RuntimeError(f"Only {len(records)} usable human assays; failures={failures[:10]}")
    data = pd.concat(records, ignore_index=True)
    audit = {
        "n_reference_human_assays": int(len(ref)),
        "n_usable_assays": int(data["DMS_id"].nunique()),
        "n_uniprot": int(data["UniProt_ID"].nunique()),
        "n_mutants": int(len(data)),
        "selected_models": models,
        "model_assay_coverage": coverage,
        "failures": failures,
    }
    return data, models, audit


def feature_columns(models: list[str]) -> list[str]:
    return [f"rank__{model}" for model in models]


def assay_performance(df: pd.DataFrame, prediction: str) -> pd.DataFrame:
    rows = []
    for assay, group in df.groupby("DMS_id", sort=False):
        rows.append({
            "DMS_id": assay,
            "UniProt_ID": group["UniProt_ID"].iloc[0],
            "selection_type": group["selection_type"].iloc[0],
            "msa_depth": group["msa_depth"].iloc[0],
            "n_mutants": int(len(group)),
            "spearman": safe_spearman(group[prediction], group["DMS_score"]),
        })
    return pd.DataFrame(rows)


def official_hierarchy_score(assay_scores: pd.DataFrame) -> float:
    clean = assay_scores.dropna(subset=["spearman"]).copy()
    if clean.empty:
        return float("nan")
    protein_function = (
        clean.groupby(["UniProt_ID", "selection_type"], as_index=False)["spearman"].mean()
    )
    function_means = protein_function.groupby("selection_type")["spearman"].mean()
    return float(function_means.mean())


def model_training_scores(train: pd.DataFrame, models: list[str]) -> dict[str, float]:
    scores = {}
    for model in models:
        column = f"rank__{model}"
        scores[model] = official_hierarchy_score(assay_performance(train, column))
    return scores


def row_weighted_average(frame: pd.DataFrame, columns: list[str], weights: np.ndarray) -> np.ndarray:
    matrix = frame[columns].to_numpy(dtype=float)
    finite = np.isfinite(matrix)
    weighted = np.where(finite, matrix * weights[None, :], 0.0).sum(axis=1)
    denominator = np.where(finite, weights[None, :], 0.0).sum(axis=1)
    fallback = np.nanmedian(matrix, axis=1)
    return np.where(denominator > 0, weighted / denominator, fallback)


def balanced_training_sample(train: pd.DataFrame, max_per_assay: int = 1500) -> tuple[pd.DataFrame, np.ndarray]:
    pieces = []
    sample_weights = []
    for _, group in train.groupby("DMS_id", sort=False):
        if len(group) > max_per_assay:
            chosen = group.sample(max_per_assay, random_state=SEED)
        else:
            chosen = group
        pieces.append(chosen)
        sample_weights.extend([1.0 / len(chosen)] * len(chosen))
    sampled = pd.concat(pieces, ignore_index=True)
    weights = np.asarray(sample_weights, dtype=float)
    weights *= len(weights) / weights.sum()
    return sampled, weights


def predict_method(
    method: str,
    train: pd.DataFrame,
    test: pd.DataFrame,
    models: list[str],
) -> np.ndarray:
    columns = feature_columns(models)
    if method == "mean_all":
        return test[columns].mean(axis=1, skipna=True).to_numpy(float)
    if method == "median_all":
        return test[columns].median(axis=1, skipna=True).to_numpy(float)

    train_scores = model_training_scores(train, models)
    ordered = sorted(models, key=lambda m: (-np.nan_to_num(train_scores[m], nan=-999), m))
    if method == "top3_mean":
        return test[[f"rank__{m}" for m in ordered[:3]]].mean(axis=1, skipna=True).to_numpy(float)
    if method == "top5_mean":
        return test[[f"rank__{m}" for m in ordered[:5]]].mean(axis=1, skipna=True).to_numpy(float)
    if method == "performance_weighted":
        raw = np.asarray([max(0.0, np.nan_to_num(train_scores[m], nan=0.0)) ** 2 for m in models])
        if raw.sum() <= 0:
            raw = np.ones(len(models), dtype=float)
        return row_weighted_average(test, columns, raw / raw.sum())
    if method.startswith("ridge_"):
        alpha = float(method.split("_", 1)[1])
        sampled, weights = balanced_training_sample(train)
        pipeline = make_pipeline(
            SimpleImputer(strategy="median", add_indicator=True),
            Ridge(alpha=alpha),
        )
        pipeline.fit(sampled[columns], sampled["target_rank"], ridge__sample_weight=weights)
        return pipeline.predict(test[columns])
    raise ValueError(f"Unknown method: {method}")


def inner_select(
    outer_train: pd.DataFrame,
    models: list[str],
    candidates: list[str],
) -> tuple[str, str, dict[str, float], dict[str, float]]:
    proteins = outer_train[["UniProt_ID"]].drop_duplicates().reset_index(drop=True)
    n_splits = min(3, len(proteins))
    splitter = GroupKFold(n_splits=n_splits)
    method_predictions = {method: [] for method in candidates}
    base_predictions = {model: [] for model in models}

    for _, (train_index, validation_index) in enumerate(
        splitter.split(proteins, groups=proteins["UniProt_ID"]), start=1
    ):
        inner_train_ids = set(proteins.iloc[train_index]["UniProt_ID"])
        inner_valid_ids = set(proteins.iloc[validation_index]["UniProt_ID"])
        inner_train = outer_train[outer_train["UniProt_ID"].isin(inner_train_ids)].copy()
        inner_valid = outer_train[outer_train["UniProt_ID"].isin(inner_valid_ids)].copy()
        for method in candidates:
            part = inner_valid[["DMS_id", "UniProt_ID", "selection_type", "msa_depth", "DMS_score"]].copy()
            part["prediction"] = predict_method(method, inner_train, inner_valid, models)
            method_predictions[method].append(part)
        for model in models:
            part = inner_valid[["DMS_id", "UniProt_ID", "selection_type", "msa_depth", "DMS_score"]].copy()
            part["prediction"] = inner_valid[f"rank__{model}"].to_numpy(float)
            base_predictions[model].append(part)

    method_scores = {
        method: official_hierarchy_score(
            assay_performance(pd.concat(parts, ignore_index=True), "prediction")
        )
        for method, parts in method_predictions.items()
    }
    base_scores = {
        model: official_hierarchy_score(
            assay_performance(pd.concat(parts, ignore_index=True), "prediction")
        )
        for model, parts in base_predictions.items()
    }
    selected_method = sorted(method_scores, key=lambda x: (-np.nan_to_num(method_scores[x], nan=-999), x))[0]
    selected_base = sorted(base_scores, key=lambda x: (-np.nan_to_num(base_scores[x], nan=-999), x))[0]
    return selected_method, selected_base, method_scores, base_scores


def stratified_bootstrap_difference(cells: pd.DataFrame, n_boot: int = 20000) -> tuple[list[float], np.ndarray]:
    rng = np.random.default_rng(SEED)
    observed_types = sorted(cells["selection_type"].dropna().unique())
    values = np.empty(n_boot, dtype=float)
    for iteration in range(n_boot):
        type_means = []
        for selection_type in observed_types:
            group = cells[cells["selection_type"] == selection_type]
            sampled = group.sample(len(group), replace=True, random_state=int(rng.integers(0, 2**31 - 1)))
            type_means.append(sampled["difference"].mean())
        values[iteration] = np.mean(type_means)
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))], values


def sign_flip_pvalue(cells: pd.DataFrame, n_perm: int = 100000) -> float:
    rng = np.random.default_rng(SEED)
    observed_types = sorted(cells["selection_type"].dropna().unique())
    observed = float(np.mean([
        cells.loc[cells["selection_type"] == kind, "difference"].mean()
        for kind in observed_types
    ]))
    differences = cells["difference"].to_numpy(float)
    labels = cells["selection_type"].to_numpy(str)
    count = 0
    for _ in range(n_perm):
        flipped = differences * rng.choice([-1.0, 1.0], size=len(differences))
        permuted = float(np.mean([
            flipped[labels == kind].mean() for kind in observed_types
        ]))
        count += int(permuted >= observed)
    return float((count + 1) / (n_perm + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-zip", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-mutants", type=int, default=100)
    parser.add_argument("--min-assay-coverage", type=float, default=0.90)
    parser.add_argument("--min-row-model-fraction", type=float, default=0.75)
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data, models, audit = load_human_data(
        Path(args.score_zip), Path(args.reference), args.min_mutants,
        args.min_assay_coverage, args.min_row_model_fraction,
    )
    proteins = data[["UniProt_ID"]].drop_duplicates().reset_index(drop=True)
    outer = GroupKFold(n_splits=min(5, len(proteins)))
    data["outer_fold"] = np.nan
    data["prediction_nested_ensemble"] = np.nan
    data["prediction_nested_base"] = np.nan
    for method in CANDIDATE_METHODS:
        data[f"prediction_fixed__{method}"] = np.nan

    fold_records = []
    for fold, (train_index, test_index) in enumerate(
        outer.split(proteins, groups=proteins["UniProt_ID"]), start=1
    ):
        train_ids = set(proteins.iloc[train_index]["UniProt_ID"])
        test_ids = set(proteins.iloc[test_index]["UniProt_ID"])
        train = data[data["UniProt_ID"].isin(train_ids)].copy()
        test = data[data["UniProt_ID"].isin(test_ids)].copy()
        selected_method, selected_base, method_scores, base_scores = inner_select(
            train, models, CANDIDATE_METHODS
        )
        data.loc[test.index, "outer_fold"] = fold
        data.loc[test.index, "prediction_nested_ensemble"] = predict_method(
            selected_method, train, test, models
        )
        data.loc[test.index, "prediction_nested_base"] = test[f"rank__{selected_base}"].to_numpy(float)
        for method in CANDIDATE_METHODS:
            data.loc[test.index, f"prediction_fixed__{method}"] = predict_method(
                method, train, test, models
            )
        fold_records.append({
            "fold": fold,
            "n_train_proteins": len(train_ids),
            "n_test_proteins": len(test_ids),
            "n_train_assays": int(train["DMS_id"].nunique()),
            "n_test_assays": int(test["DMS_id"].nunique()),
            "selected_method": selected_method,
            "selected_base": selected_base,
            "inner_method_scores": method_scores,
            "inner_base_scores": base_scores,
        })
        print(
            f"outer fold {fold}: method={selected_method}, base={selected_base}, "
            f"test proteins={len(test_ids)}", flush=True
        )

    primary_ensemble = assay_performance(data, "prediction_nested_ensemble").rename(
        columns={"spearman": "ensemble_spearman"}
    )
    primary_base = assay_performance(data, "prediction_nested_base").rename(
        columns={"spearman": "base_spearman"}
    )[["DMS_id", "base_spearman"]]
    comparison = primary_ensemble.merge(primary_base, on="DMS_id", how="inner")
    comparison["difference"] = comparison["ensemble_spearman"] - comparison["base_spearman"]

    cells = (
        comparison.groupby(["UniProt_ID", "selection_type"], as_index=False)[
            ["ensemble_spearman", "base_spearman", "difference"]
        ].mean()
    )
    ensemble_score = official_hierarchy_score(
        primary_ensemble.rename(columns={"ensemble_spearman": "spearman"})
    )
    base_score = official_hierarchy_score(
        comparison.rename(columns={"base_spearman": "spearman"})[
            ["DMS_id", "UniProt_ID", "selection_type", "msa_depth", "n_mutants", "spearman"]
        ]
    )
    ci, bootstrap_values = stratified_bootstrap_difference(cells)
    p_value = sign_flip_pvalue(cells)

    exploratory = []
    for method in CANDIDATE_METHODS:
        column = f"prediction_fixed__{method}"
        scores = assay_performance(data, column)
        exploratory.append({
            "method": method,
            "official_hierarchy_spearman": official_hierarchy_score(scores),
            "mean_assay_spearman": float(scores["spearman"].mean()),
            "median_assay_spearman": float(scores["spearman"].median()),
        })
    exploratory_df = pd.DataFrame(exploratory).sort_values(
        "official_hierarchy_spearman", ascending=False
    )

    result = {
        "benchmark": "ProteinGym v1.3 human single-substitution DMS",
        **audit,
        "outer_folds": fold_records,
        "primary_nested_ensemble_official_hierarchy_spearman": ensemble_score,
        "primary_training_selected_base_official_hierarchy_spearman": base_score,
        "primary_difference": ensemble_score - base_score,
        "stratified_protein_function_bootstrap_95ci": ci,
        "one_sided_sign_flip_p": p_value,
        "fraction_assays_improved": float((comparison["difference"] > 0).mean()),
        "fraction_protein_function_cells_improved": float((cells["difference"] > 0).mean()),
        "mean_assay_difference": float(comparison["difference"].mean()),
        "median_assay_difference": float(comparison["difference"].median()),
        "decision": (
            "ROBUST_HELD_OUT_ADVANCE"
            if ci[0] > 0 and p_value < 0.05
            else "NO_ROBUST_HELD_OUT_ADVANCE"
        ),
        "claim_boundary": (
            "This tests a nested, cross-protein meta-ensemble of existing public predictors. "
            "It does not solve all human missense effects or establish clinical pathogenicity."
        ),
    }

    prediction_columns = [
        "DMS_id", "UniProt_ID", "selection_type", "msa_depth", "mutant",
        "DMS_score", "outer_fold", "prediction_nested_ensemble",
        "prediction_nested_base",
    ] + [f"prediction_fixed__{m}" for m in CANDIDATE_METHODS]
    data[prediction_columns].to_csv(
        out / "heldout_predictions.csv.gz", index=False, compression="gzip"
    )
    comparison.to_csv(out / "assay_level_comparison.csv", index=False)
    cells.to_csv(out / "protein_function_cell_comparison.csv", index=False)
    exploratory_df.to_csv(out / "fixed_method_exploratory.csv", index=False)
    pd.DataFrame(fold_records).to_json(
        out / "fold_selection.json", orient="records", indent=2
    )
    np.save(out / "bootstrap_differences.npy", bootstrap_values)
    (out / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    report = f"""# ProteinGym human missense benchmark — nested cross-protein ensemble

## Exact question

Can a leakage-controlled ensemble of existing public zero-shot predictors rank the
functional effects of substitutions in completely held-out human proteins better than
the best individual predictor selected using training proteins only?

## Data

- Human DMS assays: **{audit['n_usable_assays']}**
- Unique UniProt proteins: **{audit['n_uniprot']}**
- Single substitutions: **{audit['n_mutants']:,}**
- Predeclared predictor families retained: **{len(models)}**

## Primary result

- Nested ensemble ProteinGym-hierarchy Spearman: **{ensemble_score:.4f}**
- Training-selected individual baseline: **{base_score:.4f}**
- Difference: **{ensemble_score-base_score:+.4f}**
- Stratified protein/function bootstrap 95% CI: **[{ci[0]:+.4f}, {ci[1]:+.4f}]**
- One-sided cell sign-flip p: **{p_value:.6f}**
- Assays improved: **{100*result['fraction_assays_improved']:.1f}%**
- Protein/function cells improved: **{100*result['fraction_protein_function_cells_improved']:.1f}%**
- Prespecified decision: **`{result['decision']}`**

## Leakage controls

Outer and inner folds were grouped by UniProt ID. No assay from a held-out protein was
used to select the ensemble method, select the comparator model, estimate weights, or fit
regression coefficients. Predictor values were converted to within-assay percentile ranks;
ProteinGym's score archive convention is higher-is-better.

## Interpretation boundary

A robust positive result would be a benchmark advance for combining public predictors,
not a solution to the complete human missense-variant problem. A null result would show
that ordinary cross-protein stacking is not the missing solution.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
