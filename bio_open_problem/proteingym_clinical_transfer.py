#!/usr/bin/env python3
"""Prospective-style transfer of a DMS-trained missense ensemble to clinical variants.

The model is trained only on ProteinGym human DMS assays. Clinical labels are not
opened until an unlabeled prediction file has been written and SHA-256 hashed.
The primary test is restricted to clinical proteins absent from the DMS training
proteins.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline

SEED = 20260817
FIXED_ALPHA = 100.0
PREFERRED_MODELS = [
    "VenusREM", "ProSST-2048", "S3F_MSA", "PoET", "ESM3", "VespaG",
    "SaProt_650M_AF2", "TranceptEVE_L", "GEMME", "ProteinMPNN",
    "ESM1v_ensemble", "MSA_Transformer_ensemble", "ProtSSN_ensemble",
    "ESM2_15B",
]


def boolish(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def members_by_basename(archive: zipfile.ZipFile) -> dict[str, str]:
    return {
        Path(name).name: name
        for name in archive.namelist()
        if name.lower().endswith(".csv") and not name.endswith("/")
    }


def member_for_row(row: pd.Series, members: dict[str, str]) -> str | None:
    candidates = [
        str(row.get("DMS_filename", "")),
        f"{row.get('DMS_id', '')}.csv",
    ]
    for candidate in candidates:
        if candidate in members:
            return members[candidate]
    return None


def rank01(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").rank(method="average", pct=True)


def balanced_sample(frame: pd.DataFrame, max_per_assay: int = 1500):
    pieces = []
    weights = []
    for _, group in frame.groupby("DMS_id", sort=False):
        chosen = (
            group.sample(max_per_assay, random_state=SEED)
            if len(group) > max_per_assay else group
        )
        pieces.append(chosen)
        weights.extend([1.0 / len(chosen)] * len(chosen))
    sampled = pd.concat(pieces, ignore_index=True)
    weight_array = np.asarray(weights, dtype=float)
    weight_array *= len(weight_array) / weight_array.sum()
    return sampled, weight_array


def select_model_intersection(
    dms_zip: zipfile.ZipFile,
    dms_ref: pd.DataFrame,
    clinical_zip: zipfile.ZipFile,
    clinical_ref: pd.DataFrame,
    config: dict,
    min_fraction: float = 0.95,
):
    dms_members = members_by_basename(dms_zip)
    clinical_members = members_by_basename(clinical_zip)
    dms_rows = dms_ref.head(min(30, len(dms_ref)))
    clinical_rows = clinical_ref.head(min(100, len(clinical_ref)))
    dms_coverage = {model: 0 for model in PREFERRED_MODELS}
    clinical_coverage = {model: 0 for model in PREFERRED_MODELS}
    dms_seen = 0
    clinical_seen = 0
    for _, row in dms_rows.iterrows():
        member = member_for_row(row, dms_members)
        if member is None:
            continue
        with dms_zip.open(member) as handle:
            columns = set(pd.read_csv(handle, nrows=0).columns)
        dms_seen += 1
        for model in PREFERRED_MODELS:
            dms_coverage[model] += int(model in columns)
    for _, row in clinical_rows.iterrows():
        member = member_for_row(row, clinical_members)
        if member is None:
            continue
        with clinical_zip.open(member) as handle:
            columns = set(pd.read_csv(handle, nrows=0).columns)
        clinical_seen += 1
        for model in PREFERRED_MODELS:
            clinical_coverage[model] += int(model in columns)
    clinical_config = config["model_list_zero_shot_substitutions_clinical"]
    selected = [
        model for model in PREFERRED_MODELS
        if model in clinical_config
        and dms_coverage[model] >= math.ceil(min_fraction * max(1, dms_seen))
        and clinical_coverage[model] >= math.ceil(min_fraction * max(1, clinical_seen))
    ]
    directions = {
        model: float(clinical_config[model].get("directionality", 1))
        for model in selected
    }
    return selected, directions, {
        "dms_probe_files": dms_seen,
        "clinical_probe_files": clinical_seen,
        "dms_probe_coverage": dms_coverage,
        "clinical_probe_coverage": clinical_coverage,
    }


def load_dms_training(
    archive: zipfile.ZipFile,
    reference: pd.DataFrame,
    models: list[str],
    min_mutants: int = 100,
):
    members = members_by_basename(archive)
    records = []
    failures = []
    for index, (_, row) in enumerate(reference.iterrows(), 1):
        member = member_for_row(row, members)
        if member is None:
            failures.append({"DMS_id": row.get("DMS_id"), "reason": "missing"})
            continue
        try:
            with archive.open(member) as handle:
                data = pd.read_csv(
                    handle,
                    usecols=lambda column: column in {"mutant", "DMS_score", *models},
                    low_memory=False,
                )
            if "mutant" not in data or "DMS_score" not in data:
                raise ValueError("required DMS columns absent")
            data["mutant"] = data["mutant"].astype(str)
            data = data[~data["mutant"].str.contains(":", regex=False)].copy()
            data["DMS_score"] = pd.to_numeric(data["DMS_score"], errors="coerce")
            for model in models:
                data[model] = pd.to_numeric(data.get(model), errors="coerce")
            data = data.dropna(subset=["DMS_score"]).drop_duplicates("mutant")
            data = data[data[models].notna().sum(axis=1) >= max(3, math.ceil(0.75 * len(models)))].copy()
            if len(data) < min_mutants:
                raise ValueError(f"only {len(data)} rows")
            data["DMS_id"] = str(row["DMS_id"])
            data["UniProt_ID"] = str(row["UniProt_ID"])
            data["target_rank"] = rank01(data["DMS_score"])
            for model in models:
                data[f"rank__{model}"] = rank01(data[model])
            records.append(
                data[["DMS_id", "UniProt_ID", "mutant", "target_rank"] + [f"rank__{m}" for m in models]]
            )
            print(f"DMS {index}/{len(reference)} {row['DMS_id']} n={len(data)}", flush=True)
        except Exception as exc:
            failures.append({"DMS_id": str(row.get("DMS_id")), "reason": repr(exc)})
    return pd.concat(records, ignore_index=True), failures


def fit_dms_model(training: pd.DataFrame, models: list[str]):
    columns = [f"rank__{model}" for model in models]
    sampled, weights = balanced_sample(training)
    pipeline = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True),
        Ridge(alpha=FIXED_ALPHA),
    )
    pipeline.fit(
        sampled[columns], sampled["target_rank"], ridge__sample_weight=weights
    )
    # DMS-only rankings used for prelabel baselines.
    model_scores = {}
    for model in models:
        correlations = []
        for _, group in training.groupby("DMS_id"):
            valid = group[[f"rank__{model}", "target_rank"]].dropna()
            if len(valid) >= 20:
                correlations.append(valid.corr(method="spearman").iloc[0, 1])
        model_scores[model] = float(np.nanmean(correlations))
    ordered = sorted(models, key=lambda m: (-model_scores[m], m))
    return pipeline, columns, model_scores, ordered


def clinical_feature_predictions(
    score_archive: zipfile.ZipFile,
    reference: pd.DataFrame,
    models: list[str],
    directions: dict[str, float],
    pipeline,
    columns: list[str],
    ordered_models: list[str],
):
    members = members_by_basename(score_archive)
    records = []
    failures = []
    for index, (_, row) in enumerate(reference.iterrows(), 1):
        member = member_for_row(row, members)
        if member is None:
            failures.append({"DMS_id": row.get("DMS_id"), "reason": "score missing"})
            continue
        try:
            with score_archive.open(member) as handle:
                frame = pd.read_csv(
                    handle,
                    usecols=lambda column: column in {"mutant", *models},
                    low_memory=False,
                )
            if "mutant" not in frame:
                raise ValueError("mutant absent")
            frame["mutant"] = frame["mutant"].astype(str)
            frame = frame[~frame["mutant"].str.contains(":", regex=False)].drop_duplicates("mutant")
            for model in models:
                raw = pd.to_numeric(frame.get(model), errors="coerce") * directions[model]
                frame[f"rank__{model}"] = rank01(raw)
            frame = frame[frame[[f"rank__{m}" for m in models]].notna().sum(axis=1) >= max(3, math.ceil(0.75 * len(models)))].copy()
            if frame.empty:
                raise ValueError("no variants after coverage filter")
            frame["DMS_id"] = str(row["DMS_id"])
            frame["UniProt_ID"] = str(row.get("UniProt_ID", ""))
            frame["prediction_ridge100"] = pipeline.predict(frame[columns])
            frame["prediction_venusrem"] = frame.get("rank__VenusREM", np.nan)
            frame["prediction_top5"] = frame[[f"rank__{m}" for m in ordered_models[:5]]].mean(axis=1, skipna=True)
            frame["prediction_mean"] = frame[[f"rank__{m}" for m in models]].mean(axis=1, skipna=True)
            records.append(
                frame[["DMS_id", "UniProt_ID", "mutant", "prediction_ridge100", "prediction_venusrem", "prediction_top5", "prediction_mean"]]
            )
            if index % 100 == 0:
                print(f"clinical scores {index}/{len(reference)}", flush=True)
        except Exception as exc:
            failures.append({"DMS_id": str(row.get("DMS_id")), "reason": repr(exc)})
    return pd.concat(records, ignore_index=True), failures


def load_clinical_labels(
    label_archive: zipfile.ZipFile,
    reference: pd.DataFrame,
):
    members = members_by_basename(label_archive)
    records = []
    failures = []
    for _, row in reference.iterrows():
        member = member_for_row(row, members)
        if member is None:
            failures.append({"DMS_id": row.get("DMS_id"), "reason": "label missing"})
            continue
        try:
            with label_archive.open(member) as handle:
                frame = pd.read_csv(handle, low_memory=False)
            label_column = next(
                (column for column in ["DMS_bin_score", "DMS_score_bin"] if column in frame),
                None,
            )
            if "mutant" not in frame or label_column is None:
                raise ValueError(f"columns={list(frame.columns)[:15]}")
            frame["mutant"] = frame["mutant"].astype(str)
            raw = frame[label_column].replace({"Pathogenic": 1, "Benign": 0})
            raw = pd.to_numeric(raw, errors="coerce")
            # ProteinGym's official clinical script flips labels so 1 means benign/fit,
            # matching the higher-is-better DMS convention.
            frame["benign_label"] = 1 - raw
            frame["DMS_id"] = str(row["DMS_id"])
            frame["UniProt_ID"] = str(row.get("UniProt_ID", ""))
            records.append(
                frame[["DMS_id", "UniProt_ID", "mutant", "benign_label"]]
                .dropna(subset=["benign_label"])
                .drop_duplicates("mutant")
            )
        except Exception as exc:
            failures.append({"DMS_id": str(row.get("DMS_id")), "reason": repr(exc)})
    return pd.concat(records, ignore_index=True), failures


def gene_aucs(frame: pd.DataFrame, prediction_columns: list[str], min_class: int = 5):
    rows = []
    for (gene, uniprot), group in frame.groupby(["DMS_id", "UniProt_ID"], sort=False):
        labels = group["benign_label"].to_numpy(float)
        n_benign = int((labels == 1).sum())
        n_pathogenic = int((labels == 0).sum())
        if min(n_benign, n_pathogenic) < min_class:
            continue
        row = {
            "DMS_id": gene,
            "UniProt_ID": uniprot,
            "n_variants": len(group),
            "n_benign": n_benign,
            "n_pathogenic": n_pathogenic,
        }
        for column in prediction_columns:
            valid = group[["benign_label", column]].dropna()
            row[column] = (
                float(roc_auc_score(valid["benign_label"], valid[column]))
                if valid["benign_label"].nunique() == 2 else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def paired_bootstrap(values: np.ndarray, n_boot: int = 20000):
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(SEED)
    sampled = np.empty(n_boot)
    for index in range(n_boot):
        sampled[index] = np.mean(rng.choice(values, size=len(values), replace=True))
    return [float(np.quantile(sampled, 0.025)), float(np.quantile(sampled, 0.975))]


def sign_flip(values: np.ndarray, n_perm: int = 100000):
    values = values[np.isfinite(values)]
    observed = float(np.mean(values))
    rng = np.random.default_rng(SEED)
    count = 0
    for _ in range(n_perm):
        permuted = float(np.mean(values * rng.choice([-1.0, 1.0], size=len(values))))
        count += int(permuted >= observed)
    return float((count + 1) / (n_perm + 1))


def summarize_auc(aucs: pd.DataFrame, ensemble: str, baseline: str):
    paired = aucs[[ensemble, baseline]].dropna()
    differences = (paired[ensemble] - paired[baseline]).to_numpy(float)
    return {
        "n_genes": int(len(paired)),
        "ensemble_mean_auc": float(paired[ensemble].mean()),
        "baseline_mean_auc": float(paired[baseline].mean()),
        "mean_difference": float(np.mean(differences)),
        "median_difference": float(np.median(differences)),
        "paired_bootstrap_95ci": paired_bootstrap(differences),
        "one_sided_sign_flip_p": sign_flip(differences),
        "fraction_genes_improved": float(np.mean(differences > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dms-scores", required=True)
    parser.add_argument("--dms-reference", required=True)
    parser.add_argument("--clinical-scores", required=True)
    parser.add_argument("--clinical-labels", required=True)
    parser.add_argument("--clinical-reference", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dms_ref = pd.read_csv(args.dms_reference)
    dms_ref = dms_ref[dms_ref["taxon"].astype(str).str.lower().eq("human")].copy()
    dms_ref = dms_ref[pd.to_numeric(dms_ref["DMS_number_single_mutants"], errors="coerce") >= 100]
    dms_ref = dms_ref[~dms_ref["includes_multiple_mutants"].map(boolish)].copy()
    clinical_ref = pd.read_csv(args.clinical_reference)
    config = json.loads(Path(args.config).read_text())

    with zipfile.ZipFile(args.dms_scores) as dms_zip, zipfile.ZipFile(args.clinical_scores) as clinical_zip:
        models, directions, coverage_audit = select_model_intersection(
            dms_zip, dms_ref, clinical_zip, clinical_ref, config
        )
        if len(models) < 6:
            raise RuntimeError(f"Only {len(models)} transferable models: {models}")
        training, dms_failures = load_dms_training(dms_zip, dms_ref, models)
        pipeline, feature_cols, model_scores, ordered_models = fit_dms_model(training, models)
        predictions, score_failures = clinical_feature_predictions(
            clinical_zip, clinical_ref, models, directions, pipeline,
            feature_cols, ordered_models,
        )

    # Freeze predictions before opening clinical labels.
    prediction_path = out / "clinical_predictions_before_labels.csv.gz"
    predictions.to_csv(prediction_path, index=False, compression="gzip")
    prediction_hash = sha256(prediction_path)
    freeze = {
        "prediction_file": prediction_path.name,
        "prediction_sha256": prediction_hash,
        "fixed_alpha": FIXED_ALPHA,
        "models": models,
        "directions": directions,
        "dms_training_assays": int(training["DMS_id"].nunique()),
        "dms_training_uniprot": int(training["UniProt_ID"].nunique()),
        "dms_training_mutants": int(len(training)),
        "created_before_clinical_labels_opened": True,
    }
    (out / "PREDICTION_FREEZE.json").write_text(json.dumps(freeze, indent=2))

    with zipfile.ZipFile(args.clinical_labels) as label_zip:
        labels, label_failures = load_clinical_labels(label_zip, clinical_ref)
    evaluated = predictions.merge(
        labels[["DMS_id", "mutant", "benign_label"]],
        on=["DMS_id", "mutant"], how="inner",
    )
    dms_uniprots = set(training["UniProt_ID"].astype(str))
    evaluated["protein_unseen_in_dms"] = ~evaluated["UniProt_ID"].astype(str).isin(dms_uniprots)
    prediction_columns = [
        "prediction_ridge100", "prediction_venusrem", "prediction_top5", "prediction_mean"
    ]
    all_auc = gene_aucs(evaluated, prediction_columns)
    unseen_auc = gene_aucs(evaluated[evaluated["protein_unseen_in_dms"]], prediction_columns)

    primary = summarize_auc(unseen_auc, "prediction_ridge100", "prediction_venusrem")
    secondary_all = summarize_auc(all_auc, "prediction_ridge100", "prediction_venusrem")
    secondary_top5 = summarize_auc(unseen_auc, "prediction_top5", "prediction_venusrem")
    result = {
        "exact_question": (
            "Does a DMS-trained ridge ensemble improve clinical benign/pathogenic "
            "classification over VenusREM on proteins absent from DMS training?"
        ),
        "prediction_freeze": freeze,
        "coverage_audit": coverage_audit,
        "dms_failures": dms_failures,
        "clinical_score_failures_n": len(score_failures),
        "clinical_label_failures_n": len(label_failures),
        "clinical_prediction_rows": int(len(predictions)),
        "clinical_evaluated_rows": int(len(evaluated)),
        "all_evaluable_genes": int(len(all_auc)),
        "unseen_protein_evaluable_genes": int(len(unseen_auc)),
        "primary_unseen_proteins_ridge_vs_venusrem": primary,
        "secondary_all_proteins_ridge_vs_venusrem": secondary_all,
        "secondary_unseen_top5_vs_venusrem": secondary_top5,
        "decision": (
            "ROBUST_CLINICAL_TRANSFER_ADVANCE"
            if primary["paired_bootstrap_95ci"][0] > 0
            and primary["one_sided_sign_flip_p"] < 0.05
            else "NO_ROBUST_CLINICAL_TRANSFER_ADVANCE"
        ),
        "interpretation_boundary": (
            "Clinical labels are ClinVar-style benign/pathogenic annotations, not new wet-lab "
            "measurements. This test evaluates medical transfer and does not prove causal "
            "mechanism or solve all variants of uncertain significance."
        ),
    }
    all_auc.to_csv(out / "clinical_auc_all_genes.csv", index=False)
    unseen_auc.to_csv(out / "clinical_auc_unseen_dms_proteins.csv", index=False)
    evaluated[["DMS_id", "UniProt_ID", "mutant", "benign_label", "protein_unseen_in_dms", *prediction_columns]].to_csv(
        out / "clinical_evaluated_variants.csv.gz", index=False, compression="gzip"
    )
    (out / "result.json").write_text(json.dumps(result, indent=2))
    report = f"""# DMS-trained ensemble → clinical missense transfer

## Frozen question

Does a ridge ensemble trained only on human DMS assays improve benign/pathogenic
classification over the DMS-training-selected VenusREM baseline on clinical proteins
that were completely absent from DMS training?

## Prediction freeze

- Unlabeled prediction SHA-256: `{prediction_hash}`
- Clinical labels were opened only after this file was written.
- Ridge alpha: **{FIXED_ALPHA:g}**
- Transferable public predictors: **{len(models)}**

## Primary unseen-protein result

- Evaluable genes: **{primary['n_genes']}**
- DMS-trained ridge mean gene AUC: **{primary['ensemble_mean_auc']:.4f}**
- VenusREM mean gene AUC: **{primary['baseline_mean_auc']:.4f}**
- Difference: **{primary['mean_difference']:+.4f}**
- Paired bootstrap 95% CI: **[{primary['paired_bootstrap_95ci'][0]:+.4f}, {primary['paired_bootstrap_95ci'][1]:+.4f}]**
- One-sided sign-flip p: **{primary['one_sided_sign_flip_p']:.6f}**
- Genes improved: **{100*primary['fraction_genes_improved']:.1f}%**
- Decision: **`{result['decision']}`**

## Boundary

This is a transfer test from DMS fitness to public clinical annotations. It is not an
experimental confirmation of individual variants and does not solve all VUS.
"""
    (out / "REPORT.md").write_text(report)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
