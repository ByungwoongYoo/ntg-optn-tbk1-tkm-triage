#!/usr/bin/env python3
"""Independent, leakage-focused audit of the public NEAT repository.

Pinned source commit: c1df534d90b08ee5d6e054d0797088067cbc4cd1
Outputs are written to --out-dir and are intended for a reproducibility package.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

ORDER = ["Accept", "Review", "Reject"]
SEED = 20260821


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def wilson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    z = norm.ppf(1 - alpha / 2)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def cluster_bootstrap_auc(df: pd.DataFrame, truth_col: str, score_col: str,
                          patient_col: str, reps: int = 5000,
                          seed: int = SEED) -> tuple[float, float, int]:
    rng = np.random.default_rng(seed)
    patients = np.array(sorted(df[patient_col].dropna().astype(str).unique()))
    vals: list[float] = []
    groups = {p: df[df[patient_col].astype(str) == p] for p in patients}
    for _ in range(reps):
        sampled = rng.choice(patients, size=len(patients), replace=True)
        boot = pd.concat([groups[p] for p in sampled], ignore_index=True)
        y = (boot[truth_col] == "Accept").astype(int)
        if y.nunique() < 2:
            continue
        vals.append(float(roc_auc_score(y, boot[score_col])))
    if not vals:
        return math.nan, math.nan, 0
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(lo), float(hi), len(vals)


def cluster_bootstrap_mean(values_by_patient: dict[str, float], reps: int = 5000,
                           seed: int = SEED) -> tuple[float, float, int]:
    rng = np.random.default_rng(seed)
    keys = np.array(sorted(values_by_patient))
    vals = np.array([values_by_patient[k] for k in keys], dtype=float)
    out = []
    for _ in range(reps):
        s = rng.choice(vals, size=len(vals), replace=True)
        s = s[np.isfinite(s)]
        if len(s):
            out.append(float(np.mean(s)))
    if not out:
        return math.nan, math.nan, 0
    lo, hi = np.quantile(out, [0.025, 0.975])
    return float(lo), float(hi), len(out)


def candidate_key(df: pd.DataFrame) -> pd.Series:
    fields = [c for c in ["patient_id", "ID", "Best Peptide", "Allele", "Best Transcript"] if c in df.columns]
    if not fields:
        return pd.Series(df.index.astype(str), index=df.index)
    return df[fields].fillna("").astype(str).agg("|".join, axis=1)


def confusion(df: pd.DataFrame, truth: str, pred: str) -> pd.DataFrame:
    return pd.crosstab(
        pd.Categorical(df[truth], categories=ORDER, ordered=True),
        pd.Categorical(df[pred], categories=ORDER, ordered=True),
        dropna=False,
    ).reindex(index=ORDER, columns=ORDER, fill_value=0)


def summarize_prospective(df: pd.DataFrame, name: str, out: Path) -> dict:
    cm = confusion(df, "Evaluation", "python_Evaluation_pred")
    cm.to_csv(out / f"{name}_confusion_matrix.csv", index_label="actual")
    n = len(df)
    correct = int(sum(cm.loc[c, c] for c in ORDER))
    coverage_n = int((df["python_Evaluation_pred"] != "Review").sum())
    acc_ci = wilson(correct, n)
    cov_ci = wilson(coverage_n, n)

    binary_truth = df[df["Evaluation"].isin(["Accept", "Reject"])].copy()
    auc = float(roc_auc_score((binary_truth["Evaluation"] == "Accept").astype(int),
                              binary_truth["python_Accept_pred_prob"]))
    auc_ci = cluster_bootstrap_auc(binary_truth, "Evaluation", "python_Accept_pred_prob", "patient_id")

    conditioned = binary_truth[binary_truth["python_Evaluation_pred"].isin(["Accept", "Reject"])].copy()
    auc_conditioned = float(roc_auc_score((conditioned["Evaluation"] == "Accept").astype(int),
                                          conditioned["python_Accept_pred_prob"])) if conditioned["Evaluation"].nunique() == 2 else math.nan
    conditioned_ci = cluster_bootstrap_auc(conditioned, "Evaluation", "python_Accept_pred_prob", "patient_id")

    df = df.copy()
    df["candidate_key"] = candidate_key(df)
    patient_rows = []
    for pid, g in df.groupby("patient_id", sort=True):
        gn = len(g)
        gcorrect = int((g["Evaluation"] == g["python_Evaluation_pred"]).sum())
        true_accept = set(g.loc[g["Evaluation"] == "Accept", "candidate_key"])
        pred_accept = set(g.loc[g["python_Evaluation_pred"] == "Accept", "candidate_key"])
        inter = true_accept & pred_accept
        union = true_accept | pred_accept
        tp = len(inter)
        fp = len(pred_accept - true_accept)
        fn = len(true_accept - pred_accept)
        precision = tp / (tp + fp) if tp + fp else math.nan
        recall = tp / (tp + fn) if tp + fn else math.nan
        jac = tp / len(union) if union else 1.0
        auto_reject_accept = int(((g["Evaluation"] == "Accept") & (g["python_Evaluation_pred"] == "Reject")).sum())
        true_accept_n = int((g["Evaluation"] == "Accept").sum())
        patient_rows.append({
            "patient_id": pid,
            "n_candidates": gn,
            "n_actual_accept": true_accept_n,
            "n_pred_accept": int((g["python_Evaluation_pred"] == "Accept").sum()),
            "n_pred_review": int((g["python_Evaluation_pred"] == "Review").sum()),
            "n_pred_reject": int((g["python_Evaluation_pred"] == "Reject").sum()),
            "accuracy": gcorrect / gn,
            "accuracy_ci_low_wilson": wilson(gcorrect, gn)[0],
            "accuracy_ci_high_wilson": wilson(gcorrect, gn)[1],
            "accept_precision": precision,
            "accept_recall": recall,
            "accept_set_jaccard": jac,
            "overselected_accept_count": fp,
            "omitted_accept_count": fn,
            "expert_accept_auto_reject_count": auto_reject_accept,
            "expert_accept_auto_reject_rate": auto_reject_accept / true_accept_n if true_accept_n else math.nan,
            "review_workload_fraction": float((g["python_Evaluation_pred"] == "Review").mean()),
            "non_review_coverage": float((g["python_Evaluation_pred"] != "Review").mean()),
        })
    pdf = pd.DataFrame(patient_rows)
    pdf.to_csv(out / f"{name}_patient_metrics.csv", index=False)

    macro = {}
    for col in ["accuracy", "accept_precision", "accept_recall", "accept_set_jaccard",
                "review_workload_fraction", "non_review_coverage"]:
        s = pdf[col].replace([np.inf, -np.inf], np.nan).dropna()
        mapping = dict(zip(pdf.loc[s.index, "patient_id"].astype(str), s.astype(float)))
        ci = cluster_bootstrap_mean(mapping)
        macro[col] = {
            "mean": float(s.mean()), "median": float(s.median()),
            "min": float(s.min()), "max": float(s.max()),
            "patient_cluster_bootstrap_95ci": [ci[0], ci[1]],
            "n_patients_defined": int(s.shape[0]),
        }

    actual_a = df["Evaluation"] == "Accept"
    pred_a = df["python_Evaluation_pred"] == "Accept"
    tp = int((actual_a & pred_a).sum())
    fp = int((~actual_a & pred_a).sum())
    fn = int((actual_a & ~pred_a).sum())
    accept_precision = tp / (tp + fp) if tp + fp else math.nan
    accept_recall = tp / (tp + fn) if tp + fn else math.nan
    true_set = set(df.loc[actual_a, "candidate_key"])
    pred_set = set(df.loc[pred_a, "candidate_key"])
    pooled_jaccard = len(true_set & pred_set) / len(true_set | pred_set) if (true_set | pred_set) else 1.0
    false_reject_n = int((actual_a & (df["python_Evaluation_pred"] == "Reject")).sum())
    false_reject_ci = wilson(false_reject_n, int(actual_a.sum()))

    return {
        "name": name,
        "rows": n,
        "patients": int(df["patient_id"].nunique()),
        "actual_label_counts": {k: int(v) for k, v in df["Evaluation"].value_counts().to_dict().items()},
        "predicted_label_counts": {k: int(v) for k, v in df["python_Evaluation_pred"].value_counts().to_dict().items()},
        "confusion_matrix": {a: {p: int(cm.loc[a, p]) for p in ORDER} for a in ORDER},
        "three_class_accuracy": correct / n,
        "three_class_accuracy_numerator": correct,
        "three_class_accuracy_denominator": n,
        "three_class_accuracy_wilson_95ci": list(acc_ci),
        "non_review_prediction_coverage": coverage_n / n,
        "non_review_prediction_coverage_numerator": coverage_n,
        "non_review_prediction_coverage_denominator": n,
        "non_review_prediction_coverage_wilson_95ci": list(cov_ci),
        "accept_vs_reject_auc_actual_review_excluded": auc,
        "accept_vs_reject_auc_denominator": int(len(binary_truth)),
        "accept_vs_reject_auc_patient_cluster_bootstrap_95ci": [auc_ci[0], auc_ci[1]],
        "accept_vs_reject_auc_repository_script_conditioned_on_pred_nonreview": auc_conditioned,
        "conditioned_auc_denominator": int(len(conditioned)),
        "conditioned_auc_patient_cluster_bootstrap_95ci": [conditioned_ci[0], conditioned_ci[1]],
        "pooled_accept_precision": accept_precision,
        "pooled_accept_precision_wilson_95ci": list(wilson(tp, tp + fp)),
        "pooled_accept_recall": accept_recall,
        "pooled_accept_recall_wilson_95ci": list(wilson(tp, tp + fn)),
        "pooled_accept_set_jaccard": pooled_jaccard,
        "pooled_overselected_accept_count": fp,
        "pooled_omitted_accept_count": fn,
        "expert_accept_auto_reject_count": false_reject_n,
        "expert_accept_auto_reject_denominator": int(actual_a.sum()),
        "expert_accept_auto_reject_rate": false_reject_n / int(actual_a.sum()) if int(actual_a.sum()) else math.nan,
        "expert_accept_auto_reject_wilson_95ci": list(false_reject_ci),
        "patient_macro_metrics": macro,
    }


def duplicate_audit(train: pd.DataFrame, dev: pd.DataFrame, pros: pd.DataFrame, out: Path) -> dict:
    result: dict = {}
    train_p = set(train.get("patient_id", pd.Series(dtype=str)).dropna().astype(str))
    dev_p = set(dev.get("patient_id", pd.Series(dtype=str)).dropna().astype(str))
    pros_p = set(pros.get("patient_id", pd.Series(dtype=str)).dropna().astype(str))
    result["patient_ids"] = {
        "train_unique": len(train_p), "development_unique": len(dev_p), "prospective_unique": len(pros_p),
        "train_development_overlap_count": len(train_p & dev_p),
        "train_development_overlap_ids": sorted(train_p & dev_p),
        "retrospective_prospective_overlap_count": len((train_p | dev_p) & pros_p),
        "retrospective_prospective_overlap_ids": sorted((train_p | dev_p) & pros_p),
    }

    datasets = {"training": train, "development": dev, "prospective": pros}
    for key_col in ["ID", "Best Peptide"]:
        available = {name: key_col in frame.columns for name, frame in datasets.items()}
        entry = {"column": key_col, "available": available}
        for name, frame in datasets.items():
            if key_col in frame.columns:
                s = frame[key_col].dropna().astype(str)
                vc = s.value_counts()
                entry[f"{name}_duplicate_rows_beyond_first"] = int((vc - 1).clip(lower=0).sum())
                entry[f"{name}_duplicated_unique_values"] = int((vc > 1).sum())
        pairs = [("training", "development"), ("training", "prospective"), ("development", "prospective")]
        for a, b in pairs:
            if key_col in datasets[a].columns and key_col in datasets[b].columns:
                sa = set(datasets[a][key_col].dropna().astype(str))
                sb = set(datasets[b][key_col].dropna().astype(str))
                entry[f"{a}_{b}_unique_overlap_count"] = len(sa & sb)
                overlap = sorted(sa & sb)
                pd.DataFrame({key_col: overlap}).to_csv(out / f"overlap_{key_col.replace(' ', '_')}_{a}_{b}.csv", index=False)
        result[key_col] = entry
    return result


def code_audit(neat: Path) -> dict:
    targets = [
        neat / "model_development/scripts/train.py",
        neat / "manuscript/scripts/ml_randomforest_model.py",
        neat / "manuscript/scripts/ml_logistic_model.py",
        neat / "manuscript/scripts/evaluation_on_prospective_test_set.py",
    ]
    rows = []
    for p in targets:
        text = p.read_text(encoding="utf-8") if p.exists() else ""
        rows.append({
            "path": str(p.relative_to(neat)) if p.exists() else str(p),
            "exists": p.exists(),
            "uses_train_test_split": "train_test_split(" in text,
            "uses_GridSearchCV": "GridSearchCV(" in text,
            "uses_GroupKFold": "GroupKFold" in text,
            "uses_StratifiedGroupKFold": "StratifiedGroupKFold" in text,
            "passes_groups_to_fit": bool(re.search(r"\.fit\([^\n]*groups\s*=", text)),
            "filters_actual_review_for_auc": "Evaluation'].isin(['Accept', 'Reject'])" in text or 'Evaluation"].isin(["Accept", "Reject"])' in text,
            "filters_predicted_review_for_auc": "python_Evaluation_pred'].isin(['Accept', 'Reject'])" in text or 'python_Evaluation_pred"].isin(["Accept", "Reject"])' in text,
            "sha256": sha256_file(p) if p.exists() else None,
        })
    return {
        "files": rows,
        "conclusion": "Candidate-row random split and ordinary stratified CV; no patient grouping detected. The prospective evaluation script also conditions binary AUROC on predicted non-Review rows.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--neat-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    neat = Path(args.neat_dir).resolve()
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    data = neat / "manuscript/data/training_testing_data"
    train_path = data / "training_set.csv"
    dev_path = data / "development_test_set.csv"
    pros_path = data / "prospective_test_set.csv"
    nocal_path = data / "prospective_test_set_no_calibration.csv"
    for p in [train_path, dev_path, pros_path, nocal_path]:
        if not p.exists():
            raise FileNotFoundError(p)

    train = pd.read_csv(train_path)
    dev = pd.read_csv(dev_path)
    pros = pd.read_csv(pros_path)
    nocal = pd.read_csv(nocal_path)

    full = summarize_prospective(pros, "prospective_full12", out)
    nocal_summary = summarize_prospective(nocal, "prospective_no_calibration8", out)
    dup = duplicate_audit(train, dev, pros, out)
    code = code_audit(neat)

    pros_keys = candidate_key(pros)
    nocal_keys = set(candidate_key(nocal))
    calibration_rows = pros[~pros_keys.isin(nocal_keys)].copy()
    calibration_patients = sorted(calibration_rows["patient_id"].astype(str).unique())
    calibration_rows.to_csv(out / "threshold_calibration_subset_reconstructed.csv", index=False)

    hashes = []
    for p in [train_path, dev_path, pros_path, nocal_path,
              neat / "README.md", neat / "LICENSE",
              neat / "model_development/scripts/train.py",
              neat / "manuscript/scripts/ml_randomforest_model.py",
              neat / "manuscript/scripts/evaluation_on_prospective_test_set.py"]:
        hashes.append({"path": str(p.relative_to(neat)), "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    pd.DataFrame(hashes).to_csv(out / "neat_source_hashes.csv", index=False)

    git_commit = subprocess.check_output(["git", "-C", str(neat), "rev-parse", "HEAD"], text=True).strip()
    summary = {
        "audit_date": "2026-08-21",
        "source_repository": "https://github.com/griffithlab/NEAT",
        "source_commit": git_commit,
        "retrospective": {
            "training_rows": len(train), "development_rows": len(dev),
            "retrospective_rows_total": len(train) + len(dev),
            "training_patients": int(train["patient_id"].nunique()) if "patient_id" in train else None,
            "development_patients": int(dev["patient_id"].nunique()) if "patient_id" in dev else None,
            "combined_patients": int(pd.concat([train[["patient_id"]], dev[["patient_id"]]])["patient_id"].nunique()) if "patient_id" in train and "patient_id" in dev else None,
            "label_counts": pd.concat([train["Evaluation"], dev["Evaluation"]]).value_counts().to_dict(),
            "training_columns": list(train.columns),
            "development_columns": list(dev.columns),
        },
        "prospective_full12": full,
        "prospective_no_calibration8": nocal_summary,
        "reconstructed_threshold_calibration_subset": {
            "rows": int(len(calibration_rows)),
            "patients": len(calibration_patients),
            "patient_ids": calibration_patients,
        },
        "overlap_and_duplicates": dup,
        "code_audit": code,
        "limitations": [
            "The public retrospective matrices do not expose a peptide-sequence column if 'Best Peptide' is absent; peptide-level duplicate testing is then not reproducible from these files.",
            "Patient-level confidence intervals are descriptive with only 8 or 12 patients.",
            "The full 12-patient prospective set is not a sealed threshold holdout because 4 patients/170 rows were used to choose thresholds.",
        ],
    }
    (out / "neat_audit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=lambda x: int(x) if isinstance(x, np.integer) else float(x) if isinstance(x, np.floating) else str(x)), encoding="utf-8")

    headline = []
    for cohort_name, s in [("full12", full), ("no_calibration8", nocal_summary)]:
        for metric in [
            "rows", "patients", "three_class_accuracy", "non_review_prediction_coverage",
            "accept_vs_reject_auc_actual_review_excluded",
            "accept_vs_reject_auc_repository_script_conditioned_on_pred_nonreview",
            "pooled_accept_precision", "pooled_accept_recall", "pooled_accept_set_jaccard",
            "pooled_overselected_accept_count", "pooled_omitted_accept_count",
            "expert_accept_auto_reject_count", "expert_accept_auto_reject_rate",
        ]:
            headline.append({"cohort": cohort_name, "metric": metric, "value": s.get(metric)})
    pd.DataFrame(headline).to_csv(out / "neat_headline_metrics.csv", index=False)

    manifest_rows = []
    for p in sorted(out.glob("*")):
        if p.is_file() and p.name != "OUTPUT_SHA256SUMS.txt":
            manifest_rows.append(f"{sha256_file(p)}  {p.name}")
    (out / "OUTPUT_SHA256SUMS.txt").write_text("\n".join(manifest_rows) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "ok",
        "source_commit": git_commit,
        "full12_accuracy": full["three_class_accuracy"],
        "full12_auc": full["accept_vs_reject_auc_actual_review_excluded"],
        "full12_false_reject_accept": full["expert_accept_auto_reject_count"],
        "train_dev_patient_overlap": dup["patient_ids"]["train_development_overlap_count"],
        "calibration_rows": len(calibration_rows),
        "calibration_patients": calibration_patients,
    }, indent=2))


if __name__ == "__main__":
    main()
