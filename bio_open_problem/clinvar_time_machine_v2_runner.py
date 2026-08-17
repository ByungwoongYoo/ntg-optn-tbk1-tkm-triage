#!/usr/bin/env python3
"""Compatibility and performance runner for the ClinVar time-machine analysis."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, log_loss

import clinvar_time_machine_v2 as tm


def find_header_line(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle):
            stripped = line.lstrip("\ufeff").lstrip("#")
            if stripped.startswith("VariationID\t"):
                return line_no
            if line_no > 200:
                break
    raise RuntimeError(f"Could not locate VariationID header in {path}")


def read_submission_summary_compatible(path: Path, target_ids: set[int], chunksize: int = 400_000) -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    use = ["VariationID", "ClinicalSignificance", "DateLastEvaluated", "ReviewStatus", "Submitter", "SCV", "CollectionMethod"]
    skiprows = find_header_line(path)
    for chunk in pd.read_csv(
        path, sep="\t", compression="gzip", dtype=str, chunksize=chunksize,
        low_memory=False, skiprows=skiprows, header=0,
    ):
        chunk = tm.normalize_columns(chunk)
        if "VariationID" not in chunk.columns:
            raise RuntimeError(f"VariationID missing from {path}; columns={chunk.columns.tolist()[:30]}")
        ids = pd.to_numeric(chunk["VariationID"], errors="coerce")
        mask = ids.isin(target_ids)
        if not mask.any():
            continue
        sub = chunk.loc[mask, [c for c in use if c in chunk.columns]].copy()
        sub["VariationID"] = pd.to_numeric(sub["VariationID"], errors="coerce").astype("Int64")
        records.append(sub)
    if not records:
        return pd.DataFrame(columns=["VariationID", "submission_n", "any_p", "any_b", "any_vus", "any_conflict_text", "max_submission_review_tier"])
    x = pd.concat(records, ignore_index=True)
    x = x[x["VariationID"].notna()].copy()
    x["VariationID"] = x["VariationID"].astype(int)
    flags = x["ClinicalSignificance"].map(tm.significance_flags)
    x["_p"] = flags.map(lambda t: t[0])
    x["_b"] = flags.map(lambda t: t[1])
    x["_v"] = flags.map(lambda t: t[2])
    x["_conflict"] = flags.map(lambda t: t[3])
    review = x["ReviewStatus"] if "ReviewStatus" in x.columns else pd.Series("", index=x.index)
    x["_tier"] = review.map(tm.review_tier)
    grouped = x.groupby("VariationID", sort=False).agg(
        submission_n=("VariationID", "size"), any_p=("_p", "max"),
        any_b=("_b", "max"), any_vus=("_v", "max"),
        any_conflict_text=("_conflict", "max"),
        max_submission_review_tier=("_tier", "max"),
    ).reset_index()
    grouped["resolved_submission_conflict"] = grouped["any_p"] & grouped["any_b"]
    grouped["pure_vus_submissions"] = grouped["any_vus"] & ~grouped["any_p"] & ~grouped["any_b"] & ~grouped["any_conflict_text"]
    return grouped


_original_load_raw_and_scores = tm.load_raw_and_scores


def load_raw_and_scores_compatible(raw_zip: Path, score_zip: Path, predictions_path: Path):
    """Translate explicit names from the authoritative family-held-out artifact."""
    merged, diagnostics = _original_load_raw_and_scores(raw_zip, score_zip, predictions_path)
    if "base_prob" not in merged.columns and "best_individual_prob" in merged.columns:
        merged["base_prob"] = merged["best_individual_prob"]
    if "base" not in merged.columns and "best_individual_raw" in merged.columns:
        merged["base"] = merged["best_individual_raw"]
    missing = {"fixed_prob", "base_prob"} - set(merged.columns)
    if missing:
        raise RuntimeError(f"Frozen prediction artifact missing required columns after compatibility mapping: {sorted(missing)}")
    diagnostics["frozen_comparator_column_source"] = "best_individual_prob"
    diagnostics["frozen_best_individual_identity"] = "PoET in every completed family-held-out outer fold"
    return merged, diagnostics


def _auc_precompute(y: np.ndarray, p: np.ndarray, cluster_index: np.ndarray):
    order = np.argsort(p, kind="mergesort")
    ps = p[order]
    ys = y[order].astype(float)
    cs = cluster_index[order]
    starts = np.r_[0, np.flatnonzero(np.diff(ps) != 0) + 1]
    return ys, cs, starts


def _weighted_auc_from_cluster_counts(counts: np.ndarray, info) -> float:
    y_sorted, cluster_sorted, starts = info
    w = counts[cluster_sorted].astype(float)
    pos = np.add.reduceat(w * y_sorted, starts)
    neg = np.add.reduceat(w * (1.0 - y_sorted), starts)
    total_pos = pos.sum(); total_neg = neg.sum()
    if total_pos <= 0 or total_neg <= 0:
        return float("nan")
    neg_before = np.cumsum(neg) - neg
    numerator = np.sum(pos * (neg_before + 0.5 * neg))
    return float(numerator / (total_pos * total_neg))


def fast_cluster_bootstrap_compare(df: pd.DataFrame, ensemble_col: str, base_col: str, n_boot: int = 5000):
    """Same cluster-resampling estimand as the original implementation, without DataFrame concatenation."""
    work = df[["cluster_id", "current_label", ensemble_col, base_col]].dropna().copy()
    cluster_codes, cluster_names = pd.factorize(work["cluster_id"], sort=True)
    n_clusters = len(cluster_names)
    if n_clusters < 3:
        return {"n_clusters": int(n_clusters), "error": "too few clusters"}
    y = work["current_label"].to_numpy(int)
    pe = np.clip(work[ensemble_col].to_numpy(float), 1e-6, 1 - 1e-6)
    pb = np.clip(work[base_col].to_numpy(float), 1e-6, 1 - 1e-6)
    observed = {
        "brier_improvement": float(brier_score_loss(y, pb) - brier_score_loss(y, pe)),
        "logloss_improvement": float(log_loss(y, pb, labels=[0, 1]) - log_loss(y, pe, labels=[0, 1])),
        "auroc_improvement": float(tm.safe_auc(y, pe) - tm.safe_auc(y, pb)),
    }

    n_by_cluster = np.bincount(cluster_codes, minlength=n_clusters).astype(float)
    brier_row = (y - pb) ** 2 - (y - pe) ** 2
    log_row = -(y * np.log(pb) + (1-y) * np.log(1-pb)) + (y * np.log(pe) + (1-y) * np.log(1-pe))
    brier_sum = np.bincount(cluster_codes, weights=brier_row, minlength=n_clusters)
    log_sum = np.bincount(cluster_codes, weights=log_row, minlength=n_clusters)
    auc_e_info = _auc_precompute(y, pe, cluster_codes)
    auc_b_info = _auc_precompute(y, pb, cluster_codes)

    rng = np.random.default_rng(tm.SEED)
    brier_vals = np.empty(n_boot); log_vals = np.empty(n_boot); auc_vals = np.full(n_boot, np.nan)
    probs = np.full(n_clusters, 1.0 / n_clusters)
    batch_size = 250
    cursor = 0
    while cursor < n_boot:
        b = min(batch_size, n_boot - cursor)
        counts_batch = rng.multinomial(n_clusters, probs, size=b)
        denom = counts_batch @ n_by_cluster
        brier_vals[cursor:cursor+b] = (counts_batch @ brier_sum) / denom
        log_vals[cursor:cursor+b] = (counts_batch @ log_sum) / denom
        for j, counts in enumerate(counts_batch):
            ae = _weighted_auc_from_cluster_counts(counts, auc_e_info)
            ab = _weighted_auc_from_cluster_counts(counts, auc_b_info)
            auc_vals[cursor+j] = ae - ab if np.isfinite(ae) and np.isfinite(ab) else np.nan
        cursor += b

    def ci(arr):
        arr = np.asarray(arr, float); arr = arr[np.isfinite(arr)]
        return [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))] if len(arr) >= 100 else [None, None]

    return {
        "n_clusters": int(n_clusters), "observed": observed,
        "cluster_bootstrap_95ci": {
            "brier_improvement": ci(brier_vals),
            "logloss_improvement": ci(log_vals),
            "auroc_improvement": ci(auc_vals),
        },
        "bootstrap_engine": "vectorized multinomial cluster resampling; weighted tie-aware AUROC",
        "n_bootstrap": int(n_boot),
    }


tm.read_submission_summary = read_submission_summary_compatible
tm.load_raw_and_scores = load_raw_and_scores_compatible
tm.cluster_bootstrap_compare = fast_cluster_bootstrap_compare

if __name__ == "__main__":
    tm.main()
