#!/usr/bin/env python3
"""Compatibility runner for ClinVar submission files and family-artifact column names."""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

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
        path,
        sep="\t",
        compression="gzip",
        dtype=str,
        chunksize=chunksize,
        low_memory=False,
        skiprows=skiprows,
        header=0,
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
        submission_n=("VariationID", "size"),
        any_p=("_p", "max"),
        any_b=("_b", "max"),
        any_vus=("_v", "max"),
        any_conflict_text=("_conflict", "max"),
        max_submission_review_tier=("_tier", "max"),
    ).reset_index()
    grouped["resolved_submission_conflict"] = grouped["any_p"] & grouped["any_b"]
    grouped["pure_vus_submissions"] = grouped["any_vus"] & ~grouped["any_p"] & ~grouped["any_b"] & ~grouped["any_conflict_text"]
    return grouped


_original_load_raw_and_scores = tm.load_raw_and_scores


def load_raw_and_scores_compatible(raw_zip: Path, score_zip: Path, predictions_path: Path):
    """Translate explicit column names from the authoritative family-held-out artifact.

    That artifact stores the training-selected comparator as `best_individual_prob`.
    All five outer folds selected PoET, so this is the frozen PoET comparator requested
    by the protocol. The shorthand aliases are added without changing predictions.
    """
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


tm.read_submission_summary = read_submission_summary_compatible
tm.load_raw_and_scores = load_raw_and_scores_compatible

if __name__ == "__main__":
    tm.main()
