#!/usr/bin/env python3
"""ClinVar historical VUS time-machine validation for a frozen ProteinGym ensemble.

The script evaluates two predeclared predictors on variants that were pure VUS at a
historical ClinVar cutoff and were subsequently resolved by the 2026-08 release:

1. The already-frozen family-held-out out-of-fold fixed ensemble and its PoET comparator.
2. A temporal-clean weighted5_p1 recipe whose orientations, weights, and Platt
   calibration are fit only on variants already resolved at the historical cutoff.

All score comparisons use the identical complete-case variants. ClinVar submission
summaries are used to enforce pure historical VUS status and to remove current
pathogenic-versus-benign submission conflicts. Uncertainty is resampled by the
pre-existing MMseqs global30 homology component.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score

SEED = 20260817
MODELS = ["PoET", "TranceptEVE_L", "GEMME", "EVE", "ESM1b"]
CUTOFFS = ["2021-08", "2022-08", "2023-08"]
TRUTH_TIERS = ["primary_high_confidence", "secondary_one_star_or_higher", "secondary_all_resolved"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).lstrip("#").strip() for c in df.columns]
    return df


def significance_flags(value: Any) -> tuple[bool, bool, bool, bool]:
    text = str(value or "").strip().lower().replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    explicit_conflict = "conflicting" in text
    has_p = bool(re.search(r"(?:^|[;/|,])\s*(?:likely\s+)?pathogenic(?:\s*,?\s*low\s+penetrance)?(?:$|[;/|,])", text))
    has_b = bool(re.search(r"(?:^|[;/|,])\s*(?:likely\s+)?benign(?:$|[;/|,])", text))
    has_p = has_p or text in {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic", "likely pathogenic/pathogenic"}
    has_b = has_b or text in {"benign", "likely benign", "benign/likely benign", "likely benign/benign"}
    has_v = "uncertain significance" in text or re.search(r"\bvus(?:\b|[- ])", text) is not None
    return has_p, has_b, has_v, explicit_conflict


def classification_side(value: Any) -> str:
    p, b, v, conflict = significance_flags(value)
    if conflict or (p and b):
        return "CONFLICT"
    if p:
        return "P"
    if b:
        return "B"
    if v:
        return "VUS"
    return "OTHER"


def review_tier(value: Any) -> int:
    text = str(value or "").strip().lower().replace("_", " ")
    if "practice guideline" in text:
        return 4
    if "reviewed by expert panel" in text:
        return 3
    if "criteria provided, multiple submitters, no conflicts" in text:
        return 2
    if "criteria provided, single submitter" in text:
        return 1
    return 0


def choose_variant_row(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    x = frame.copy()
    assembly = x.get("Assembly", pd.Series("", index=x.index)).astype(str)
    x["_assembly_rank"] = np.where(assembly.eq("GRCh38"), 0, np.where(assembly.eq("GRCh37"), 1, 2))
    review = x["ReviewStatus"] if "ReviewStatus" in x.columns else pd.Series("", index=x.index)
    x["_review_tier"] = review.map(review_tier)
    x = x.sort_values(["VariationID", "_assembly_rank", "_review_tier"], ascending=[True, True, False])
    return x.drop_duplicates("VariationID", keep="first").drop(columns=["_assembly_rank"], errors="ignore")


def read_variant_summary(path: Path, target_ids: set[int], chunksize: int = 250_000) -> pd.DataFrame:
    wanted = [
        "AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
        "LastEvaluated", "OriginSimple", "Assembly", "Chromosome", "Start", "Stop",
        "ReferenceAllele", "AlternateAllele", "ReviewStatus", "NumberSubmitters",
        "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
    ]
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", dtype=str, chunksize=chunksize, low_memory=False):
        chunk = normalize_columns(chunk)
        if "VariationID" not in chunk.columns:
            raise RuntimeError(f"VariationID missing from {path}; columns={chunk.columns.tolist()[:30]}")
        ids = pd.to_numeric(chunk["VariationID"], errors="coerce")
        mask = ids.isin(target_ids)
        if not mask.any():
            continue
        sub = chunk.loc[mask, [c for c in wanted if c in chunk.columns]].copy()
        sub["VariationID"] = pd.to_numeric(sub["VariationID"], errors="coerce").astype("Int64")
        pieces.append(sub)
    if not pieces:
        return pd.DataFrame(columns=wanted)
    out = pd.concat(pieces, ignore_index=True)
    out = out[out["VariationID"].notna()].copy()
    out["VariationID"] = out["VariationID"].astype(int)
    disagreement = out.groupby("VariationID")["ClinicalSignificance"].nunique(dropna=False)
    chosen = choose_variant_row(out)
    chosen["assembly_significance_disagreement"] = chosen["VariationID"].map(disagreement).fillna(0).astype(int) > 1
    chosen["aggregate_side"] = chosen["ClinicalSignificance"].map(classification_side)
    review = chosen["ReviewStatus"] if "ReviewStatus" in chosen.columns else pd.Series("", index=chosen.index)
    chosen["review_tier"] = review.map(review_tier)
    return chosen


def read_submission_summary(path: Path, target_ids: set[int], chunksize: int = 400_000) -> pd.DataFrame:
    records: list[pd.DataFrame] = []
    use = ["VariationID", "ClinicalSignificance", "DateLastEvaluated", "ReviewStatus", "Submitter", "SCV", "CollectionMethod"]
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", dtype=str, chunksize=chunksize, low_memory=False):
        chunk = normalize_columns(chunk)
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
    flags = x["ClinicalSignificance"].map(significance_flags)
    x["_p"] = flags.map(lambda t: t[0])
    x["_b"] = flags.map(lambda t: t[1])
    x["_v"] = flags.map(lambda t: t[2])
    x["_conflict"] = flags.map(lambda t: t[3])
    review = x["ReviewStatus"] if "ReviewStatus" in x.columns else pd.Series("", index=x.index)
    x["_tier"] = review.map(review_tier)
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


def load_raw_and_scores(raw_zip: Path, score_zip: Path, predictions_path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    pred = pd.read_csv(predictions_path)
    needed_files = set(pred["protein_file"].astype(str))
    with zipfile.ZipFile(raw_zip) as z:
        members = [n for n in z.namelist() if n.lower().endswith(".csv")]
        if len(members) != 1:
            raise RuntimeError(f"Expected one raw clinical CSV, found {len(members)}")
        with z.open(members[0]) as f:
            raw = pd.read_csv(f, low_memory=False)
    raw = raw.reset_index(names="raw_index")
    raw["VariationID"] = pd.to_numeric(raw["ID"], errors="coerce").astype("Int64")
    raw["GeneSymbol"] = raw.get("GENEINFO", "").astype(str).str.split(":").str[0]

    pieces = []
    with zipfile.ZipFile(score_zip) as z:
        member_by_base = {Path(n).name: n for n in z.namelist() if n.lower().endswith(".csv")}
        missing_files = sorted(needed_files - set(member_by_base))
        if missing_files:
            raise RuntimeError(f"Missing {len(missing_files)} prediction files, examples={missing_files[:10]}")
        for i, fname in enumerate(sorted(needed_files), 1):
            with z.open(member_by_base[fname]) as f:
                d = pd.read_csv(f, usecols=["Unnamed: 0", "protein", "mutant", "DMS_bin_score"] + MODELS, low_memory=False)
            d = d.rename(columns={"Unnamed: 0": "raw_index"})
            d["protein_file"] = fname
            d["raw_index"] = pd.to_numeric(d["raw_index"], errors="coerce").astype("Int64")
            for m in MODELS:
                d[m] = pd.to_numeric(d[m], errors="coerce")
            pieces.append(d)
            if i % 500 == 0:
                print(f"scores {i}/{len(needed_files)}", flush=True)
    scores = pd.concat(pieces, ignore_index=True)
    scores = scores[scores["raw_index"].notna()].copy()
    scores["raw_index"] = scores["raw_index"].astype(int)
    if scores["raw_index"].duplicated().any():
        raise RuntimeError("Score raw_index is not unique")
    raw_lookup = raw.set_index("raw_index")
    if not set(scores["raw_index"]).issubset(set(raw_lookup.index)):
        raise RuntimeError("Score raw_index exceeds raw table")
    mapped = raw_lookup.loc[scores["raw_index"]].reset_index()
    protein_match = mapped["protein"].astype(str).to_numpy() == scores["protein"].astype(str).to_numpy()
    mutant_match = mapped["mutant"].astype(str).to_numpy() == scores["mutant"].astype(str).to_numpy()
    if protein_match.mean() < 0.999 or mutant_match.mean() < 0.999:
        raise RuntimeError(f"Raw/score mapping mismatch: protein={protein_match.mean()}, mutant={mutant_match.mean()}")
    score_raw = scores.merge(
        raw[["raw_index", "VariationID", "#CHROM", "POS", "REF", "ALT", "GeneSymbol", "Feature", "HGVSc", "HGVSp", "protein", "mutant", "CLNSIG", "CLNREVSTAT", "CLNVC", "Consequence"]],
        on=["raw_index", "protein", "mutant"], how="left", validate="one_to_one"
    )
    merged = pred.merge(score_raw, on=["protein_file", "mutant"], how="left", validate="one_to_one", suffixes=("", "_score"))
    missing = int(merged["VariationID"].isna().sum())
    merged = merged[merged["VariationID"].notna()].copy()
    merged["VariationID"] = merged["VariationID"].astype(int)
    merged["proteinGym_side"] = np.where(merged["label"].astype(int).eq(1), "P", "B")
    merged["raw_side"] = merged["CLNSIG"].map(classification_side)
    diagnostics = {
        "raw_rows": int(len(raw)),
        "score_rows": int(len(scores)),
        "frozen_prediction_rows": int(len(pred)),
        "mapped_rows": int(len(merged)),
        "mapping_missing": missing,
        "protein_match_fraction": float(protein_match.mean()),
        "mutant_match_fraction": float(mutant_match.mean()),
        "raw_label_agreement_fraction": float((merged["proteinGym_side"] == merged["raw_side"]).mean()),
        "unique_variation_ids": int(merged["VariationID"].nunique()),
    }
    return merged, diagnostics


def rank_scores_within_protein(df: pd.DataFrame) -> pd.DataFrame:
    x = df.copy()
    for m in MODELS:
        x[m + "_rank"] = x.groupby("protein_file")[m].rank(method="average", pct=True)
    return x


def safe_auc(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y, int); p = np.asarray(p, float)
    mask = np.isfinite(p)
    if mask.sum() < 4 or np.unique(y[mask]).size < 2:
        return float("nan")
    return float(roc_auc_score(y[mask], p[mask]))


def safe_ap(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y, int); p = np.asarray(p, float)
    mask = np.isfinite(p)
    if mask.sum() < 4 or np.unique(y[mask]).size < 2:
        return float("nan")
    return float(average_precision_score(y[mask], p[mask]))


def fit_temporal_recipe(train: pd.DataFrame) -> dict[str, Any]:
    signs: dict[str, int] = {}
    perf: dict[str, float] = {}
    usable_proteins = 0
    for m in MODELS:
        vals = []
        for _, g in train.groupby("protein_file", sort=False):
            auc = safe_auc(g["historical_label"].to_numpy(), g[m + "_rank"].to_numpy())
            if np.isfinite(auc):
                vals.append(auc)
        mean_auc = float(np.mean(vals)) if vals else float("nan")
        signs[m] = 1 if not np.isfinite(mean_auc) or mean_auc >= 0.5 else -1
        perf[m] = mean_auc if signs[m] > 0 else 1.0 - mean_auc
        usable_proteins = max(usable_proteins, len(vals))
    weights = np.array([max(perf[m] - 0.5, 0.001) for m in MODELS], dtype=float)
    weights /= weights.sum()

    def raw_scores(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        matrix = np.column_stack([
            frame[m + "_rank"].to_numpy(float) if signs[m] > 0 else 1.0 - frame[m + "_rank"].to_numpy(float)
            for m in MODELS
        ])
        ensemble = matrix @ weights
        poet = matrix[:, MODELS.index("PoET")]
        return ensemble, poet

    ens_train, poet_train = raw_scores(train)
    y = train["historical_label"].to_numpy(int)
    counts = train.groupby(["protein_file", "historical_label"]).size().to_dict()
    sw = np.array([0.5 / counts[(p, int(lbl))] for p, lbl in zip(train["protein_file"], train["historical_label"])], dtype=float)
    sw /= sw.mean()

    def fit_platt(score: np.ndarray) -> LogisticRegression:
        model = LogisticRegression(C=1e4, solver="lbfgs", max_iter=3000)
        model.fit(score.reshape(-1, 1), y, sample_weight=sw)
        return model

    return {
        "signs": signs,
        "performance": perf,
        "weights": {m: float(w) for m, w in zip(MODELS, weights)},
        "ensemble_platt": fit_platt(ens_train),
        "poet_platt": fit_platt(poet_train),
        "raw_score_function": raw_scores,
        "training_variants": int(len(train)),
        "training_proteins": int(train["protein_file"].nunique()),
        "training_auc_eligible_proteins": int(usable_proteins),
    }


def apply_temporal_recipe(model: dict[str, Any], test: pd.DataFrame) -> pd.DataFrame:
    out = test.copy()
    ens_raw, poet_raw = model["raw_score_function"](test)
    out["temporal_ensemble_raw"] = ens_raw
    out["temporal_poet_raw"] = poet_raw
    out["temporal_ensemble_prob"] = model["ensemble_platt"].predict_proba(ens_raw.reshape(-1, 1))[:, 1]
    out["temporal_poet_prob"] = model["poet_platt"].predict_proba(poet_raw.reshape(-1, 1))[:, 1]
    return out


def calibration_intercept_slope(y: np.ndarray, p: np.ndarray) -> dict[str, float | None]:
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    y = np.asarray(y, int)
    if len(y) < 20 or np.unique(y).size < 2:
        return {"intercept": None, "slope": None}
    x = np.log(p / (1 - p)).reshape(-1, 1)
    try:
        model = LogisticRegression(C=1e4, solver="lbfgs", max_iter=3000)
        model.fit(x, y)
        return {"intercept": float(model.intercept_[0]), "slope": float(model.coef_[0, 0])}
    except Exception:
        return {"intercept": None, "slope": None}


def metric_set(y: np.ndarray, p: np.ndarray) -> dict[str, Any]:
    p = np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6)
    y = np.asarray(y, int)
    return {
        "n": int(len(y)),
        "n_pathogenic": int(y.sum()),
        "n_benign": int((1 - y).sum()),
        "auroc": safe_auc(y, p),
        "auprc": safe_ap(y, p),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
        "accuracy_at_0_5": float(np.mean((p >= 0.5).astype(int) == y)),
        "calibration": calibration_intercept_slope(y, p),
    }


def cluster_bootstrap_compare(df: pd.DataFrame, ensemble_col: str, base_col: str, n_boot: int = 5000) -> dict[str, Any]:
    work = df[["cluster_id", "current_label", ensemble_col, base_col]].dropna().copy()
    groups = {cid: g for cid, g in work.groupby("cluster_id", sort=False)}
    ids = np.array(list(groups), dtype=object)
    if len(ids) < 3:
        return {"n_clusters": int(len(ids)), "error": "too few clusters"}
    y = work["current_label"].to_numpy(int)
    pe = np.clip(work[ensemble_col].to_numpy(float), 1e-6, 1 - 1e-6)
    pb = np.clip(work[base_col].to_numpy(float), 1e-6, 1 - 1e-6)
    observed = {
        "brier_improvement": float(brier_score_loss(y, pb) - brier_score_loss(y, pe)),
        "logloss_improvement": float(log_loss(y, pb, labels=[0,1]) - log_loss(y, pe, labels=[0,1])),
        "auroc_improvement": float(safe_auc(y, pe) - safe_auc(y, pb)),
    }
    rng = np.random.default_rng(SEED)
    brier = []
    logloss = []
    aucdiff = []
    for _ in range(n_boot):
        sampled = rng.choice(ids, size=len(ids), replace=True)
        boot = pd.concat([groups[cid] for cid in sampled], ignore_index=True)
        yy = boot["current_label"].to_numpy(int)
        ee = np.clip(boot[ensemble_col].to_numpy(float), 1e-6, 1-1e-6)
        bb = np.clip(boot[base_col].to_numpy(float), 1e-6, 1-1e-6)
        brier.append(brier_score_loss(yy, bb) - brier_score_loss(yy, ee))
        logloss.append(log_loss(yy, bb, labels=[0,1]) - log_loss(yy, ee, labels=[0,1]))
        a1, a0 = safe_auc(yy, ee), safe_auc(yy, bb)
        if np.isfinite(a1) and np.isfinite(a0):
            aucdiff.append(a1 - a0)
    def ci(values: list[float]) -> list[float | None]:
        arr = np.asarray(values, float)
        arr = arr[np.isfinite(arr)]
        if len(arr) < 100:
            return [None, None]
        return [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))]
    return {
        "n_clusters": int(len(ids)),
        "observed": observed,
        "cluster_bootstrap_95ci": {
            "brier_improvement": ci(brier),
            "logloss_improvement": ci(logloss),
            "auroc_improvement": ci(aucdiff),
        },
    }


def abstention_table(df: pd.DataFrame, probability_col: str) -> pd.DataFrame:
    rows = []
    x = df[["current_label", probability_col]].dropna().copy()
    x["confidence"] = (x[probability_col] - 0.5).abs()
    for coverage in [1.0, 0.9, 0.8, 0.7, 0.5]:
        n = max(1, int(math.floor(len(x) * coverage)))
        sub = x.nlargest(n, "confidence")
        y = sub["current_label"].to_numpy(int)
        p = np.clip(sub[probability_col].to_numpy(float), 1e-6, 1 - 1e-6)
        pred = (p >= 0.5).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        tn = int(((pred == 0) & (y == 0)).sum())
        fn = int(((pred == 0) & (y == 1)).sum())
        rows.append({
            "coverage": coverage, "n": len(sub), "error_rate": float(np.mean(pred != y)),
            "brier": float(brier_score_loss(y, p)),
            "ppv": float(tp / (tp + fp)) if tp + fp else None,
            "npv": float(tn / (tn + fn)) if tn + fn else None,
        })
    return pd.DataFrame(rows)


def make_status_table(base: pd.DataFrame, variant: pd.DataFrame, submission: pd.DataFrame, prefix: str) -> pd.DataFrame:
    vcols = [
        "VariationID", "ClinicalSignificance", "ReviewStatus", "review_tier", "aggregate_side",
        "NumberSubmitters", "LastEvaluated", "GeneSymbol", "Assembly", "Chromosome", "Start", "Stop",
        "ReferenceAllele", "AlternateAllele", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
    ]
    v = variant[[c for c in vcols if c in variant.columns]].copy()
    v = v.rename(columns={c: f"{prefix}_{c}" for c in v.columns if c != "VariationID"})
    s = submission.copy()
    s = s.rename(columns={c: f"{prefix}_{c}" for c in s.columns if c != "VariationID"})
    return base.merge(v, on="VariationID", how="left", validate="many_to_one").merge(s, on="VariationID", how="left", validate="many_to_one")


def historical_train_test(scored: pd.DataFrame, cutoff: str, current_prefix: str = "current") -> tuple[pd.DataFrame, pd.DataFrame]:
    hist = cutoff
    test_mask = (
        scored[f"{hist}_aggregate_side"].eq("VUS")
        & scored[f"{hist}_pure_vus_submissions"].fillna(False)
        & scored[f"{current_prefix}_aggregate_side"].isin(["P", "B"])
        & ~scored[f"{current_prefix}_resolved_submission_conflict"].fillna(True)
        & scored["proteinGym_side"].eq(scored[f"{current_prefix}_aggregate_side"])
    )
    test = scored[test_mask].copy()
    test["current_label"] = test[f"{current_prefix}_aggregate_side"].map({"B": 0, "P": 1}).astype(int)
    train_mask = (
        scored[f"{hist}_aggregate_side"].isin(["P", "B"])
        & ~scored[f"{hist}_resolved_submission_conflict"].fillna(True)
        & scored["proteinGym_side"].eq(scored[f"{hist}_aggregate_side"])
    )
    train = scored[train_mask].copy()
    train["historical_label"] = train[f"{hist}_aggregate_side"].map({"B": 0, "P": 1}).astype(int)
    return train, test


def truth_subset(test: pd.DataFrame, tier: str) -> pd.DataFrame:
    if tier == "primary_high_confidence":
        return test[test["current_review_tier"].fillna(0).astype(int) >= 2].copy()
    if tier == "secondary_one_star_or_higher":
        return test[test["current_review_tier"].fillna(0).astype(int) >= 1].copy()
    return test.copy()


def evaluate_cohort(df: pd.DataFrame, ensemble_prob: str, base_prob: str) -> dict[str, Any]:
    y = df["current_label"].to_numpy(int)
    return {
        "n_variants": int(len(df)),
        "n_proteins": int(df["protein_file"].nunique()),
        "n_clusters": int(df["cluster_id"].nunique()),
        "n_pathogenic": int(y.sum()),
        "n_benign": int((1-y).sum()),
        "ensemble": metric_set(y, df[ensemble_prob].to_numpy(float)),
        "poet": metric_set(y, df[base_prob].to_numpy(float)),
        "comparison": cluster_bootstrap_compare(df, ensemble_prob, base_prob),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-clinical-zip", required=True)
    ap.add_argument("--score-zip", required=True)
    ap.add_argument("--frozen-predictions", required=True)
    ap.add_argument("--protocol-lock", required=True)
    ap.add_argument("--variant", action="append", required=True, help="RELEASE=path")
    ap.add_argument("--submission", action="append", required=True, help="RELEASE=path")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    protocol = json.loads(Path(args.protocol_lock).read_text())
    if [x["release"] for x in protocol["historical_cutoffs"]] != CUTOFFS:
        raise RuntimeError("Protocol cutoff mismatch")

    variant_paths = dict(item.split("=", 1) for item in args.variant)
    submission_paths = dict(item.split("=", 1) for item in args.submission)
    required_releases = CUTOFFS + ["2026-08"]
    if set(required_releases) - set(variant_paths) or set(required_releases) - set(submission_paths):
        raise RuntimeError("All historical and current files are required")

    scored, mapping_diag = load_raw_and_scores(Path(args.raw_clinical_zip), Path(args.score_zip), Path(args.frozen_predictions))
    scored = rank_scores_within_protein(scored)
    target_ids = set(scored["VariationID"].astype(int))
    release_diag = {}
    for release in required_releases:
        print(f"Reading ClinVar {release} variant summary", flush=True)
        v = read_variant_summary(Path(variant_paths[release]), target_ids)
        ids_for_submission = set(v["VariationID"].astype(int)) if not v.empty else target_ids
        print(f"Reading ClinVar {release} submission summary", flush=True)
        s = read_submission_summary(Path(submission_paths[release]), ids_for_submission)
        prefix = "current" if release == "2026-08" else release
        scored = make_status_table(scored, v, s, prefix)
        release_diag[release] = {
            "variant_ids_found": int(v["VariationID"].nunique()) if not v.empty else 0,
            "submission_ids_found": int(s["VariationID"].nunique()) if not s.empty else 0,
            "aggregate_side_counts": v["aggregate_side"].value_counts(dropna=False).to_dict() if not v.empty else {},
            "assembly_significance_disagreements": int(v.get("assembly_significance_disagreement", pd.Series(dtype=bool)).sum()) if not v.empty else 0,
        }
        print(json.dumps({release: release_diag[release]}, indent=2), flush=True)

    all_results: dict[str, Any] = {
        "protocol": protocol,
        "mapping_diagnostics": mapping_diag,
        "release_diagnostics": release_diag,
        "cutoffs": {},
        "literature_boundary": "Historical VUS-to-future-resolution validation predates this study; this analysis tests the frozen ProteinGym ensemble, not the novelty of temporal validation itself."
    }
    cohort_frames = []
    abstention_frames = []

    for cutoff in CUTOFFS:
        train, test = historical_train_test(scored, cutoff)
        cutoff_result: dict[str, Any] = {
            "historical_resolved_training_variants": int(len(train)),
            "historical_resolved_training_proteins": int(train["protein_file"].nunique()),
            "pure_historical_vus_resolved_by_2026_all_tiers": int(len(test)),
            "truth_tiers": {},
        }
        temporal_model = None
        if len(train) >= 1000 and train["historical_label"].nunique() == 2:
            temporal_model = fit_temporal_recipe(train)
            test = apply_temporal_recipe(temporal_model, test)
            cutoff_result["temporal_model"] = {
                "signs": temporal_model["signs"],
                "performance": temporal_model["performance"],
                "weights": temporal_model["weights"],
                "training_variants": temporal_model["training_variants"],
                "training_proteins": temporal_model["training_proteins"],
            }
        else:
            cutoff_result["temporal_model_error"] = "Insufficient historical resolved training data"

        for tier in TRUTH_TIERS:
            cohort = truth_subset(test, tier)
            cohort["cutoff"] = cutoff; cohort["truth_tier"] = tier
            cohort_frames.append(cohort)
            tier_result: dict[str, Any] = {
                "frozen_existing": evaluate_cohort(cohort, "fixed_prob", "base_prob") if len(cohort) >= 2 and cohort["current_label"].nunique() == 2 else {"n_variants": int(len(cohort)), "error": "insufficient two-class cohort"}
            }
            if temporal_model is not None and len(cohort) >= 2 and cohort["current_label"].nunique() == 2:
                tier_result["temporal_clean"] = evaluate_cohort(cohort, "temporal_ensemble_prob", "temporal_poet_prob")
                abst = abstention_table(cohort, "temporal_ensemble_prob")
                abst["cutoff"] = cutoff; abst["truth_tier"] = tier; abst["model"] = "temporal_clean_ensemble"
                abstention_frames.append(abst)
            cutoff_result["truth_tiers"][tier] = tier_result
        all_results["cutoffs"][cutoff] = cutoff_result

    primary = all_results["cutoffs"]["2021-08"]["truth_tiers"]["primary_high_confidence"].get("temporal_clean", {})
    primary_n = int(primary.get("n_variants", 0))
    primary_cmp = primary.get("comparison", {})
    primary_obs = primary_cmp.get("observed", {})
    primary_ci = primary_cmp.get("cluster_bootstrap_95ci", {})
    brier_ci = primary_ci.get("brier_improvement", [None, None])
    log_ci = primary_ci.get("logloss_improvement", [None, None])
    proper_positive = (
        (brier_ci[0] is not None and brier_ci[0] > 0 and primary_obs.get("logloss_improvement", -1) >= 0)
        or (log_ci[0] is not None and log_ci[0] > 0 and primary_obs.get("brier_improvement", -1) >= 0)
    )
    cutoff_directions = []
    for cutoff in CUTOFFS:
        tr = all_results["cutoffs"][cutoff]["truth_tiers"]["secondary_one_star_or_higher"].get("temporal_clean", {})
        obs = tr.get("comparison", {}).get("observed", {})
        cutoff_directions.append(bool(obs.get("brier_improvement", 0) > 0 and obs.get("logloss_improvement", 0) > 0))
    replicated = sum(cutoff_directions) >= 2
    high_conf_direction = bool(primary_obs.get("brier_improvement", 0) > 0 and primary_obs.get("logloss_improvement", 0) > 0)
    success = bool(primary_n >= 100 and proper_positive and replicated and high_conf_direction)
    all_results["decision"] = {
        "primary_high_confidence_n": primary_n,
        "proper_score_superiority_confirmed": proper_positive,
        "positive_direction_in_at_least_two_cutoffs": replicated,
        "cutoff_directions": dict(zip(CUTOFFS, cutoff_directions)),
        "historical_time_machine_success": success,
        "level": "LEVEL_3_CANDIDATE" if success else "LEVEL_1_OR_2_ONLY",
        "grand_problem": "NOT_SOLVED",
    }

    cohorts = pd.concat(cohort_frames, ignore_index=True) if cohort_frames else pd.DataFrame()
    cohorts.to_csv(out / "historical_vus_cohorts.csv.gz", index=False, compression="gzip")
    if abstention_frames:
        pd.concat(abstention_frames, ignore_index=True).to_csv(out / "abstention_results.csv", index=False)
    rows = []
    for cutoff, cr in all_results["cutoffs"].items():
        for tier, tr in cr["truth_tiers"].items():
            for model_kind, res in tr.items():
                if "error" in res:
                    rows.append({"cutoff": cutoff, "truth_tier": tier, "model_kind": model_kind, "n": res.get("n_variants", 0), "error": res["error"]})
                    continue
                cmp = res["comparison"]
                rows.append({
                    "cutoff": cutoff, "truth_tier": tier, "model_kind": model_kind,
                    "n": res["n_variants"], "proteins": res["n_proteins"], "clusters": res["n_clusters"],
                    "ensemble_auc": res["ensemble"]["auroc"], "poet_auc": res["poet"]["auroc"],
                    "ensemble_brier": res["ensemble"]["brier"], "poet_brier": res["poet"]["brier"],
                    "ensemble_logloss": res["ensemble"]["log_loss"], "poet_logloss": res["poet"]["log_loss"],
                    "brier_improvement": cmp["observed"]["brier_improvement"],
                    "logloss_improvement": cmp["observed"]["logloss_improvement"],
                    "auc_improvement": cmp["observed"]["auroc_improvement"],
                    "brier_ci_low": cmp["cluster_bootstrap_95ci"]["brier_improvement"][0],
                    "brier_ci_high": cmp["cluster_bootstrap_95ci"]["brier_improvement"][1],
                    "logloss_ci_low": cmp["cluster_bootstrap_95ci"]["logloss_improvement"][0],
                    "logloss_ci_high": cmp["cluster_bootstrap_95ci"]["logloss_improvement"][1],
                    "auc_ci_low": cmp["cluster_bootstrap_95ci"]["auroc_improvement"][0],
                    "auc_ci_high": cmp["cluster_bootstrap_95ci"]["auroc_improvement"][1],
                })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "time_machine_summary.csv", index=False)
    (out / "FINAL_RESULTS.json").write_text(json.dumps(all_results, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    plot = summary[(summary["truth_tier"] == "secondary_one_star_or_higher") & (summary["model_kind"] == "temporal_clean")].copy()
    if not plot.empty and "auc_improvement" in plot:
        fig, ax = plt.subplots(figsize=(8, 5))
        y = np.arange(len(plot))
        lo = plot["auc_improvement"] - plot["auc_ci_low"]
        hi = plot["auc_ci_high"] - plot["auc_improvement"]
        ax.errorbar(plot["auc_improvement"], y, xerr=[lo, hi], fmt="o", capsize=4)
        ax.axvline(0, linewidth=1)
        ax.set_yticks(y, plot["cutoff"])
        ax.set_xlabel("AUROC improvement: temporal-clean ensemble minus PoET")
        ax.set_title("ClinVar historical VUS time-machine validation")
        fig.tight_layout(); fig.savefig(out / "FIG_PRIMARY_AUC.png", dpi=200); plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(plot)); width = 0.35
        ax.bar(x - width/2, plot["brier_improvement"], width, label="Brier improvement")
        ax.bar(x + width/2, plot["logloss_improvement"], width, label="Log-loss improvement")
        ax.axhline(0, linewidth=1)
        ax.set_xticks(x, plot["cutoff"])
        ax.set_ylabel("PoET loss minus ensemble loss; positive favors ensemble")
        ax.set_title("Proper-score temporal replication")
        ax.legend(); fig.tight_layout(); fig.savefig(out / "FIG_SENSITIVITY_PROPER_SCORES.png", dpi=200); plt.close(fig)

    report_lines = [
        "# ClinVar historical VUS time-machine validation", "",
        "- Primary cutoff: 2021-08", "- Replication cutoffs: 2022-08, 2023-08",
        "- Outcome release: 2026-08", "- Grand problem: **NOT SOLVED**", "",
        "## Decision", "", f"`historical_time_machine_success = {success}`", "",
        f"Primary high-confidence resolved historical VUS: **{primary_n}**", "",
        "## Interpretation boundary", "",
        "This analysis evaluates variants that were pure VUS in archived ClinVar releases and later became nonconflicting P/LP or B/LB records. The temporal-clean recipe is trained only on variants already resolved at each historical cutoff. The pre-existing frozen out-of-fold ensemble is reported as a supportive, not temporally pure, analysis.", "",
    ]
    if not summary.empty:
        report_lines += ["## Summary table", "", summary.to_markdown(index=False), ""]
    (out / "RESULT_REPORT.md").write_text("\n".join(report_lines), encoding="utf-8")
    print(json.dumps(all_results["decision"], indent=2), flush=True)


if __name__ == "__main__":
    main()
