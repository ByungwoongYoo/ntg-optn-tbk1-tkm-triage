#!/usr/bin/env python3
"""Efficient exact-subset comparison against all ProteinGym clinical baselines.

Reads each clinical score file once, computes per-gene AUCs for every available
zero-shot and clinically supervised predictor, and compares the unchanged frozen
DMS-trained ensemble with the strongest high-coverage zero-shot baseline.
"""
from __future__ import annotations

import argparse
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SEED = 20260817


def member_map(archive: zipfile.ZipFile) -> dict[str, str]:
    return {
        Path(name).name: name
        for name in archive.namelist()
        if name.lower().endswith(".csv") and not name.endswith("/")
    }


def score_specs(config: dict) -> dict[str, dict]:
    specs: dict[str, dict] = {}
    for group, category in [
        ("model_list_zero_shot_substitutions_clinical", "zero_shot"),
        ("model_list_supervised_substitutions_clinical", "clinically_supervised"),
    ]:
        for model, detail in config.get(group, {}).items():
            specs[model] = {
                "category": category,
                "direction": float(detail.get("directionality", 1)),
            }
    return specs


def rank01(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").rank(method="average", pct=True)


def bootstrap(values: np.ndarray, n: int = 50000) -> list[float]:
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(SEED)
    draws = np.empty(n, dtype=float)
    for index in range(n):
        draws[index] = np.mean(rng.choice(values, size=len(values), replace=True))
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def sign_flip(values: np.ndarray, n: int = 200000) -> float:
    values = values[np.isfinite(values)]
    observed = float(np.mean(values))
    rng = np.random.default_rng(SEED)
    count = 0
    for _ in range(n):
        null_value = float(np.mean(values * rng.choice([-1.0, 1.0], size=len(values))))
        count += int(null_value >= observed)
    return float((count + 1) / (n + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinical-scores", required=True)
    parser.add_argument("--clinical-reference", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--frozen-variants", required=True)
    parser.add_argument("--corrected-overlap", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(args.config).read_text())
    specs = score_specs(config)
    reference = pd.read_csv(args.clinical_reference)
    file_by_gene = dict(
        zip(reference["DMS_id"].astype(str), reference["DMS_filename"].astype(str))
    )
    frozen = pd.read_csv(args.frozen_variants, compression="gzip")
    frozen["DMS_id"] = frozen["DMS_id"].astype(str)
    frozen["mutant"] = frozen["mutant"].astype(str)
    overlap = pd.read_csv(args.corrected_overlap)
    unseen = set(
        overlap.loc[~overlap["seen_conservative_union"], "DMS_id"].astype(str)
    )
    frozen = frozen[frozen["DMS_id"].isin(unseen)].copy()

    rows: list[dict] = []
    failures: list[dict] = []
    with zipfile.ZipFile(args.clinical_scores) as archive:
        members = member_map(archive)
        for index, gene in enumerate(sorted(unseen), 1):
            target = frozen[frozen["DMS_id"] == gene][
                ["mutant", "benign_label", "prediction_ridge100"]
            ].copy()
            if target.empty:
                continue
            if min((target["benign_label"] == 0).sum(), (target["benign_label"] == 1).sum()) < 5:
                continue
            filename = file_by_gene.get(gene, f"{gene}.csv")
            member = members.get(filename) or members.get(f"{gene}.csv")
            if member is None:
                failures.append({"DMS_id": gene, "reason": "score file absent"})
                continue
            try:
                with archive.open(member) as handle:
                    header = pd.read_csv(handle, nrows=0).columns.tolist()
                model_columns = [model for model in specs if model in header]
                with archive.open(member) as handle:
                    scores = pd.read_csv(
                        handle,
                        usecols=["mutant"] + model_columns,
                        low_memory=False,
                    )
                scores["mutant"] = scores["mutant"].astype(str)
                scores = scores.drop_duplicates("mutant")
                merged = target.merge(scores, on="mutant", how="left")
                ensemble_auc = float(
                    roc_auc_score(merged["benign_label"], merged["prediction_ridge100"])
                )
                for model in model_columns:
                    raw = pd.to_numeric(merged[model], errors="coerce") * specs[model]["direction"]
                    valid = pd.DataFrame(
                        {"label": merged["benign_label"], "score": rank01(raw)}
                    ).dropna()
                    if min((valid["label"] == 0).sum(), (valid["label"] == 1).sum()) < 5:
                        continue
                    baseline_auc = float(roc_auc_score(valid["label"], valid["score"]))
                    rows.append(
                        {
                            "DMS_id": gene,
                            "model": model,
                            "category": specs[model]["category"],
                            "n_frozen_variants": int(len(merged)),
                            "n_model_variants": int(len(valid)),
                            "variant_coverage": float(len(valid) / len(merged)),
                            "ensemble_auc": ensemble_auc,
                            "baseline_auc": baseline_auc,
                            "difference": ensemble_auc - baseline_auc,
                        }
                    )
                if index % 50 == 0:
                    print(f"processed {index}/{len(unseen)} genes", flush=True)
            except Exception as exc:
                failures.append({"DMS_id": gene, "reason": repr(exc)})

    pairwise = pd.DataFrame(rows)
    summary_rows = []
    for (model, category), group in pairwise.groupby(["model", "category"]):
        summary_rows.append(
            {
                "model": model,
                "category": category,
                "n_genes": int(group["DMS_id"].nunique()),
                "gene_coverage": float(group["DMS_id"].nunique() / 700),
                "weighted_variant_coverage": float(
                    group["n_model_variants"].sum() / group["n_frozen_variants"].sum()
                ),
                "ensemble_mean_auc_on_common_genes": float(group["ensemble_auc"].mean()),
                "baseline_mean_auc": float(group["baseline_auc"].mean()),
                "ensemble_minus_baseline": float(group["difference"].mean()),
                "median_difference": float(group["difference"].median()),
                "fraction_genes_ensemble_better": float((group["difference"] > 0).mean()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    eligible = summary[
        (summary["category"] == "zero_shot")
        & (summary["gene_coverage"] >= 0.90)
        & (summary["weighted_variant_coverage"] >= 0.90)
    ].copy()
    if eligible.empty:
        raise RuntimeError("No zero-shot model passed the 90% gene/variant coverage gate")
    strongest = eligible.sort_values(
        ["baseline_mean_auc", "gene_coverage"], ascending=False
    ).iloc[0]
    primary_data = pairwise[pairwise["model"] == strongest["model"]].copy()
    differences = primary_data["difference"].to_numpy(float)
    primary = {
        "strongest_eligible_zero_shot_baseline": strongest["model"],
        "n_common_genes": int(len(primary_data)),
        "ensemble_mean_auc": float(primary_data["ensemble_auc"].mean()),
        "baseline_mean_auc": float(primary_data["baseline_auc"].mean()),
        "mean_difference": float(differences.mean()),
        "median_difference": float(np.median(differences)),
        "paired_bootstrap_95ci": bootstrap(differences),
        "one_sided_sign_flip_p": sign_flip(differences),
        "fraction_genes_ensemble_better": float((differences > 0).mean()),
    }
    primary["decision"] = (
        "ENSEMBLE_BEATS_STRONGEST_AVAILABLE_ZERO_SHOT_BASELINE"
        if primary["paired_bootstrap_95ci"][0] > 0
        and primary["one_sided_sign_flip_p"] < 0.05
        else "NO_ROBUST_ADVANCE_OVER_STRONGEST_AVAILABLE_ZERO_SHOT_BASELINE"
    )
    result = {
        "scope": "exact 700-gene conservative DMS-unseen frozen clinical set",
        "frozen_variant_rows": int(len(frozen)),
        "models_evaluated": int(len(summary)),
        "eligible_zero_shot_models": eligible["model"].tolist(),
        "primary_comparison": primary,
        "important_boundary": (
            "Clinically supervised predictors are descriptive only. This is a post-label "
            "benchmark of an unchanged frozen prediction file, not a new prospective clinical test."
        ),
        "failures_n": int(len(failures)),
    }
    pairwise.to_csv(out / "all_baseline_gene_pairwise.csv.gz", index=False, compression="gzip")
    summary.sort_values(["category", "baseline_mean_auc"], ascending=[True, False]).to_csv(
        out / "all_baseline_summary.csv", index=False
    )
    pd.DataFrame(failures).to_csv(out / "failures.csv", index=False)
    (out / "result.json").write_text(json.dumps(result, indent=2))
    p = primary
    (out / "REPORT.md").write_text(
        f"""# Frozen ensemble versus all available ProteinGym clinical baselines

- Scope: exact 700-gene conservative DMS-unseen set
- Strongest eligible zero-shot baseline: **{p['strongest_eligible_zero_shot_baseline']}**
- Common genes: **{p['n_common_genes']}**
- Frozen ensemble mean gene AUC: **{p['ensemble_mean_auc']:.4f}**
- Baseline mean gene AUC: **{p['baseline_mean_auc']:.4f}**
- Difference: **{p['mean_difference']:+.4f}**
- Paired bootstrap 95% CI: **[{p['paired_bootstrap_95ci'][0]:+.4f}, {p['paired_bootstrap_95ci'][1]:+.4f}]**
- One-sided sign-flip p: **{p['one_sided_sign_flip_p']:.6f}**
- Genes improved: **{100*p['fraction_genes_ensemble_better']:.1f}%**
- Decision: **`{p['decision']}`**

Clinically supervised models are not fair zero-shot comparators and are reported only descriptively.
"""
    )
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
