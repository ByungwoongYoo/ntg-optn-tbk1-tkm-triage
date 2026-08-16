#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

from proteingym_wide_v2 import (
    bootstrap_mean_ci,
    family,
    load_data,
    model_train_stats,
    rank01,
    stable_bucket,
)

SEED = 20260817
LABEL_CANDIDATES = [
    "DMS_score_bin", "label", "clinical_label", "ClinVar_label",
    "pathogenic", "target", "y"
]
META_EXCLUDE = {
    "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin",
    "DMS_score_bin_manual", "target_seq", "sequence", "wildtype",
    "wild_type", "wt", "fitness", "score", "label", "clinical_label",
    "ClinVar_label", "pathogenic", "target", "y"
}


def safe_auc(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y)
    p = np.asarray(p, dtype=float)
    m = np.isfinite(p) & pd.notna(y)
    y = y[m].astype(int)
    p = p[m]
    if len(y) < 4 or len(np.unique(y)) != 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def normalize_name(name: str) -> str:
    x = name.lower()
    x = re.sub(r"[^a-z0-9]+", "", x)
    aliases = {
        "tranceptevel": "tranceptevel",
        "tranceptevelarge": "tranceptevel",
        "eve": "eveensemble",
        "eveensemble": "eveensemble",
        "poet200m": "poet",
        "poet": "poet",
        "esm1b": "esm1b",
        "gemme": "gemme",
    }
    return aliases.get(x, x)


def zip_csv_map(z: zipfile.ZipFile) -> dict[str, str]:
    return {Path(n).name: n for n in z.namelist() if n.lower().endswith(".csv")}


def choose_label_column(df: pd.DataFrame) -> str | None:
    for c in LABEL_CANDIDATES:
        if c in df.columns:
            return c
    for c in df.columns:
        lc = str(c).lower()
        if "bin" in lc or "pathogen" in lc or "clinvar" in lc:
            vals = pd.to_numeric(df[c], errors="coerce").dropna().unique()
            if len(vals) and set(np.unique(vals)).issubset({0, 1}):
                return str(c)
    return None


def select_common_models(
    dms_models: list[str],
    clinical_columns: list[str],
    train_perf: dict[str, float],
    max_models: int = 10,
) -> tuple[list[tuple[str, str]], dict[str, Any]]:
    dms_norm: dict[str, list[str]] = {}
    for m in dms_models:
        dms_norm.setdefault(normalize_name(m), []).append(m)
    clinical_norm: dict[str, list[str]] = {}
    for c in clinical_columns:
        if c in META_EXCLUDE or str(c).lower().startswith("dms_"):
            continue
        clinical_norm.setdefault(normalize_name(str(c)), []).append(str(c))

    pairs = []
    for key in sorted(set(dms_norm).intersection(clinical_norm)):
        dm = max(dms_norm[key], key=lambda m: np.nan_to_num(train_perf.get(m, np.nan), nan=-9))
        cm = clinical_norm[key][0]
        pairs.append((dm, cm))

    # Keep the strongest DMS-trained representative per model family.
    best_by_family: dict[str, tuple[str, str, float]] = {}
    for dm, cm in pairs:
        f = family(dm)
        p = float(train_perf.get(dm, np.nan))
        if f not in best_by_family or p > best_by_family[f][2]:
            best_by_family[f] = (dm, cm, p)
    selected = sorted(best_by_family.values(), key=lambda x: (-np.nan_to_num(x[2], nan=-9), x[0]))[:max_models]
    return [(a, b) for a, b, _ in selected], {
        "all_exact_or_alias_pairs": pairs,
        "selected_with_performance": selected,
        "clinical_columns": clinical_columns,
    }


def load_clinical(
    data_zip_path: Path,
    score_zip_path: Path,
    selected_pairs: list[tuple[str, str]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    frames = []
    failures = []
    header_inventory = {}
    with zipfile.ZipFile(data_zip_path) as dz, zipfile.ZipFile(score_zip_path) as sz:
        dmap = zip_csv_map(dz)
        smap = zip_csv_map(sz)
        common = sorted(set(dmap).intersection(smap))
        for i, fn in enumerate(common, 1):
            try:
                with dz.open(dmap[fn]) as f:
                    lab = pd.read_csv(f, low_memory=False)
                with sz.open(smap[fn]) as f:
                    scr = pd.read_csv(f, low_memory=False)
                label_col = choose_label_column(lab)
                if label_col is None and choose_label_column(scr) is not None:
                    label_col = choose_label_column(scr)
                    lab = scr
                if label_col is None or "mutant" not in lab.columns or "mutant" not in scr.columns:
                    failures.append({"file": fn, "reason": "missing mutant or binary label"})
                    continue
                header_inventory.setdefault("data_columns", list(lab.columns))
                header_inventory.setdefault("score_columns", list(scr.columns))
                base = lab[["mutant", label_col]].copy().rename(columns={label_col: "label"})
                base["mutant"] = base["mutant"].astype(str)
                base["label"] = pd.to_numeric(base["label"], errors="coerce")
                base = base[base["label"].isin([0, 1])].drop_duplicates("mutant")
                score_cols = [cm for _, cm in selected_pairs if cm in scr.columns]
                if len(score_cols) < 3:
                    failures.append({"file": fn, "reason": f"only {len(score_cols)} selected score columns"})
                    continue
                s = scr[["mutant"] + score_cols].copy()
                s["mutant"] = s["mutant"].astype(str)
                s = s.drop_duplicates("mutant")
                merged = base.merge(s, on="mutant", how="inner")
                if merged["label"].nunique() != 2 or len(merged) < 10:
                    failures.append({"file": fn, "reason": "fewer than 10 binary-labeled merged variants or one class"})
                    continue
                merged["protein_file"] = fn
                frames.append(merged)
            except Exception as exc:
                failures.append({"file": fn, "reason": repr(exc)})
            if i % 250 == 0:
                print(f"clinical files {i}/{len(common)} accepted={len(frames)}", flush=True)
    if not frames:
        raise RuntimeError(f"No clinical proteins loaded. Failures: {failures[:10]}")
    return pd.concat(frames, ignore_index=True), {
        "data_files": len(dmap), "score_files": len(smap), "common_files": len(common),
        "accepted_proteins": len(frames), "failures": failures, "header_inventory": header_inventory,
    }


def protein_auc_table(df: pd.DataFrame, prediction_cols: list[str]) -> pd.DataFrame:
    rows = []
    for protein, g in df.groupby("protein_file", sort=False):
        row = {
            "protein_file": protein,
            "n_variants": len(g),
            "n_pathogenic": int((g["label"] == 1).sum()),
            "n_benign": int((g["label"] == 0).sum()),
        }
        for c in prediction_cols:
            row[c] = safe_auc(g["label"].to_numpy(), g[c].to_numpy())
        rows.append(row)
    return pd.DataFrame(rows)


def comparison_stats(table: pd.DataFrame, a: str, b: str) -> dict[str, Any]:
    d = (table[a] - table[b]).to_numpy(float)
    d = d[np.isfinite(d)]
    try:
        p = float(wilcoxon(d, alternative="greater").pvalue) if len(d) >= 5 and np.any(d != 0) else None
    except Exception:
        p = None
    return {
        "n_proteins": int(len(d)),
        "ensemble_mean_protein_auc": float(table[a].mean()),
        "comparator_mean_protein_auc": float(table[b].mean()),
        "paired_mean_auc_gain": float(np.mean(d)),
        "paired_gain_95ci": bootstrap_mean_ci(d),
        "fraction_proteins_improved": float(np.mean(d > 0)),
        "wilcoxon_one_sided_p": p,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dms-score-zip", required=True)
    ap.add_argument("--dms-reference", required=True)
    ap.add_argument("--clinical-data-zip", required=True)
    ap.add_argument("--clinical-score-zip", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dms, dms_models, dms_diag = load_data(Path(args.dms_score_zip), Path(args.dms_reference), human_only=False)
    human = dms[dms["taxon"].astype(str).str.lower().eq("human")].copy()
    proteins = sorted(human["UniProt_ID"].unique())
    bucket = {p: stable_bucket(p) for p in proteins}
    human["split"] = human["UniProt_ID"].map(lambda x: "train" if bucket[x] <= 5 else ("tune" if bucket[x] <= 7 else "final"))
    split_counts = human[["DMS_id", "UniProt_ID", "split"]].drop_duplicates().groupby("split").agg(
        assays=("DMS_id", "nunique"), proteins=("UniProt_ID", "nunique")
    ).to_dict("index")
    if any(split_counts.get(s, {}).get("proteins", 0) < 6 for s in ["train", "tune", "final"]):
        assignments = {}
        for i, p in enumerate(proteins):
            frac = i / max(1, len(proteins))
            assignments[p] = "train" if frac < 0.6 else ("tune" if frac < 0.8 else "final")
        human["split"] = human["UniProt_ID"].map(assignments)
    train_tune = human[human["split"].isin(["train", "tune"])].copy()
    signs, perf = model_train_stats(train_tune, dms_models)

    with zipfile.ZipFile(args.clinical_score_zip) as sz:
        smap = zip_csv_map(sz)
        first = next(iter(smap.values()))
        with sz.open(first) as f:
            clinical_header = pd.read_csv(f, nrows=0).columns.tolist()

    selected_pairs, pair_diag = select_common_models(dms_models, clinical_header, perf, max_models=10)
    if len(selected_pairs) < 3:
        raise RuntimeError(f"Only {len(selected_pairs)} cross-benchmark model pairs: {selected_pairs}; clinical header={clinical_header}")

    clinical, clinical_diag = load_clinical(
        Path(args.clinical_data_zip), Path(args.clinical_score_zip), selected_pairs
    )

    # Convert every model score to an oriented within-protein fitness percentile using only DMS-trained signs.
    selected_dms = [dm for dm, _ in selected_pairs]
    selected_clin = [cm for _, cm in selected_pairs]
    weights_raw = np.array([max(float(perf[dm]), 0.001) ** 2 for dm in selected_dms])
    weights = weights_raw / weights_raw.sum()

    for dm, cm in selected_pairs:
        oriented_col = f"fit__{dm}"
        clinical[oriented_col] = clinical.groupby("protein_file")[cm].transform(rank01)
        if signs[dm] < 0:
            clinical[oriented_col] = 1.0 - clinical[oriented_col]
        clinical[f"risk__{dm}"] = 1.0 - clinical[oriented_col]

    fit_cols = [f"fit__{dm}" for dm in selected_dms]
    X = clinical[fit_cols].to_numpy(float)
    ok = np.isfinite(X)
    W = ok * weights.reshape(1, -1)
    denom = W.sum(axis=1)
    ensemble_fit = np.full(len(clinical), np.nan)
    valid = denom > 0
    ensemble_fit[valid] = np.nansum(X[valid] * weights.reshape(1, -1), axis=1) / denom[valid]
    clinical["ensemble_risk"] = 1.0 - ensemble_fit

    pred_cols = ["ensemble_risk"] + [f"risk__{dm}" for dm in selected_dms]
    aucs = protein_auc_table(clinical, pred_cols)
    summary_rows = []
    for c in pred_cols:
        summary_rows.append({
            "model": c,
            "mean_protein_auc": float(aucs[c].mean()),
            "median_protein_auc": float(aucs[c].median()),
            "n_proteins": int(aucs[c].notna().sum()),
            "pooled_auc": safe_auc(clinical["label"].to_numpy(), clinical[c].to_numpy()),
        })
    summary = pd.DataFrame(summary_rows).sort_values("mean_protein_auc", ascending=False)
    best_posthoc = str(summary[summary["model"] != "ensemble_risk"].iloc[0]["model"])
    comp = comparison_stats(aucs, "ensemble_risk", best_posthoc)

    # Compare against every accessible individual on the exact same clinical variant rows.
    pairwise = {}
    for c in pred_cols[1:]:
        pairwise[c] = comparison_stats(aucs, "ensemble_risk", c)

    result = {
        "question": "Does a DMS-trained frozen rank ensemble improve clinical pathogenicity discrimination?",
        "dms_training_scope": {
            "human_train_tune_assays": int(train_tune["DMS_id"].nunique()),
            "human_train_tune_proteins": int(train_tune["UniProt_ID"].nunique()),
        },
        "selected_model_pairs": selected_pairs,
        "weights": {dm: float(w) for dm, w in zip(selected_dms, weights)},
        "dms_signs": {dm: int(signs[dm]) for dm in selected_dms},
        "clinical_diagnostics": clinical_diag,
        "clinical_proteins_scored": int(aucs["protein_file"].nunique()),
        "clinical_variants_scored": int(len(clinical)),
        "ensemble_mean_protein_auc": float(aucs["ensemble_risk"].mean()),
        "posthoc_best_accessible_individual": best_posthoc,
        "posthoc_best_mean_protein_auc": float(aucs[best_posthoc].mean()),
        "posthoc_comparison": comp,
        "pairwise_comparisons": pairwise,
        "clinical_transfer_advance_confirmed": bool(comp["paired_gain_95ci"][0] is not None and comp["paired_gain_95ci"][0] > 0),
        "open_problem_fully_solved": False,
        "pair_diagnostics": pair_diag,
        "dms_diagnostics": dms_diag,
    }

    aucs.to_csv(out / "clinical_protein_auc.csv", index=False)
    summary.to_csv(out / "clinical_model_summary.csv", index=False)
    clinical[["protein_file", "mutant", "label", "ensemble_risk"] + [f"risk__{dm}" for dm in selected_dms]].to_csv(
        out / "clinical_predictions.csv.gz", index=False, compression="gzip"
    )
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    ci = comp["paired_gain_95ci"]
    report = f"""# ProteinGym clinical substitution transfer test v5

## Locked question

A model-weighting rule learned only from human DMS assays was transferred to the independent ProteinGym clinical substitution benchmark. No clinical label was used to select constituent models, signs, or weights.

## Data

- Clinical proteins with both classes and usable scores: {result['clinical_proteins_scored']}
- Clinical missense variants: {result['clinical_variants_scored']:,}
- Cross-benchmark constituent models: {', '.join(selected_dms)}

## Result

- Ensemble mean protein AUC: {result['ensemble_mean_protein_auc']:.4f}
- Best accessible individual selected after viewing clinical labels: `{best_posthoc}`
- Best individual mean protein AUC: {result['posthoc_best_mean_protein_auc']:.4f}
- Paired mean AUC gain: {comp['paired_mean_auc_gain']:+.4f}
- Paired bootstrap 95% CI: {ci}
- Proteins improved: {100 * comp['fraction_proteins_improved']:.1f}%
- One-sided Wilcoxon p: {comp['wilcoxon_one_sided_p']}

`clinical_transfer_advance_confirmed = {result['clinical_transfer_advance_confirmed']}`

This test evaluates clinical discrimination, not calibration, penetrance, inheritance, or patient-level diagnosis. A positive result is a benchmark advance; it is not a complete solution to all human missense variants.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
