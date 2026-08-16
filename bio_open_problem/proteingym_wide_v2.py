#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20260817
META_EXCLUDE = {
    "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin",
    "DMS_score_bin_manual", "target_seq", "sequence", "wildtype",
    "wild_type", "wt", "fitness", "score"
}


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 20 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return float("nan")
    return float(spearmanr(x[mask], y[mask]).statistic)


def rank01(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").rank(method="average", pct=True)


def stable_bucket(text: str, modulo: int = 10) -> int:
    h = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(h[:12], 16) % modulo


def family(name: str) -> str:
    n = name.lower()
    rules = [
        (r"aido", "AIDO"), (r"venusrem", "VenusREM"),
        (r"prosst", "ProSST"), (r"s[23]f", "S2F_S3F"),
        (r"trancept|tranception", "Tranception"),
        (r"deepsequence|eve|evmutation", "EvolutionaryVAE"),
        (r"esm1v|esm1b|esm2|esm3|esmc", "ESM"),
        (r"saprot", "SaProt"), (r"poet", "PoET"),
        (r"gemme", "GEMME"), (r"vespa", "VESPA"),
        (r"rsalor", "RSALOR"), (r"s s?ite|siterm", "SiteRM"),
        (r"proteinmpnn|esm-if", "InverseFolding"),
        (r"protssn", "ProtSSN"), (r"progen", "ProGen"),
        (r"xtrimopglm", "xTrimoPGLM"), (r"carp", "CARP"),
        (r"mif", "MIF"), (r"unirep", "UniRep"),
        (r"msa_transformer", "MSATransformer"),
        (r"site_independent", "SiteIndependent"),
        (r"wavenet", "WaveNet"), (r"rita", "RITA"),
        (r"escott", "ESCOTT"), (r"protgpt", "ProtGPT"),
    ]
    for pat, fam in rules:
        if re.search(pat, n):
            return fam
    return re.split(r"[_\- ]", name)[0]


def bootstrap_mean_ci(values: np.ndarray, n: int = 20000) -> list[float | None]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 3:
        return [None, None]
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(values), size=(n, len(values)))
    means = values[idx].mean(axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def assay_mean_performance(df: pd.DataFrame, model: str, target: str = "target_rank") -> float:
    vals = []
    for _, g in df.groupby("DMS_id", sort=False):
        vals.append(safe_spearman(g[model].to_numpy(), g[target].to_numpy()))
    return float(np.nanmean(vals))


def model_train_stats(df: pd.DataFrame, models: list[str]) -> tuple[dict[str, int], dict[str, float]]:
    signs: dict[str, int] = {}
    perf: dict[str, float] = {}
    for m in models:
        vals = []
        for _, g in df.groupby("DMS_id", sort=False):
            rho = safe_spearman(g[m].to_numpy(), g["target_rank"].to_numpy())
            vals.append(rho)
        med = float(np.nanmedian(vals)) if np.isfinite(vals).any() else float("nan")
        signs[m] = 1 if not np.isfinite(med) or med >= 0 else -1
        oriented = [v * signs[m] if np.isfinite(v) else np.nan for v in vals]
        perf[m] = float(np.nanmean(oriented))
    return signs, perf


def select_diverse(perf: dict[str, float], k: int) -> list[str]:
    by_family: dict[str, tuple[str, float]] = {}
    for m, p in perf.items():
        f = family(m)
        if f not in by_family or p > by_family[f][1]:
            by_family[f] = (m, p)
    ranked = sorted(by_family.values(), key=lambda x: (-np.nan_to_num(x[1], nan=-9), x[0]))
    return [m for m, _ in ranked[:k]]


def oriented_matrix(df: pd.DataFrame, selected: list[str], signs: dict[str, int]) -> np.ndarray:
    cols = []
    for m in selected:
        x = pd.to_numeric(df[m], errors="coerce").to_numpy(float)
        if signs[m] < 0:
            x = 1.0 - x
        cols.append(x)
    return np.column_stack(cols)


def weighted_row_mean(X: np.ndarray, weights: np.ndarray) -> np.ndarray:
    ok = np.isfinite(X)
    w = ok * weights.reshape(1, -1)
    denom = w.sum(axis=1)
    out = np.full(len(X), np.nan)
    valid = denom > 0
    out[valid] = np.nansum(X[valid] * weights.reshape(1, -1), axis=1) / denom[valid]
    return out


def feature_matrix(df: pd.DataFrame, selected: list[str], signs: dict[str, int]) -> np.ndarray:
    X = oriented_matrix(df, selected, signs)
    miss = (~np.isfinite(X)).astype(float)
    mean = np.nanmean(X, axis=1, keepdims=True)
    median = np.nanmedian(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, keepdims=True)
    n = np.isfinite(X).sum(axis=1, keepdims=True) / max(1, X.shape[1])
    return np.column_stack([X, miss, mean, median, sd, n])


def balanced_sample(df: pd.DataFrame, max_per_assay: int = 1200) -> pd.DataFrame:
    chunks = []
    for _, g in df.groupby("DMS_id", sort=False):
        if len(g) > max_per_assay:
            chunks.append(g.sample(max_per_assay, random_state=SEED))
        else:
            chunks.append(g)
    return pd.concat(chunks, ignore_index=False)


@dataclass(frozen=True)
class Spec:
    name: str
    kind: str
    k: int
    power: float = 1.0
    alpha: float = 25.0


SPECS = [
    Spec("median5", "median", 5), Spec("median10", "median", 10), Spec("median20", "median", 20),
    Spec("mean5", "mean", 5), Spec("mean10", "mean", 10), Spec("mean20", "mean", 20),
    Spec("weighted5_p1", "weighted", 5, 1.0), Spec("weighted10_p1", "weighted", 10, 1.0),
    Spec("weighted20_p1", "weighted", 20, 1.0), Spec("weighted10_p2", "weighted", 10, 2.0),
    Spec("weighted20_p2", "weighted", 20, 2.0), Spec("weighted10_p4", "weighted", 10, 4.0),
    Spec("ridge10_a10", "ridge", 10, alpha=10.0), Spec("ridge20_a10", "ridge", 20, alpha=10.0),
    Spec("ridge20_a100", "ridge", 20, alpha=100.0), Spec("hgb10", "hgb", 10),
    Spec("hgb20", "hgb", 20), Spec("extra10", "extra", 10), Spec("extra20", "extra", 20),
]


class FittedMethod:
    def __init__(self, spec: Spec, selected: list[str], signs: dict[str, int], perf: dict[str, float], model: Any = None):
        self.spec = spec
        self.selected = selected
        self.signs = signs
        self.perf = perf
        self.model = model

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        X = oriented_matrix(df, self.selected, self.signs)
        if self.spec.kind == "median":
            return np.nanmedian(X, axis=1)
        if self.spec.kind == "mean":
            return np.nanmean(X, axis=1)
        if self.spec.kind == "weighted":
            raw = np.array([max(self.perf[m], 0.001) ** self.spec.power for m in self.selected])
            return weighted_row_mean(X, raw / raw.sum())
        F = feature_matrix(df, self.selected, self.signs)
        return self.model.predict(F)


def fit_method(train: pd.DataFrame, all_models: list[str], spec: Spec) -> FittedMethod:
    signs, perf = model_train_stats(train, all_models)
    selected = select_diverse(perf, spec.k)
    model = None
    if spec.kind in {"ridge", "hgb", "extra"}:
        sample = balanced_sample(train, max_per_assay=1200)
        F = feature_matrix(sample, selected, signs)
        y = sample["target_rank"].to_numpy(float)
        if spec.kind == "ridge":
            model = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), Ridge(alpha=spec.alpha))
        elif spec.kind == "hgb":
            model = make_pipeline(SimpleImputer(strategy="median"), HistGradientBoostingRegressor(
                max_iter=250, learning_rate=0.04, max_leaf_nodes=15,
                min_samples_leaf=80, l2_regularization=3.0, random_state=SEED))
        else:
            model = make_pipeline(SimpleImputer(strategy="median"), ExtraTreesRegressor(
                n_estimators=300, min_samples_leaf=40, max_features=0.75,
                n_jobs=-1, random_state=SEED))
        model.fit(F, y)
    return FittedMethod(spec, selected, signs, perf, model)


def per_assay_scores(df: pd.DataFrame, prediction: np.ndarray, name: str) -> pd.DataFrame:
    tmp = df[["DMS_id", "UniProt_ID", "taxon", "selection_type", "msa_depth", "target_rank"]].copy()
    tmp["prediction"] = prediction
    rows = []
    for assay, g in tmp.groupby("DMS_id", sort=False):
        rows.append({
            "DMS_id": assay,
            "UniProt_ID": g["UniProt_ID"].iloc[0],
            "taxon": g["taxon"].iloc[0],
            "selection_type": g["selection_type"].iloc[0],
            "msa_depth": g["msa_depth"].iloc[0],
            "n_mutants": len(g),
            name: safe_spearman(g["prediction"].to_numpy(), g["target_rank"].to_numpy()),
        })
    return pd.DataFrame(rows)


def evaluate(df: pd.DataFrame, fitted: FittedMethod, name: str) -> tuple[pd.DataFrame, float]:
    p = fitted.predict(df)
    scores = per_assay_scores(df, p, name)
    return scores, float(scores[name].mean())


def load_data(score_zip_path: Path, reference_path: Path, human_only: bool) -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    ref = pd.read_csv(reference_path)
    if human_only:
        ref = ref[ref["taxon"].astype(str).str.lower().eq("human")].copy()
    ref = ref[(ref["DMS_number_single_mutants"] >= 100) & (~ref["includes_multiple_mutants"].astype(bool))].copy()
    assay_names = set(ref["DMS_filename"].astype(str))
    ref_by_file = ref.set_index("DMS_filename")
    with zipfile.ZipFile(score_zip_path) as z:
        members = {Path(n).name: n for n in z.namelist() if n.lower().endswith(".csv")}
        usable_files = sorted(assay_names.intersection(members))
        if len(usable_files) < 15:
            raise RuntimeError(f"Only {len(usable_files)} official score files matched the reference table")
        header_counts: dict[str, int] = {}
        for fn in usable_files:
            with z.open(members[fn]) as f:
                cols = pd.read_csv(f, nrows=0).columns.tolist()
            for c in cols:
                if c not in META_EXCLUDE and not c.lower().startswith("dms_"):
                    header_counts[c] = header_counts.get(c, 0) + 1
        min_assays = max(10, math.ceil(0.80 * len(usable_files)))
        models = sorted([c for c, n in header_counts.items() if n >= min_assays])
        frames = []
        failures = []
        for i, fn in enumerate(usable_files, 1):
            meta = ref_by_file.loc[fn]
            try:
                with z.open(members[fn]) as f:
                    d = pd.read_csv(f, low_memory=False)
                keep_models = [m for m in models if m in d.columns]
                keep = ["mutant", "DMS_score"] + keep_models
                d = d[keep].copy()
                d["mutant"] = d["mutant"].astype(str)
                d = d[~d["mutant"].str.contains(":", regex=False)]
                d["DMS_score"] = pd.to_numeric(d["DMS_score"], errors="coerce")
                d = d.dropna(subset=["DMS_score"]).drop_duplicates("mutant")
                if len(d) < 100:
                    failures.append({"file": fn, "reason": "fewer than 100 usable single substitutions"})
                    continue
                d["target_rank"] = rank01(d["DMS_score"])
                for m in models:
                    d[m] = rank01(d[m]) if m in d.columns else np.nan
                d["DMS_id"] = meta["DMS_id"]
                d["UniProt_ID"] = meta["UniProt_ID"]
                d["taxon"] = meta["taxon"]
                d["selection_type"] = meta.get("coarse_selection_type", meta.get("selection_type", "unknown"))
                d["msa_depth"] = meta.get("MSA_Neff_L_category", "unknown")
                frames.append(d[["DMS_id", "UniProt_ID", "taxon", "selection_type", "msa_depth", "mutant", "DMS_score", "target_rank"] + models])
                print(f"loaded {i}/{len(usable_files)} assays; accepted={len(frames)}", flush=True)
            except Exception as exc:
                failures.append({"file": fn, "reason": repr(exc)})
    data = pd.concat(frames, ignore_index=True)
    # Drop score columns with very poor mutation-level coverage after merging.
    coverage = data[models].notna().mean().sort_values(ascending=False)
    models = coverage[coverage >= 0.55].index.tolist()
    data = data[["DMS_id", "UniProt_ID", "taxon", "selection_type", "msa_depth", "mutant", "DMS_score", "target_rank"] + models]
    diagnostics = {
        "reference_assays": int(len(ref)), "loaded_assays": int(data["DMS_id"].nunique()),
        "unique_uniprot": int(data["UniProt_ID"].nunique()), "mutants": int(len(data)),
        "candidate_model_columns": len(models), "model_coverage": coverage.loc[models].to_dict(),
        "failures": failures,
    }
    return data, models, diagnostics


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-zip", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--human-only", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    data, models, diagnostics = load_data(Path(args.score_zip), Path(args.reference), args.human_only)
    proteins = sorted(data["UniProt_ID"].unique())
    bucket = {p: stable_bucket(p) for p in proteins}
    data["split"] = data["UniProt_ID"].map(lambda x: "train" if bucket[x] <= 5 else ("tune" if bucket[x] <= 7 else "final"))
    split_counts = data[["DMS_id", "UniProt_ID", "split"]].drop_duplicates().groupby("split").agg(assays=("DMS_id", "nunique"), proteins=("UniProt_ID", "nunique")).to_dict("index")
    if any(split_counts.get(s, {}).get("proteins", 0) < 6 for s in ["train", "tune", "final"]):
        # Deterministic fallback with roughly 60/20/20 proteins.
        assignments = {}
        for i, p in enumerate(proteins):
            frac = i / max(1, len(proteins))
            assignments[p] = "train" if frac < 0.6 else ("tune" if frac < 0.8 else "final")
        data["split"] = data["UniProt_ID"].map(assignments)
        split_counts = data[["DMS_id", "UniProt_ID", "split"]].drop_duplicates().groupby("split").agg(assays=("DMS_id", "nunique"), proteins=("UniProt_ID", "nunique")).to_dict("index")

    train = data[data["split"] == "train"].copy()
    tune = data[data["split"] == "tune"].copy()
    final = data[data["split"] == "final"].copy()

    tuning_rows = []
    fitted_by_name: dict[str, FittedMethod] = {}
    for spec in SPECS:
        fitted = fit_method(train, models, spec)
        scores, mean_score = evaluate(tune, fitted, spec.name)
        tuning_rows.append({
            "method": spec.name, "kind": spec.kind, "k": spec.k,
            "tune_mean_assay_spearman": mean_score,
            "tune_median_assay_spearman": float(scores[spec.name].median()),
            "selected_models": ";".join(fitted.selected),
        })
        fitted_by_name[spec.name] = fitted
        print(f"tune {spec.name}: {mean_score:.5f}", flush=True)
    tune_table = pd.DataFrame(tuning_rows).sort_values(["tune_mean_assay_spearman", "method"], ascending=[False, True])
    chosen_name = str(tune_table.iloc[0]["method"])
    chosen_spec = next(s for s in SPECS if s.name == chosen_name)

    train_tune = data[data["split"].isin(["train", "tune"])].copy()
    final_fit = fit_method(train_tune, models, chosen_spec)
    final_scores, final_mean = evaluate(final, final_fit, "ensemble")

    # Fair nested individual baseline: choose one model on train+tune, orient it there, then evaluate final.
    base_signs, base_perf = model_train_stats(train_tune, models)
    best_base = max(base_perf, key=lambda m: np.nan_to_num(base_perf[m], nan=-9))
    base_pred = pd.to_numeric(final[best_base], errors="coerce").to_numpy(float)
    if base_signs[best_base] < 0:
        base_pred = 1.0 - base_pred
    base_scores = per_assay_scores(final, base_pred, "best_base")
    compare = final_scores.merge(base_scores[["DMS_id", "best_base"]], on="DMS_id", how="inner")
    diffs = (compare["ensemble"] - compare["best_base"]).to_numpy(float)
    diffs = diffs[np.isfinite(diffs)]
    try:
        wilcox = float(wilcoxon(diffs, alternative="greater").pvalue) if len(diffs) >= 5 and np.any(diffs != 0) else None
    except Exception:
        wilcox = None

    # Named leaderboard models on the same sealed final subset.
    named_rows = []
    for m in ["VenusREM", "ProSST-2048", "ProSST-4096", "S3F_MSA", "S2F_MSA", "ESM3", "VespaG", "SaProt_650M_AF2", "GEMME", "TranceptEVE_L", "ProteinMPNN"]:
        if m in models:
            p = pd.to_numeric(final[m], errors="coerce").to_numpy(float)
            if base_signs[m] < 0:
                p = 1.0 - p
            s = per_assay_scores(final, p, m)
            named_rows.append({"model": m, "final_mean_assay_spearman": float(s[m].mean()), "n_assays": int(s[m].notna().sum())})

    result = {
        "scope": "human single substitutions" if args.human_only else "all taxa single substitutions",
        "diagnostics": diagnostics,
        "split_counts": split_counts,
        "chosen_method": chosen_name,
        "chosen_kind": chosen_spec.kind,
        "chosen_models": final_fit.selected,
        "final_mean_assay_spearman": final_mean,
        "final_median_assay_spearman": float(final_scores["ensemble"].median()),
        "best_individual_selected_without_final_labels": best_base,
        "best_individual_final_mean_assay_spearman": float(compare["best_base"].mean()),
        "paired_mean_improvement": float(np.mean(diffs)),
        "paired_improvement_95ci": bootstrap_mean_ci(diffs),
        "fraction_final_assays_improved": float(np.mean(diffs > 0)),
        "paired_wilcoxon_one_sided_p": wilcox,
        "benchmark_advance_confirmed": bool(bootstrap_mean_ci(diffs)[0] is not None and bootstrap_mean_ci(diffs)[0] > 0),
        "open_problem_fully_solved": False,
        "named_model_final_results": named_rows,
    }

    tune_table.to_csv(out / "tuning_results.csv", index=False)
    compare.to_csv(out / "sealed_final_assay_scores.csv", index=False)
    pd.DataFrame(named_rows).sort_values("final_mean_assay_spearman", ascending=False).to_csv(out / "named_model_final_results.csv", index=False)
    data[["DMS_id", "UniProt_ID", "mutant", "split"]].to_csv(out / "split_assignment.csv.gz", index=False, compression="gzip")
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8")

    ci = result["paired_improvement_95ci"]
    report = f"""# ProteinGym sealed-protein missense ensemble experiment

## Question

Can an ensemble of public zero-shot mutation-effect predictors rank held-out human missense effects more accurately than the best single predictor selected without access to the final proteins?

## Leakage control

Human proteins were deterministically partitioned before model selection into train, tuning, and sealed-final groups. The final protein labels were not used to select the method, score orientation, constituent models, weights, or the best-individual baseline.

## Data

- Human assays loaded: {diagnostics['loaded_assays']}
- Unique human proteins: {diagnostics['unique_uniprot']}
- Single substitutions: {diagnostics['mutants']:,}
- Public model-score columns with adequate coverage: {diagnostics['candidate_model_columns']}
- Split: {json.dumps(split_counts)}

## Sealed-final result

- Chosen method from tuning set: `{chosen_name}`
- Constituent model families: {', '.join(final_fit.selected)}
- Ensemble mean assay Spearman: **{final_mean:.4f}**
- Best individual selected on train+tune: `{best_base}`
- Best-individual mean assay Spearman: **{result['best_individual_final_mean_assay_spearman']:.4f}**
- Paired improvement: **{result['paired_mean_improvement']:+.4f}**
- Paired bootstrap 95% CI: **[{ci[0]:+.4f}, {ci[1]:+.4f}]**
- Final assays improved: **{100 * result['fraction_final_assays_improved']:.1f}%**
- One-sided paired Wilcoxon p: **{wilcox if wilcox is not None else 'NA'}**

## Decision

`benchmark_advance_confirmed = {result['benchmark_advance_confirmed']}`

A confirmed positive result is a reproducible benchmark advance in cross-protein DMS ranking. It is not a complete solution to all 82 million possible human missense variants and does not by itself establish clinical pathogenicity.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
