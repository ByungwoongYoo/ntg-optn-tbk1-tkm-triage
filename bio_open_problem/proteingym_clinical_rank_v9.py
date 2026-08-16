#!/usr/bin/env python3
"""Leakage-controlled, exact-sequence-grouped ProteinGym clinical rank ensemble."""
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
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

SEED = 20260817
META = {
    "Unnamed: 0", "protein", "protein_sequence", "mutant", "mutated_sequence",
    "DMS_bin_score", "DMS_score_bin", "label"
}


def label01(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().map({
        "pathogenic": 1, "likely pathogenic": 1,
        "benign": 0, "likely benign": 0,
    })


def rank01(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").rank(method="average", pct=True)


def stable_bucket(text: str) -> int:
    return int(hashlib.sha256(text.encode()).hexdigest()[:12], 16) % 10


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(str(seq).strip().upper().encode()).hexdigest()


def safe_auc(y: pd.Series | np.ndarray, p: pd.Series | np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    m = np.isfinite(y) & np.isfinite(p)
    if m.sum() < 4 or len(np.unique(y[m])) != 2:
        return float("nan")
    return float(roc_auc_score(y[m].astype(int), p[m]))


def bootstrap_ci(v: np.ndarray, n: int = 20000) -> list[float | None]:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return [None, None]
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(v), size=(n, len(v)))
    means = v[idx].mean(axis=1)
    return [float(np.quantile(means, .025)), float(np.quantile(means, .975))]


def family(name: str) -> str:
    n = str(name).lower()
    rules = [
        (r"polyphen", "PolyPhen"), (r"sift", "SIFT"), (r"varity", "VARITY"),
        (r"bayesdel", "BayesDel"), (r"trancepteve", "TranceptEVE"),
        (r"^eve$", "EVE"), (r"poet", "PoET"), (r"gemme", "GEMME"),
        (r"esm1b", "ESM1b"), (r"mutationtaster", "MutationTaster"),
        (r"mutationassessor", "MutationAssessor"), (r"clinpred", "ClinPred"),
        (r"dann", "DANN"), (r"cadd", "CADD"), (r"fathmm", "FATHMM"),
        (r"list", "LIST"), (r"primate", "PrimateAI"), (r"mutpred", "MutPred"),
        (r"metarnn", "MetaRNN"), (r"revel", "REVEL"), (r"vest", "VEST"),
        (r"provean", "PROVEAN"), (r"deogen", "DEOGEN"), (r"gmvp", "gMVP"),
        (r"mpc", "MPC"), (r"^lrt", "LRT"),
    ]
    for pat, fam in rules:
        if re.search(pat, n):
            return fam
    return re.sub(r"[^a-z0-9]+", "", n)


def load_scores(path: Path) -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    frames: list[pd.DataFrame] = []
    header_counts: dict[str, int] = {}
    failures: list[dict[str, str]] = []
    with zipfile.ZipFile(path) as z:
        members = [n for n in z.namelist() if n.lower().endswith(".csv")]
        for i, member in enumerate(members, 1):
            try:
                with z.open(member) as f:
                    d = pd.read_csv(f, low_memory=False)
                required = {"mutant", "DMS_bin_score", "protein_sequence"}
                if not required.issubset(d.columns):
                    failures.append({"file": member, "reason": "missing required fields"})
                    continue
                model_cols = [c for c in d.columns if c not in META]
                for c in model_cols:
                    header_counts[c] = header_counts.get(c, 0) + 1
                d["label"] = label01(d["DMS_bin_score"])
                d = d[d["label"].isin([0, 1])].drop_duplicates("mutant")
                if len(d) < 4 or d["label"].nunique() != 2:
                    continue
                d["protein_file"] = Path(member).name
                d["sequence_hash"] = sequence_hash(d["protein_sequence"].iloc[0])
                frames.append(d[["protein_file", "sequence_hash", "mutant", "label"] + model_cols])
            except Exception as exc:
                failures.append({"file": member, "reason": repr(exc)})
            if i % 500 == 0:
                print(f"load {i}/{len(members)} accepted={len(frames)}", flush=True)
    if not frames:
        raise RuntimeError("No clinical score files loaded")
    min_files = math.ceil(.95 * len(frames))
    models = sorted([c for c, n in header_counts.items() if n >= min_files])
    data = pd.concat(frames, ignore_index=True)
    for c in models:
        if c not in data.columns:
            data[c] = np.nan
        data[c] = data.groupby("protein_file")[c].transform(rank01)
    data = data[["protein_file", "sequence_hash", "mutant", "label"] + models]
    return data, models, {
        "accepted_protein_files": len(frames), "variants": len(data),
        "sequence_groups": int(data["sequence_hash"].nunique()),
        "models": models, "failures": failures,
    }


def protein_auc(df: pd.DataFrame, prediction: np.ndarray, name: str) -> pd.DataFrame:
    t = df[["protein_file", "sequence_hash", "label"]].copy()
    t[name] = prediction
    rows = []
    for p, g in t.groupby("protein_file", sort=False):
        rows.append({
            "protein_file": p, "sequence_hash": g["sequence_hash"].iloc[0],
            "n_variants": len(g), "n_pathogenic": int(g["label"].sum()),
            name: safe_auc(g["label"], g[name]),
        })
    return pd.DataFrame(rows)


def orientation_performance(df: pd.DataFrame, models: list[str]) -> tuple[dict[str, int], dict[str, float]]:
    signs: dict[str, int] = {}
    perf: dict[str, float] = {}
    for m in models:
        vals = []
        for _, g in df.groupby("protein_file", sort=False):
            vals.append(safe_auc(g["label"], g[m]))
        mean = float(np.nanmean(vals))
        signs[m] = 1 if not np.isfinite(mean) or mean >= .5 else -1
        perf[m] = mean if signs[m] > 0 else 1.0 - mean
    return signs, perf


def select_diverse(perf: dict[str, float], k: int) -> list[str]:
    best: dict[str, tuple[str, float]] = {}
    for m, p in perf.items():
        f = family(m)
        if f not in best or p > best[f][1]:
            best[f] = (m, p)
    ranked = sorted(best.values(), key=lambda x: (-np.nan_to_num(x[1], nan=-9), x[0]))
    return [m for m, _ in ranked[:min(k, len(ranked))]]


def oriented_matrix(df: pd.DataFrame, models: list[str], signs: dict[str, int]) -> np.ndarray:
    x = df[models].to_numpy(dtype=float, copy=True)
    for j, m in enumerate(models):
        if signs[m] < 0:
            x[:, j] = 1.0 - x[:, j]
    return x


@dataclass(frozen=True)
class Spec:
    name: str
    kind: str
    k: int
    power: float = 1.0


SPECS = [
    Spec("mean3", "mean", 3), Spec("mean5", "mean", 5),
    Spec("mean10", "mean", 10), Spec("mean15", "mean", 15),
    Spec("mean20", "mean", 20), Spec("mean_all", "mean", 99),
    Spec("median3", "median", 3), Spec("median5", "median", 5),
    Spec("median10", "median", 10), Spec("median20", "median", 20),
    Spec("median_all", "median", 99),
    Spec("weighted5_p1", "weighted", 5, 1), Spec("weighted10_p1", "weighted", 10, 1),
    Spec("weighted20_p1", "weighted", 20, 1), Spec("weighted_all_p1", "weighted", 99, 1),
    Spec("weighted5_p2", "weighted", 5, 2), Spec("weighted10_p2", "weighted", 10, 2),
    Spec("weighted20_p2", "weighted", 20, 2), Spec("weighted_all_p2", "weighted", 99, 2),
    Spec("weighted10_p4", "weighted", 10, 4), Spec("weighted20_p4", "weighted", 20, 4),
]


def predict(df: pd.DataFrame, spec: Spec, models: list[str], signs: dict[str, int], perf: dict[str, float]) -> np.ndarray:
    selected = select_diverse(perf, spec.k)
    x = oriented_matrix(df, selected, signs)
    if spec.kind == "mean":
        return np.nanmean(x, axis=1)
    if spec.kind == "median":
        return np.nanmedian(x, axis=1)
    w = np.array([max(perf[m] - .5, .001) ** spec.power for m in selected], dtype=float)
    w /= w.sum()
    ok = np.isfinite(x)
    ww = ok * w.reshape(1, -1)
    den = ww.sum(axis=1)
    out = np.full(len(df), np.nan)
    valid = den > 0
    out[valid] = np.nansum(x[valid] * w.reshape(1, -1), axis=1) / den[valid]
    return out


def compare(tab: pd.DataFrame, a: str, b: str) -> dict[str, Any]:
    d = (tab[a] - tab[b]).to_numpy(float)
    d = d[np.isfinite(d)]
    try:
        p = float(wilcoxon(d, alternative="greater").pvalue) if len(d) >= 5 and np.any(d != 0) else None
    except Exception:
        p = None
    return {
        "n_proteins": int(len(d)), "ensemble_mean_auc": float(tab[a].mean()),
        "comparator_mean_auc": float(tab[b].mean()), "mean_gain": float(np.mean(d)),
        "gain_95ci": bootstrap_ci(d), "fraction_improved": float(np.mean(d > 0)),
        "wilcoxon_one_sided_p": p,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical-score-zip", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    data, models, diagnostics = load_scores(Path(args.clinical_score_zip))
    data["split"] = data["sequence_hash"].map(
        lambda h: "train" if stable_bucket(h) <= 5 else ("tune" if stable_bucket(h) <= 7 else "final")
    )
    split = data[["protein_file", "sequence_hash", "split"]].drop_duplicates().groupby("split").agg(
        protein_files=("protein_file", "nunique"), sequence_groups=("sequence_hash", "nunique")
    ).to_dict("index")
    train = data[data["split"].eq("train")].copy()
    tune = data[data["split"].eq("tune")].copy()
    final = data[data["split"].eq("final")].copy()

    signs_train, perf_train = orientation_performance(train, models)
    tuning_rows = []
    for spec in SPECS:
        p = predict(tune, spec, models, signs_train, perf_train)
        tab = protein_auc(tune, p, spec.name)
        tuning_rows.append({
            "method": spec.name, "kind": spec.kind, "k": spec.k,
            "power": spec.power, "mean_tune_auc": float(tab[spec.name].mean()),
            "median_tune_auc": float(tab[spec.name].median()),
            "selected_models": ";".join(select_diverse(perf_train, spec.k)),
        })
        print(spec.name, tuning_rows[-1]["mean_tune_auc"], flush=True)
    tuning = pd.DataFrame(tuning_rows).sort_values(["mean_tune_auc", "method"], ascending=[False, True])
    chosen_name = str(tuning.iloc[0]["method"])
    chosen = next(s for s in SPECS if s.name == chosen_name)

    train_tune = data[data["split"].isin(["train", "tune"])].copy()
    signs, perf = orientation_performance(train_tune, models)
    final_pred = predict(final, chosen, models, signs, perf)
    ens = protein_auc(final, final_pred, "ensemble")

    best_selected = max(perf, key=lambda m: np.nan_to_num(perf[m], nan=-9))
    bp = final[best_selected].to_numpy(float, copy=True)
    if signs[best_selected] < 0:
        bp = 1.0 - bp
    base = protein_auc(final, bp, "best_selected")
    comp = ens.merge(base[["protein_file", "best_selected"]], on="protein_file")

    individual_rows = []
    individual_tabs: dict[str, pd.DataFrame] = {}
    for m in models:
        q = final[m].to_numpy(float, copy=True)
        if signs[m] < 0:
            q = 1.0 - q
        t = protein_auc(final, q, "auc")
        individual_tabs[m] = t
        individual_rows.append({"model": m, "mean_final_auc": float(t["auc"].mean()), "n_proteins": int(t["auc"].notna().sum())})
    individuals = pd.DataFrame(individual_rows).sort_values("mean_final_auc", ascending=False)
    oracle = str(individuals.iloc[0]["model"])
    comp = comp.merge(individual_tabs[oracle][["protein_file", "auc"]].rename(columns={"auc": "oracle_best"}), on="protein_file")

    vs_selected = compare(comp, "ensemble", "best_selected")
    vs_oracle = compare(comp, "ensemble", "oracle_best")
    result = {
        "question": "Can an exact-sequence-grouped, protein-held-out rank ensemble improve ProteinGym clinical missense classification?",
        "diagnostics": diagnostics, "split": split,
        "chosen_method": chosen_name, "chosen_models": select_diverse(perf, chosen.k),
        "sealed_final_protein_files": int(final["protein_file"].nunique()),
        "sealed_final_sequence_groups": int(final["sequence_hash"].nunique()),
        "sealed_final_variants": int(len(final)),
        "best_individual_selected_without_final_labels": best_selected,
        "posthoc_best_final_individual": oracle,
        "vs_selected_best": vs_selected, "vs_posthoc_oracle": vs_oracle,
        "benchmark_advance_confirmed": bool(vs_selected["gain_95ci"][0] is not None and vs_selected["gain_95ci"][0] > 0),
        "beats_posthoc_oracle_confirmed": bool(vs_oracle["gain_95ci"][0] is not None and vs_oracle["gain_95ci"][0] > 0),
        "open_problem_fully_solved": False,
    }
    tuning.to_csv(out / "tuning_results.csv", index=False)
    individuals.to_csv(out / "final_individual_models.csv", index=False)
    comp.to_csv(out / "sealed_final_protein_auc.csv", index=False)
    pd.DataFrame({
        "protein_file": final["protein_file"], "sequence_hash": final["sequence_hash"],
        "mutant": final["mutant"], "label": final["label"], "ensemble_prediction": final_pred,
    }).to_csv(out / "sealed_final_predictions.csv.gz", index=False, compression="gzip")
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "REPORT.md").write_text(f"""# ProteinGym clinical rank ensemble v9

Exact protein sequences were assigned before method selection to train, tuning, and sealed-final sets. The sealed-final labels were not used to choose the ensemble form, constituent models, score directions, weights, or comparator.

- Protein files: {diagnostics['accepted_protein_files']}
- Variants: {diagnostics['variants']:,}
- Official predictor columns: {len(models)}
- Split: {json.dumps(split)}
- Chosen method: `{chosen_name}`
- Sealed-final sequence groups: {result['sealed_final_sequence_groups']}
- Sealed-final variants: {result['sealed_final_variants']:,}
- Ensemble mean protein AUC: {vs_selected['ensemble_mean_auc']:.4f}
- Best individual selected on train+tune: `{best_selected}` = {vs_selected['comparator_mean_auc']:.4f}
- Gain: {vs_selected['mean_gain']:+.4f}; 95% CI {vs_selected['gain_95ci']}
- Post-hoc best final individual: `{oracle}`
- Gain versus post-hoc oracle: {vs_oracle['mean_gain']:+.4f}; 95% CI {vs_oracle['gain_95ci']}
- Benchmark advance confirmed: {result['benchmark_advance_confirmed']}

This is a fixed benchmark result, not a complete solution to all human missense effects or clinical pathogenicity.
""", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
