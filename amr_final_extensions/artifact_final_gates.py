#!/usr/bin/env python3
"""Final artifact-only gates for the K. pneumoniae-colistin discovery project.

Candidates are frozen from completed upstream analyses before these tests. The script:
1) maps the 61 portal-cohort replicated unitigs to the pinned reference;
2) evaluates the same frozen unitigs in the Kleborate-filtered deep-residual validation set;
3) retests frozen targeted and unitig candidates in the deposited numeric broth-dilution subset using phenotype-blind Mash-PC sensitivity models.

A positive result remains an association, not a causal resistance mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact, mannwhitneyu, norm


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--targeted-root", required=True)
    p.add_argument("--unitig-root", required=True)
    p.add_argument("--deep-root", required=True)
    p.add_argument("--phenotype-root", required=True)
    p.add_argument("--lineage-root", required=True)
    p.add_argument("--rare-root", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def find_one(root: str | Path, name: str) -> Path:
    hits = sorted(Path(root).rglob(name), key=lambda x: (len(x.parts), str(x)))
    if not hits:
        raise FileNotFoundError(f"{name} not found under {root}")
    return hits[0]


def bh(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    v = p[ok]
    order = np.argsort(v)
    ranked = v[order]
    q = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    out[np.flatnonzero(ok)[order]] = np.clip(q, 0, 1)
    return out


def safe_min(values: list[float]) -> float:
    x = np.asarray(values, dtype=float)
    return float(np.nanmin(x)) if np.isfinite(x).any() else float("nan")


def safe_max(values: list[float]) -> float:
    x = np.asarray(values, dtype=float)
    return float(np.nanmax(x)) if np.isfinite(x).any() else float("nan")


def contingency(x: np.ndarray, y: np.ndarray) -> tuple[int, int, int, int]:
    x = np.asarray(x).astype(bool)
    y = np.asarray(y).astype(bool)
    return int((x & y).sum()), int((x & ~y).sum()), int((~x & y).sum()), int((~x & ~y).sum())


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    log_or = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return math.exp(log_or), math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se)


def random_effects(x: np.ndarray, y: np.ndarray, groups: pd.Series) -> dict:
    effects = []
    gv = groups.fillna("UNKNOWN").astype(str).to_numpy()
    for group in sorted(set(gv)):
        idx = gv == group
        a, b, c, d = contingency(x[idx], y[idx])
        if min(a + c, b + d, a + b, c + d) == 0:
            continue
        aa, bb, cc, dd = map(float, (a, b, c, d))
        if min(aa, bb, cc, dd) == 0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5
        log_or = math.log((aa * dd) / (bb * cc))
        var = 1 / aa + 1 / bb + 1 / cc + 1 / dd
        effects.append((str(group), log_or, var, a, b, c, d))
    if not effects:
        return {"n_groups": 0, "details": []}
    yi = np.array([e[1] for e in effects])
    vi = np.array([e[2] for e in effects])
    wi = 1 / vi
    fixed = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * (yi - fixed) ** 2))
    df = len(effects) - 1
    cval = float(np.sum(wi) - np.sum(wi * wi) / np.sum(wi))
    tau2 = max(0.0, (q - df) / cval) if df > 0 and cval > 0 else 0.0
    wr = 1 / (vi + tau2)
    pooled = float(np.sum(wr * yi) / np.sum(wr))
    se = math.sqrt(1 / float(np.sum(wr)))
    return {
        "n_groups": len(effects),
        "odds_ratio": math.exp(pooled),
        "ci_low": math.exp(pooled - 1.96 * se),
        "ci_high": math.exp(pooled + 1.96 * se),
        "one_sided_p": float(norm.sf(pooled / se)) if se else 1.0,
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 and df > 0 else 0.0,
        "positive_groups": int(sum(e[1] > 0 for e in effects)),
        "details": [
            {"group": e[0], "log_or": e[1], "variance": e[2], "a": e[3], "b": e[4], "c": e[5], "d": e[6]}
            for e in effects
        ],
    }


def read_rtab(path: Path) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = x.index.astype(str).str.upper()
    x.columns = x.columns.astype(str)
    return x.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)


def pcoa(distance: pd.DataFrame, n_components: int) -> np.ndarray:
    d = distance.to_numpy(float)
    d = (d + d.T) / 2
    np.fill_diagonal(d, 0)
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d * d) @ j
    evals, evecs = np.linalg.eigh(b)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    keep = evals > max(float(evals[0]) * 1e-12, 1e-15)
    evals = evals[keep]
    evecs = evecs[:, keep]
    k = min(n_components, len(evals))
    if k < 1:
        raise ValueError("No positive PCoA dimension")
    z = evecs[:, :k] * np.sqrt(evals[:k])
    sd = z.std(0)
    sd[sd == 0] = 1
    return (z - z.mean(0)) / sd


def fit_null_logistic(c: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c = np.asarray(c, float)
    y = np.asarray(y, float)
    beta = np.zeros(c.shape[1])
    prevalence = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
    beta[0] = math.log(prevalence / (1 - prevalence))
    for _ in range(200):
        eta = np.clip(c @ beta, -30, 30)
        mu = 1 / (1 + np.exp(-eta))
        w = np.clip(mu * (1 - mu), 1e-8, None)
        score = c.T @ (y - mu)
        hessian = c.T @ (w[:, None] * c)
        ridge = np.eye(hessian.shape[0]) * 1e-8
        ridge[0, 0] = 0
        try:
            delta = np.linalg.solve(hessian + ridge, score)
        except np.linalg.LinAlgError:
            delta = np.linalg.pinv(hessian + ridge) @ score
        beta += delta
        if np.max(np.abs(delta)) < 1e-8:
            break
    mu = 1 / (1 + np.exp(-np.clip(c @ beta, -30, 30)))
    w = np.clip(mu * (1 - mu), 1e-8, None)
    return mu, w


def logistic_score_matrix(x: np.ndarray, y: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu, w = fit_null_logistic(c, y)
    inv = np.linalg.pinv(c.T @ (w[:, None] * c))
    x = x.astype(float)
    xres = x - c @ (inv @ (c.T @ (w[:, None] * x)))
    u = xres.T @ (y - mu)
    v = np.sum(w[:, None] * xres * xres, axis=0)
    beta = np.full(x.shape[1], np.nan)
    p = np.full(x.shape[1], np.nan)
    ok = v > 1e-12
    beta[ok] = u[ok] / v[ok]
    z = np.full(x.shape[1], np.nan)
    z[ok] = u[ok] / np.sqrt(v[ok])
    p[ok] = 2 * norm.sf(np.abs(z[ok]))
    return beta, p


def linear_score_matrix(x: np.ndarray, y: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = x.astype(float)
    y = np.asarray(y, float)
    q, _ = np.linalg.qr(c, mode="reduced")
    yr = y - q @ (q.T @ y)
    xr = x - q @ (q.T @ x)
    denom = np.sum(xr * xr, axis=0)
    beta = np.full(x.shape[1], np.nan)
    p = np.full(x.shape[1], np.nan)
    ok = denom > 1e-12
    beta[ok] = (xr[:, ok].T @ yr) / denom[ok]
    residual = yr[:, None] - xr[:, ok] * beta[ok]
    df = max(1, len(y) - c.shape[1] - 1)
    sigma2 = np.sum(residual * residual, axis=0) / df
    se = np.sqrt(sigma2 / denom[ok])
    z = beta[ok] / np.where(se > 0, se, np.nan)
    p[ok] = 2 * norm.sf(np.abs(z))
    return beta, p


def reference_unitig_context(strict: pd.DataFrame, rtab: pd.DataFrame, gbff: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = list(SeqIO.parse(str(gbff), "genbank"))
    contexts = []
    pattern_map: dict[str, str] = {}
    next_pattern = 1
    for _, row in strict.iterrows():
        seq = str(row["canonical_sequence"]).upper()
        rc = str(Seq(seq).reverse_complement())
        if seq not in rtab.index:
            continue
        bit_hash = hashlib.sha256(rtab.loc[seq].to_numpy(np.uint8).tobytes()).hexdigest()
        if bit_hash not in pattern_map:
            pattern_map[bit_hash] = f"PATTERN_{next_pattern:03d}"
            next_pattern += 1
        exact = []
        for rec in records:
            refseq = str(rec.seq).upper()
            for query, orientation in [(seq, "+"), (rc, "-")]:
                start = 0
                while True:
                    pos = refseq.find(query, start)
                    if pos < 0:
                        break
                    annotations = []
                    for feat in rec.features:
                        fs = int(feat.location.start)
                        fe = int(feat.location.end)
                        if fs < pos + len(query) and fe > pos and feat.type in {"CDS", "gene", "rRNA", "tRNA", "ncRNA"}:
                            annotations.append({
                                "type": feat.type,
                                "gene": ";".join(feat.qualifiers.get("gene", [])),
                                "product": ";".join(feat.qualifiers.get("product", [])),
                                "location": str(feat.location),
                            })
                    exact.append({
                        "record": rec.id,
                        "description": rec.description,
                        "position_1based": pos + 1,
                        "orientation": orientation,
                        "annotations": annotations,
                    })
                    start = pos + 1
        products = " | ".join(a.get("product", "") for hit in exact for a in hit["annotations"]).lower()
        descriptions = " | ".join(hit["description"] for hit in exact).lower()
        if not exact:
            category = "NO_EXACT_PINNED_REFERENCE_MATCH"
        elif "beta-lactamase" in products or "shv-" in products:
            category = "BETA_LACTAMASE_LOCUS"
        elif "transposase" in products or "insertion sequence" in products or "is26" in products:
            category = "MOBILE_ELEMENT_OR_TRANSPOSASE"
        elif "plasmid" in descriptions:
            category = "PLASMID_OTHER"
        elif any(hit["annotations"] for hit in exact):
            category = "CHROMOSOMAL_ANNOTATED_FEATURE"
        else:
            category = "REFERENCE_INTERGENIC_OR_UNANNOTATED"
        contexts.append({
            "candidate_id": row["candidate_id"],
            "canonical_sequence": seq,
            "sequence_length": len(seq),
            "occurrence_pattern": pattern_map[bit_hash],
            "reference_exact_hit_count": len(exact),
            "reference_context_category": category,
            "reference_hits_json": json.dumps(exact, ensure_ascii=False),
        })
    context = pd.DataFrame(contexts)
    patterns = context.groupby(["occurrence_pattern", "reference_context_category"], dropna=False).agg(
        n_unitigs=("candidate_id", "size"),
        candidate_ids=("candidate_id", lambda s: ";".join(s.astype(str))),
    ).reset_index()
    return context, patterns


def evaluate_frozen_unitigs_in_deep_residual(strict: pd.DataFrame, rtab: pd.DataFrame, deep_manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    validation = deep_manifest[
        deep_manifest["split"].astype(str).eq("validation")
        & deep_manifest["assembly_ID"].astype(str).isin(rtab.columns)
    ].copy()
    y = validation["phenotype"].astype(str).eq("R").to_numpy()
    rows = []
    details = {}
    for _, candidate in strict.iterrows():
        seq = str(candidate["canonical_sequence"]).upper()
        if seq not in rtab.index:
            continue
        x = validation["assembly_ID"].map(rtab.loc[seq]).fillna(0).astype(int).to_numpy()
        a, b, c, d = contingency(x, y)
        odds_ratio, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
        p = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        source = random_effects(x, y, validation.get("BioProject", validation.get("source_group", pd.Series(["UNKNOWN"] * len(validation)))))
        country = random_effects(x, y, validation.get("ISO_country_code", pd.Series(["UNKNOWN"] * len(validation))))
        details[str(candidate["candidate_id"])] = {"source": source, "country": country}
        rows.append({
            "candidate_id": candidate["candidate_id"],
            "canonical_sequence": seq,
            "validation_R_present": a,
            "validation_S_present": b,
            "validation_R_absent": c,
            "validation_S_absent": d,
            "validation_or": odds_ratio,
            "validation_ci_low": ci_low,
            "validation_ci_high": ci_high,
            "validation_p": p,
            "source_n": source.get("n_groups", 0),
            "source_or": source.get("odds_ratio"),
            "source_ci_low": source.get("ci_low"),
            "source_p": source.get("one_sided_p"),
            "country_n": country.get("n_groups", 0),
            "country_or": country.get("odds_ratio"),
            "country_ci_low": country.get("ci_low"),
            "country_p": country.get("one_sided_p"),
        })
    evidence = pd.DataFrame(rows)
    evidence["validation_q"] = bh(evidence["validation_p"])
    evidence["heldout_replication"] = (
        (evidence.validation_R_present + evidence.validation_S_present >= 5)
        & (evidence.validation_or > 1)
        & (evidence.validation_ci_low > 1)
        & (evidence.validation_q <= 0.05)
    )
    evidence["source_replication"] = (
        (evidence.source_n >= 3)
        & (pd.to_numeric(evidence.source_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(evidence.source_p, errors="coerce") <= 0.05)
    )
    evidence["country_replication"] = (
        (evidence.country_n >= 3)
        & (pd.to_numeric(evidence.country_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(evidence.country_p, errors="coerce") <= 0.05)
    )
    evidence["strict_deep_residual_replication"] = evidence.heldout_replication & evidence.source_replication & evidence.country_replication
    evidence = evidence.sort_values(["strict_deep_residual_replication", "validation_q"], ascending=[False, True])
    summary = {
        "n_deep_validation": len(validation),
        "deep_validation_counts": validation.phenotype.value_counts().to_dict(),
        "n_tested": len(evidence),
        "n_heldout_replicated": int(evidence.heldout_replication.sum()),
        "n_source_replicated": int(evidence.source_replication.sum()),
        "n_country_replicated": int(evidence.country_replication.sum()),
        "n_strict_deep_residual_replicated": int(evidence.strict_deep_residual_replication.sum()),
        "strict_candidate_ids": evidence.loc[evidence.strict_deep_residual_replication, "candidate_id"].astype(str).tolist(),
        "boundary": "The original 61 portal-cohort unitigs were frozen before this Kleborate-filtered validation test. Statistical persistence would not establish causality or novelty.",
    }
    return evidence, {"summary": summary, "details": details}


def aggregate_numeric_broth(rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    x = rows.copy()
    method = x["laboratory_typing_method"].fillna("").astype(str)
    units = x["measurement_units"].fillna("").astype(str)
    x = x[
        method.str.contains("broth", case=False, regex=False)
        & units.str.contains("mg", case=False, regex=False)
    ].copy()
    x["measurement_numeric"] = pd.to_numeric(x["measurement_numeric"], errors="coerce")
    x = x[x.measurement_numeric.gt(0) & x.raw_sir.astype(str).isin(["R", "S"])].copy()
    records = []
    conflicts = []
    for accession, group in x.groupby("assembly_ID"):
        labels = sorted(set(group.raw_sir.astype(str)))
        if len(labels) != 1:
            conflicts.append({"assembly_ID": accession, "reason": "conflicting_raw_sir", "values": ";".join(labels)})
            continue
        median_mic = float(group.measurement_numeric.median())
        records.append({
            "assembly_ID": str(accession),
            "raw_sir": labels[0],
            "reported_mic_mg_l": median_mic,
            "log2_reported_mic": float(np.log2(median_mic)),
            "measurement_signs": ";".join(sorted(set(group.measurement_sign.fillna("").astype(str)))),
            "platforms": ";".join(sorted(set(group.platform.fillna("").astype(str)) - {""})),
            "methods": ";".join(sorted(set(group.laboratory_typing_method.fillna("").astype(str)) - {""})),
            "n_rows": len(group),
            "censored_measurement": bool(~group.measurement_sign.astype(str).isin(["=", "=="]).all()),
            "sensititre": bool(
                group.platform.fillna("").astype(str).str.contains("Sensititre", case=False).any()
                or group.laboratory_typing_method.fillna("").astype(str).str.contains("Sensititre", case=False).any()
            ),
        })
    return pd.DataFrame(records), pd.DataFrame(conflicts)


def panel_tests(
    panel_name: str,
    candidates: pd.DataFrame,
    matrix: pd.DataFrame,
    manifest: pd.DataFrame,
    broth: pd.DataFrame,
    distance: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    merged = manifest.merge(broth, on="assembly_ID", how="inner", validate="one_to_one")
    candidate_rows = []
    for _, candidate in candidates.iterrows():
        key = str(candidate["matrix_key"]).upper()
        if key in matrix.index:
            candidate_rows.append((str(candidate["candidate_id"]), key))
    if not candidate_rows:
        return pd.DataFrame(), {"panel": panel_name, "n_samples": len(merged), "status": "NO_TESTABLE_CANDIDATE"}
    keys = [x[1] for x in candidate_rows]
    evidence_rows = []
    for split in ["discovery", "validation"]:
        subset = merged[merged["split"].astype(str).eq(split)].copy()
        subset = subset[subset.assembly_ID.isin(matrix.columns) & subset.assembly_ID.isin(distance.index)].copy()
        if subset.empty:
            continue
        sample_ids = subset.assembly_ID.astype(str).tolist()
        x = matrix.loc[keys, sample_ids].T.to_numpy(np.uint8)
        y = subset.raw_sir.eq("R").to_numpy(np.uint8)
        mic = subset.log2_reported_mic.to_numpy(float)
        dimensions = [d for d in ([5, 10] if split == "discovery" else [3, 5]) if d < len(subset) - 4]
        if not dimensions:
            dimensions = [max(1, min(2, len(subset) - 4))]
        distance_subset = distance.loc[sample_ids, sample_ids]
        binary_models = []
        mic_models = []
        for dimension in dimensions:
            pcs = pcoa(distance_subset, dimension)
            covariates = np.column_stack([np.ones(len(subset)), pcs])
            binary_beta, binary_p = logistic_score_matrix(x, y, covariates)
            mic_beta, mic_p = linear_score_matrix(x, mic, covariates)
            binary_models.append((dimension, binary_beta, binary_p, bh(binary_p)))
            mic_models.append((dimension, mic_beta, mic_p, bh(mic_p)))
        for j, (candidate_id, key) in enumerate(candidate_rows):
            present = x[:, j]
            a, b, c, d = contingency(present, y)
            odds_ratio, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
            fisher_p = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
            carrier = mic[present == 1]
            noncarrier = mic[present == 0]
            if len(carrier) >= 3 and len(noncarrier) >= 3:
                rank_p = float(mannwhitneyu(carrier, noncarrier, alternative="greater", method="asymptotic").pvalue)
                carrier_median = float(np.median(carrier))
                noncarrier_median = float(np.median(noncarrier))
            else:
                rank_p = float("nan")
                carrier_median = float(np.median(carrier)) if len(carrier) else float("nan")
                noncarrier_median = float(np.median(noncarrier)) if len(noncarrier) else float("nan")
            row = {
                "panel": panel_name,
                "candidate_id": candidate_id,
                "matrix_key": key,
                "split": split,
                "n": len(subset),
                "R": int(y.sum()),
                "S": int((1 - y).sum()),
                "R_present": a,
                "S_present": b,
                "R_absent": c,
                "S_absent": d,
                "binary_or": odds_ratio,
                "binary_ci_low": ci_low,
                "binary_ci_high": ci_high,
                "binary_p": fisher_p,
                "carrier_n": len(carrier),
                "noncarrier_n": len(noncarrier),
                "carrier_median_log2_mic": carrier_median,
                "noncarrier_median_log2_mic": noncarrier_median,
                "mic_rank_p": rank_p,
                "pc_dims": ";".join(map(str, dimensions)),
            }
            for dimension, beta, pvalue, qvalue in binary_models:
                row[f"binary_score_beta_pc{dimension}"] = beta[j]
                row[f"binary_score_p_pc{dimension}"] = pvalue[j]
                row[f"binary_score_q_pc{dimension}"] = qvalue[j]
            for dimension, beta, pvalue, qvalue in mic_models:
                row[f"mic_score_beta_pc{dimension}"] = beta[j]
                row[f"mic_score_p_pc{dimension}"] = pvalue[j]
                row[f"mic_score_q_pc{dimension}"] = qvalue[j]
            binary_betas = [beta[j] for _, beta, _, _ in binary_models]
            binary_qs = [qvalue[j] for _, _, _, qvalue in binary_models]
            mic_betas = [beta[j] for _, beta, _, _ in mic_models]
            mic_qs = [qvalue[j] for _, _, _, qvalue in mic_models]
            row["binary_structure_beta_min"] = safe_min(binary_betas)
            row["binary_structure_q_max"] = safe_max(binary_qs)
            row["mic_structure_beta_min"] = safe_min(mic_betas)
            row["mic_structure_q_max"] = safe_max(mic_qs)
            row["binary_structure_stable"] = bool(
                np.isfinite(row["binary_structure_beta_min"])
                and np.isfinite(row["binary_structure_q_max"])
                and row["binary_structure_beta_min"] > 0
                and row["binary_structure_q_max"] <= 0.05
            )
            row["mic_structure_stable"] = bool(
                np.isfinite(row["mic_structure_beta_min"])
                and np.isfinite(row["mic_structure_q_max"])
                and row["mic_structure_beta_min"] > 0
                and row["mic_structure_q_max"] <= 0.05
            )
            evidence_rows.append(row)
    evidence = pd.DataFrame(evidence_rows)
    evidence["binary_q"] = np.nan
    evidence["mic_rank_q"] = np.nan
    for split, idx in evidence.groupby("split").groups.items():
        evidence.loc[idx, "binary_q"] = bh(evidence.loc[idx, "binary_p"])
        evidence.loc[idx, "mic_rank_q"] = bh(evidence.loc[idx, "mic_rank_p"])
    evidence["refined_binary_replication"] = (
        evidence.split.eq("validation")
        & (evidence.R_present + evidence.S_present >= 5)
        & (evidence.binary_or > 1)
        & (evidence.binary_ci_low > 1)
        & (evidence.binary_q <= 0.05)
        & evidence.binary_structure_stable
    )
    evidence["refined_mic_direction"] = (
        evidence.split.eq("validation")
        & (evidence.carrier_median_log2_mic > evidence.noncarrier_median_log2_mic)
        & (evidence.mic_rank_q <= 0.05)
        & evidence.mic_structure_stable
    )
    evidence["strict_numeric_broth_replication"] = evidence.refined_binary_replication & evidence.refined_mic_direction
    summary = {
        "panel": panel_name,
        "n_numeric_broth_samples": len(merged),
        "sample_counts": {
            f"{key[0]}|{key[1]}": int(value)
            for key, value in merged.groupby(["split", "raw_sir"]).size().to_dict().items()
        },
        "n_bioprojects": int(merged.BioProject.nunique()) if "BioProject" in merged else None,
        "n_candidates": int(evidence.candidate_id.nunique()),
        "n_validation_binary_replicated": int(evidence.refined_binary_replication.sum()),
        "n_validation_mic_direction": int(evidence.refined_mic_direction.sum()),
        "n_strict_numeric_broth_replicated": int(evidence.strict_numeric_broth_replication.sum()),
        "strict_candidate_ids": evidence.loc[evidence.strict_numeric_broth_replication, "candidate_id"].astype(str).tolist(),
        "boundary": "Numeric broth-dilution is a deposited-method sensitivity subset, not independently verified reference BMD. MIC signs include censoring; rank/linear tests use the reported numeric boundary and are sensitivity analyses. Population structure is controlled with phenotype-blind Mash PCoA sensitivity models.",
    }
    return evidence, summary


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    targeted_evidence = pd.read_csv(find_one(a.targeted_root, "ALL_MARKER_EVIDENCE.csv"))
    manifest_candidates = sorted(Path(a.targeted_root).rglob("gwas_sample_manifest.csv"))
    targeted_manifest = pd.read_csv(next((p for p in manifest_candidates if "gwas_inputs_bioproject" in str(p)), manifest_candidates[0]))
    bioproject_input = Path(a.targeted_root) / "gwas_inputs_bioproject"
    targeted_rtab = read_rtab(find_one(bioproject_input if bioproject_input.exists() else a.targeted_root, "all_variants.Rtab"))
    unitig_strict = pd.read_csv(find_one(a.unitig_root, "STRICT_STATISTICALLY_REPLICATED_UNITIGS.csv"))
    unitig_manifest = pd.read_csv(find_one(a.unitig_root, "gwas_sample_manifest.csv"))
    unitig_rtab = read_rtab(find_one(a.unitig_root, "all.rtab"))
    deep_manifest = pd.read_csv(find_one(a.deep_root, "gwas_sample_manifest.csv"))
    lineage_frozen = pd.read_csv(find_one(a.lineage_root, "FROZEN_LINEAGE_STRATIFIED_CANDIDATES.csv"))
    lineage_strict = pd.read_csv(find_one(a.lineage_root, "STRICT_LINEAGE_STRATIFIED_REPLICATES.csv"))
    rare = pd.read_csv(find_one(a.rare_root, "RARE_BURDEN_DISCOVERY_VALIDATION_EVIDENCE.csv"))
    phenotype_rows = pd.read_csv(find_one(a.phenotype_root, "all_target_phenotype_rows.csv"))
    gbff = find_one(a.targeted_root, "reference.gbff")
    distance_hits = sorted(Path(a.targeted_root).rglob("all_mash.tsv"))
    distance_path = next((p for p in distance_hits if "structure_bioproject" in str(p)), distance_hits[0])
    distance = pd.read_csv(distance_path, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)

    targeted_features = sorted(
        set(targeted_evidence.loc[targeted_evidence.discovery_stable.eq(True), "feature"].astype(str))
        | set(lineage_frozen.feature.astype(str))
    )
    targeted_panel = pd.DataFrame({
        "candidate_id": ["TARGET_" + feature for feature in targeted_features],
        "matrix_key": targeted_features,
        "candidate_type": "targeted_variant",
    })
    unitig_panel = unitig_strict[["candidate_id", "canonical_sequence"]].rename(columns={"canonical_sequence": "matrix_key"}).copy()
    unitig_panel["candidate_type"] = "unitig"
    frozen = pd.concat([targeted_panel, unitig_panel], ignore_index=True)
    frozen.to_csv(out / "FROZEN_CANDIDATE_PANEL.csv", index=False)
    lineage_strict.to_csv(out / "LINEAGE_STRICT_INPUT.csv", index=False)
    rare.to_csv(out / "RARE_BURDEN_INPUT.csv", index=False)

    context, pattern_summary = reference_unitig_context(unitig_strict, unitig_rtab, gbff)
    context.to_csv(out / "UNITIG_PINNED_REFERENCE_CONTEXT.csv", index=False)
    pattern_summary.to_csv(out / "UNITIG_OCCURRENCE_PATTERN_SUMMARY.csv", index=False)

    deep_evidence, deep_object = evaluate_frozen_unitigs_in_deep_residual(unitig_strict, unitig_rtab, deep_manifest)
    deep_evidence.to_csv(out / "ORIGINAL_UNITIGS_IN_DEEP_RESIDUAL_VALIDATION.csv", index=False)
    (out / "ORIGINAL_UNITIGS_DEEP_RESIDUAL_DETAILS.json").write_text(
        json.dumps(deep_object, indent=2, ensure_ascii=False, default=str) + "\n"
    )

    broth, conflicts = aggregate_numeric_broth(phenotype_rows)
    broth.to_csv(out / "NUMERIC_BROTH_ASSEMBLY_COHORT.csv", index=False)
    conflicts.to_csv(out / "NUMERIC_BROTH_CONFLICTS.csv", index=False)
    targeted_tests, targeted_summary = panel_tests(
        "targeted", targeted_panel, targeted_rtab, targeted_manifest, broth, distance
    )
    unitig_tests, unitig_summary = panel_tests(
        "unitig", unitig_panel, unitig_rtab, unitig_manifest, broth, distance
    )
    targeted_tests.to_csv(out / "TARGETED_NUMERIC_BROTH_TESTS.csv", index=False)
    unitig_tests.to_csv(out / "UNITIG_NUMERIC_BROTH_TESTS.csv", index=False)

    sensititre = broth[broth.sensititre].copy()
    sensititre_summary = {
        "n_all_target_sensititre_assemblies": len(sensititre),
        "label_counts": sensititre.raw_sir.value_counts().to_dict(),
        "targeted_manifest_overlap": int(sensititre.assembly_ID.isin(targeted_manifest.assembly_ID).sum()),
        "unitig_manifest_overlap": int(sensititre.assembly_ID.isin(unitig_manifest.assembly_ID).sum()),
        "boundary": "The Sensititre overlap is descriptive because it is too small for a reliable multivariable discovery/validation claim.",
    }

    artifact_promoted = (
        int(deep_evidence.strict_deep_residual_replication.sum())
        + int(targeted_tests.strict_numeric_broth_replication.sum() if len(targeted_tests) else 0)
        + int(unitig_tests.strict_numeric_broth_replication.sum() if len(unitig_tests) else 0)
    )
    summary = {
        "frozen_targeted_features": len(targeted_panel),
        "frozen_unitigs": len(unitig_panel),
        "frozen_total": len(frozen),
        "unitig_reference_context_counts": context.reference_context_category.value_counts().to_dict(),
        "unitig_unique_occurrence_patterns": int(context.occurrence_pattern.nunique()),
        "deep_residual_cross_method": deep_object["summary"],
        "numeric_broth_targeted": targeted_summary,
        "numeric_broth_unitig": unitig_summary,
        "sensititre": sensititre_summary,
        "n_artifact_only_promoted_candidates": artifact_promoted,
        "status": "CANDIDATE_REQUIRES_EXTERNAL_AND_STRUCTURAL_AUDIT" if artifact_promoted else "NO_NOVEL_CANDIDATE_SURVIVED_ARTIFACT_ONLY_FINAL_GATES",
        "claim_boundary": "No candidate may be called novel or causal from these analyses. Known mgrB disruption is an established positive-control mechanism. Structural context, independent external phenotype-linked genomes, current database/literature audit, and laboratory validation remain separate gates.",
    }
    (out / "ARTIFACT_FINAL_GATES_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str) + "\n"
    )
    report = [
        "# K. pneumoniae-colistin artifact-only final gates",
        "",
        f"- Frozen targeted features: **{len(targeted_panel)}**",
        f"- Frozen unitigs: **{len(unitig_panel)}**",
        f"- Unique unitig occurrence patterns: **{context.occurrence_pattern.nunique()}**",
        f"- Original unitigs surviving complete deep-residual validation: **{int(deep_evidence.strict_deep_residual_replication.sum())}**",
        f"- Targeted candidates surviving population-structure-adjusted numeric-broth binary+MIC gates: **{int(targeted_tests.strict_numeric_broth_replication.sum()) if len(targeted_tests) else 0}**",
        f"- Unitigs surviving population-structure-adjusted numeric-broth binary+MIC gates: **{int(unitig_tests.strict_numeric_broth_replication.sum()) if len(unitig_tests) else 0}**",
        "",
        f"**Status: {summary['status']}**",
        "",
        summary["claim_boundary"],
    ]
    (out / "ARTIFACT_FINAL_GATES_REPORT.md").write_text("\n".join(report) + "\n")
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
