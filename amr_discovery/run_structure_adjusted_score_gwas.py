#!/usr/bin/env python3
"""Population-structure-adjusted binary-feature GWAS with untouched validation.

This implements a conditional logistic score test under a null model containing
classical-MDS coordinates derived from a genomic distance matrix. It is intended
as an independent cross-check of pyseer, not as a substitute for biological
validation. Feature selection uses discovery samples only; validation samples are
not examined until the discovery list is frozen.
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
from scipy.special import expit
from scipy.stats import fisher_exact, norm


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--rtab", required=True)
    p.add_argument("--distance", required=True)
    p.add_argument("--feature-meta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--pc-dims", default="10,20,30")
    p.add_argument("--min-af", type=float, default=0.01)
    p.add_argument("--max-af", type=float, default=0.99)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5)
    p.add_argument("--min-validation-odds-ratio", type=float, default=2.0)
    return p.parse_args()


def bh(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    v = p[ok]
    order = np.argsort(v)
    ranked = v[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out[np.flatnonzero(ok)[order]] = q
    return out


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


def pcoa(distance: pd.DataFrame, n_components: int) -> np.ndarray:
    d = distance.to_numpy(dtype=float)
    if d.shape[0] != d.shape[1]:
        raise ValueError("Distance matrix is not square")
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d * d) @ j
    evals, evecs = np.linalg.eigh(b)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    positive = evals > max(float(evals[0]) * 1e-12, 1e-15)
    evals = evals[positive]
    evecs = evecs[:, positive]
    k = min(n_components, len(evals))
    if k < 1:
        raise ValueError("No positive PCoA eigenvalue")
    coords = evecs[:, :k] * np.sqrt(evals[:k])
    means = coords.mean(axis=0)
    sds = coords.std(axis=0, ddof=0)
    sds[sds == 0] = 1.0
    return (coords - means) / sds


def fit_null_logistic(covariates: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c = np.asarray(covariates, dtype=float)
    y = np.asarray(y, dtype=float)
    beta = np.zeros(c.shape[1], dtype=float)
    prevalence = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
    beta[0] = math.log(prevalence / (1 - prevalence))
    for _ in range(200):
        eta = np.clip(c @ beta, -30, 30)
        mu = expit(eta)
        w = np.clip(mu * (1 - mu), 1e-8, None)
        score = c.T @ (y - mu)
        hessian = c.T @ (w[:, None] * c)
        ridge = np.eye(hessian.shape[0]) * 1e-8
        ridge[0, 0] = 0.0
        try:
            delta = np.linalg.solve(hessian + ridge, score)
        except np.linalg.LinAlgError:
            delta = np.linalg.pinv(hessian + ridge) @ score
        beta += delta
        if float(np.max(np.abs(delta))) < 1e-8:
            break
    mu = expit(np.clip(c @ beta, -30, 30))
    w = np.clip(mu * (1 - mu), 1e-8, None)
    return beta, mu, w


def score_test_matrix(x: np.ndarray, y: np.ndarray, c: np.ndarray, batch_size: int = 512) -> pd.DataFrame:
    _, mu, w = fit_null_logistic(c, y)
    ctwc = c.T @ (w[:, None] * c)
    inv = np.linalg.pinv(ctwc)
    residual = y - mu
    rows: list[pd.DataFrame] = []
    for start in range(0, x.shape[1], batch_size):
        end = min(start + batch_size, x.shape[1])
        xb = x[:, start:end].astype(float, copy=False)
        projection = inv @ (c.T @ (w[:, None] * xb))
        xres = xb - c @ projection
        u = xres.T @ residual
        v = np.sum(w[:, None] * xres * xres, axis=0)
        valid = v > 1e-12
        beta = np.full(end - start, np.nan)
        se = np.full(end - start, np.nan)
        z = np.full(end - start, np.nan)
        beta[valid] = u[valid] / v[valid]
        se[valid] = 1.0 / np.sqrt(v[valid])
        z[valid] = u[valid] / np.sqrt(v[valid])
        p = 2 * norm.sf(np.abs(z))
        rows.append(pd.DataFrame({"feature_index": np.arange(start, end), "beta_score": beta, "se_score": se, "z_score": z, "p_score": p}))
    return pd.concat(rows, ignore_index=True)


def contingency(x: np.ndarray, y: np.ndarray) -> tuple[int, int, int, int]:
    present = x.astype(bool)
    resistant = y.astype(bool)
    return (
        int(np.sum(present & resistant)),
        int(np.sum(present & ~resistant)),
        int(np.sum(~present & resistant)),
        int(np.sum(~present & ~resistant)),
    )


def random_effects(x: np.ndarray, y: np.ndarray, groups: pd.Series) -> dict:
    effects = []
    group_values = groups.fillna("UNKNOWN").astype(str).to_numpy()
    for group in sorted(set(group_values)):
        idx = group_values == group
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
        effects.append((group, log_or, var, a, b, c, d))
    if not effects:
        return {"n_groups": 0}
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
    z = pooled / se if se > 0 else 0.0
    return {
        "n_groups": len(effects),
        "odds_ratio": math.exp(pooled),
        "ci_low": math.exp(pooled - 1.96 * se),
        "ci_high": math.exp(pooled + 1.96 * se),
        "one_sided_p": float(norm.sf(z)),
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 and df > 0 else 0.0,
        "positive_groups": int(sum(e[1] > 0 for e in effects)),
    }


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    dims = [int(v) for v in a.pc_dims.split(",") if v.strip()]

    manifest = pd.read_csv(a.manifest, dtype={"assembly_ID": str})
    required = {"assembly_ID", "phenotype", "split"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    manifest = manifest.drop_duplicates("assembly_ID").copy()
    manifest["y"] = manifest["phenotype"].astype(str).eq("R").astype(int)

    rtab = pd.read_csv(a.rtab, sep="\t", index_col=0)
    rtab.index = rtab.index.astype(str)
    rtab.columns = rtab.columns.astype(str)
    rtab = rtab.apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.uint8)

    distance = pd.read_csv(a.distance, sep="\t", index_col=0)
    distance.index = distance.index.astype(str)
    distance.columns = distance.columns.astype(str)

    ids = manifest["assembly_ID"].tolist()
    if set(ids) != set(rtab.columns) or set(ids) != set(distance.index) or set(distance.index) != set(distance.columns):
        raise ValueError("Manifest, Rtab, and distance sample IDs are not identical")
    rtab = rtab.loc[:, ids]
    distance = distance.loc[ids, ids]

    feature_meta = pd.read_csv(a.feature_meta, dtype=str)
    feature_name_column = next((c for c in ["variant", "feature", "feature_id", "name"] if c in feature_meta.columns), None)
    if feature_name_column is not None:
        feature_meta = feature_meta.drop_duplicates(feature_name_column).set_index(feature_name_column)
    else:
        feature_meta = pd.DataFrame(index=rtab.index)

    discovery_ids = manifest.loc[manifest["split"].eq("discovery"), "assembly_ID"].tolist()
    validation_ids = manifest.loc[manifest["split"].eq("validation"), "assembly_ID"].tolist()
    if set(discovery_ids) & set(validation_ids):
        raise ValueError("Discovery and validation overlap")
    if not discovery_ids or not validation_ids:
        raise ValueError("Empty discovery or validation split")

    x_all = rtab.T.to_numpy(dtype=np.uint8)
    y_all = manifest["y"].to_numpy(dtype=np.uint8)
    id_to_pos = {sample: i for i, sample in enumerate(ids)}
    dpos = np.array([id_to_pos[s] for s in discovery_ids])
    vpos = np.array([id_to_pos[s] for s in validation_ids])
    x_disc = x_all[dpos]
    y_disc = y_all[dpos]
    af = x_disc.mean(axis=0)
    eligible = (af >= a.min_af) & (af <= a.max_af)
    eligible_idx = np.flatnonzero(eligible)
    if len(eligible_idx) == 0:
        raise ValueError("No eligible features")

    model_tables = []
    disc_distance = distance.loc[discovery_ids, discovery_ids]
    for dim in dims:
        pcs = pcoa(disc_distance, dim)
        cov = np.column_stack([np.ones(len(discovery_ids)), pcs])
        stats = score_test_matrix(x_disc[:, eligible_idx], y_disc, cov)
        stats["feature"] = rtab.index.to_numpy()[eligible_idx[stats["feature_index"].to_numpy()]]
        stats["q_score"] = bh(stats["p_score"])
        stats["pc_dim"] = dim
        stats.to_csv(out / f"DISCOVERY_SCORE_GWAS_PC{dim}.csv", index=False)
        model_tables.append(stats[["feature", "pc_dim", "beta_score", "se_score", "z_score", "p_score", "q_score"]])

    long = pd.concat(model_tables, ignore_index=True)
    wide = long.pivot(index="feature", columns="pc_dim")
    wide.columns = [f"{metric}_pc{dim}" for metric, dim in wide.columns]
    wide = wide.reset_index()
    beta_cols = [f"beta_score_pc{d}" for d in dims]
    p_cols = [f"p_score_pc{d}" for d in dims]
    q_cols = [f"q_score_pc{d}" for d in dims]
    wide["beta_min"] = wide[beta_cols].min(axis=1)
    wide["beta_max"] = wide[beta_cols].max(axis=1)
    wide["p_max"] = wide[p_cols].max(axis=1)
    wide["q_max"] = wide[q_cols].max(axis=1)
    wide["all_positive"] = wide["beta_min"] > 0
    bonf = a.alpha / len(eligible_idx)
    wide["discovery_bonferroni_stable"] = wide["all_positive"] & (wide["p_max"] <= bonf)
    wide["discovery_fdr_stable"] = wide["all_positive"] & (wide["q_max"] <= a.alpha)

    selected = wide.loc[wide["discovery_fdr_stable"]].copy()
    selected = selected.sort_values(["discovery_bonferroni_stable", "q_max", "p_max"], ascending=[False, True, True])
    selected.to_csv(out / "FROZEN_DISCOVERY_CANDIDATES.csv", index=False)

    validation_manifest = manifest.iloc[vpos].reset_index(drop=True)
    validation_y = y_all[vpos]
    evidence_rows = []
    meta_details = {}
    feature_to_index = {f: i for i, f in enumerate(rtab.index)}
    for _, row in selected.iterrows():
        feature = str(row["feature"])
        idx = feature_to_index[feature]
        xd = x_disc[:, idx]
        xv = x_all[vpos, idx]
        xa = x_all[:, idx]
        da, db, dc, dd = contingency(xd, y_disc)
        va, vb, vc, vd = contingency(xv, validation_y)
        aa, ab, ac, ad = contingency(xa, y_all)
        dor, dlo, dhi = odds_ratio_ci(da, db, dc, dd)
        vor, vlo, vhi = odds_ratio_ci(va, vb, vc, vd)
        aor, alo, ahi = odds_ratio_ci(aa, ab, ac, ad)
        validation_p = float(fisher_exact([[va, vb], [vc, vd]], alternative="greater").pvalue)
        source = random_effects(xv, validation_y, validation_manifest.get("source_group", pd.Series(["UNKNOWN"] * len(validation_manifest))))
        country = random_effects(xv, validation_y, validation_manifest.get("ISO_country_code", pd.Series(["UNKNOWN"] * len(validation_manifest))))
        meta_details[feature] = {"source": source, "country": country}
        item = row.to_dict()
        item.update({
            "discovery_R_present": da, "discovery_S_present": db, "discovery_R_absent": dc, "discovery_S_absent": dd,
            "discovery_or": dor, "discovery_ci_low": dlo, "discovery_ci_high": dhi,
            "validation_R_present": va, "validation_S_present": vb, "validation_R_absent": vc, "validation_S_absent": vd,
            "validation_or": vor, "validation_ci_low": vlo, "validation_ci_high": vhi, "validation_p": validation_p,
            "all_or": aor, "all_ci_low": alo, "all_ci_high": ahi,
            "source_n_groups": source.get("n_groups", 0), "source_or": source.get("odds_ratio"), "source_ci_low": source.get("ci_low"), "source_ci_high": source.get("ci_high"), "source_p": source.get("one_sided_p"), "source_I2": source.get("I2_percent"),
            "country_n_groups": country.get("n_groups", 0), "country_or": country.get("odds_ratio"), "country_ci_low": country.get("ci_low"), "country_ci_high": country.get("ci_high"), "country_p": country.get("one_sided_p"), "country_I2": country.get("I2_percent"),
        })
        if feature in feature_meta.index:
            for c, value in feature_meta.loc[feature].items():
                item[f"meta_{c}"] = value
        evidence_rows.append(item)

    evidence = pd.DataFrame(evidence_rows)
    if len(evidence):
        evidence["validation_q"] = bh(evidence["validation_p"])
        evidence["heldout_replication"] = (
            (evidence["validation_R_present"] + evidence["validation_S_present"] >= a.min_validation_present)
            & (evidence["validation_or"] >= a.min_validation_odds_ratio)
            & (evidence["validation_ci_low"] > 1)
            & (evidence["validation_q"] <= a.alpha)
        )
        evidence["source_replication"] = (
            (evidence["source_n_groups"] >= 3)
            & (pd.to_numeric(evidence["source_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["source_p"], errors="coerce") <= a.alpha)
        )
        evidence["country_replication"] = (
            (evidence["country_n_groups"] >= 3)
            & (pd.to_numeric(evidence["country_ci_low"], errors="coerce") > 1)
            & (pd.to_numeric(evidence["country_p"], errors="coerce") <= a.alpha)
        )
        evidence["strict_statistical_gate"] = evidence["heldout_replication"] & evidence["source_replication"] & evidence["country_replication"]
        evidence = evidence.sort_values(["strict_statistical_gate", "heldout_replication", "validation_q", "q_max"], ascending=[False, False, True, True])
    else:
        for column in ["validation_q", "heldout_replication", "source_replication", "country_replication", "strict_statistical_gate"]:
            evidence[column] = pd.Series(dtype=float if column == "validation_q" else bool)

    evidence.to_csv(out / "ALL_DISCOVERY_TO_VALIDATION_EVIDENCE.csv", index=False)
    strict = evidence.loc[evidence["strict_statistical_gate"].eq(True)].copy() if len(evidence) else evidence.copy()
    strict.to_csv(out / "STRICT_STATISTICALLY_REPLICATED_TARGETED_FEATURES.csv", index=False)
    (out / "META_ANALYSIS_DETAILS.json").write_text(json.dumps(meta_details, indent=2, ensure_ascii=False, default=str) + "\n")

    summary = {
        "n_all_samples": int(len(ids)),
        "n_discovery": int(len(discovery_ids)),
        "n_validation": int(len(validation_ids)),
        "n_total_features": int(len(rtab)),
        "n_eligible_features": int(len(eligible_idx)),
        "pc_dimensions": dims,
        "bonferroni_threshold": bonf,
        "n_discovery_fdr_stable": int(len(selected)),
        "n_discovery_bonferroni_stable": int(selected["discovery_bonferroni_stable"].sum()) if len(selected) else 0,
        "n_heldout_replicated": int(evidence["heldout_replication"].sum()) if len(evidence) else 0,
        "n_source_replicated": int(evidence["source_replication"].sum()) if len(evidence) else 0,
        "n_country_replicated": int(evidence["country_replication"].sum()) if len(evidence) else 0,
        "n_strict_statistically_replicated": int(len(strict)),
        "strict_features": strict["feature"].astype(str).tolist() if len(strict) else [],
        "status": "TARGETED_STATISTICAL_CANDIDATES_REQUIRE_PYSEER_CONCORDANCE_AND_KNOWN_MECHANISM_AUDIT" if len(strict) else "NO_TARGETED_FEATURE_SURVIVED_COMPLETE_STATISTICAL_GATE",
        "claim_boundary": "A statistically replicated feature is not a novel resistance determinant. It must agree with an independent GWAS implementation, survive known-mechanism, lineage, assembly-context, database and literature audits, and still requires biological validation for causality.",
    }
    (out / "STRUCTURE_ADJUSTED_SCORE_GWAS_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    report = [
        "# Structure-adjusted targeted colistin GWAS",
        "",
        f"- All samples: **{summary['n_all_samples']:,}**",
        f"- Discovery: **{summary['n_discovery']:,}**",
        f"- Untouched validation: **{summary['n_validation']:,}**",
        f"- Eligible targeted features: **{summary['n_eligible_features']:,}**",
        f"- Stable discovery FDR candidates: **{summary['n_discovery_fdr_stable']:,}**",
        f"- Held-out replicates: **{summary['n_heldout_replicated']:,}**",
        f"- Complete statistical gate: **{summary['n_strict_statistically_replicated']:,}**",
        "",
        "## Claim boundary",
        "",
        summary["claim_boundary"],
    ]
    if len(evidence):
        cols = [c for c in ["feature", "q_max", "validation_or", "validation_ci_low", "validation_ci_high", "validation_q", "source_n_groups", "source_or", "source_ci_low", "country_n_groups", "country_or", "country_ci_low", "strict_statistical_gate"] if c in evidence.columns]
        report.extend(["", "## Evidence-ranked candidates", "", evidence.head(50)[cols].to_markdown(index=False)])
    (out / "STRUCTURE_ADJUSTED_SCORE_GWAS_REPORT.md").write_text("\n".join(report) + "\n")

    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
