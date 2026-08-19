#!/usr/bin/env python3
"""Evaluate discovery-selected unitigs in an untouched source-held-out cohort.

This is a statistical replication gate only. It does not establish causality, biological
novelty, clinical validity, or resistance mechanism status. Candidates must subsequently be
intersected with the independent sequence-level known-mechanism audit and context review.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--selection", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--validation-rtab", required=True)
    p.add_argument("--all-rtab", required=True)
    p.add_argument("--whole-pyseer", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min-validation-present", type=int, default=5)
    return p.parse_args()


def bh(s):
    p = pd.to_numeric(s, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    if len(v):
        order = np.argsort(v)
        ranked = v[order]
        q = np.minimum.accumulate(
            (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
        )[::-1]
        out[np.flatnonzero(ok)[order]] = np.clip(q, 0, 1)
    return pd.Series(out, index=s.index)


def rtab(path):
    x = pd.read_csv(path, sep="\t", index_col=0)
    x.index = x.index.astype(str).str.upper()
    x.columns = x.columns.astype(str)
    if x.index.duplicated().any():
        raise RuntimeError(f"Duplicate unitig sequences in {path}")
    if x.columns.duplicated().any():
        raise RuntimeError(f"Duplicate assembly IDs in {path}")
    return x.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def cont(row, meta):
    z = meta[["assembly_ID", "phenotype"]].copy()
    z["x"] = z.assembly_ID.map(row).fillna(0).astype(int)
    resistant = z.phenotype.astype(str).eq("R")
    return (
        int((resistant & (z.x == 1)).sum()),
        int((~resistant & (z.x == 1)).sum()),
        int((resistant & (z.x == 0)).sum()),
        int((~resistant & (z.x == 0)).sum()),
    )


def orci(a, b, c, d):
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    log_or = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return (
        math.exp(log_or),
        math.exp(log_or - 1.96 * se),
        math.exp(log_or + 1.96 * se),
    )


def random_meta(row, meta, col):
    effects = []
    for group, sub in meta.groupby(col, dropna=False):
        a, b, c, d = cont(row, sub)
        if (a + c) == 0 or (b + d) == 0 or (a + b) == 0 or (c + d) == 0:
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
        return {"n": 0}

    y = np.array([e[1] for e in effects])
    w = 1 / np.array([e[2] for e in effects])
    fixed = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fixed) ** 2))
    df = len(y) - 1
    cval = float(np.sum(w) - np.sum(w * w) / np.sum(w))
    tau2 = max(0, (q - df) / cval) if df > 0 and cval > 0 else 0
    wr = 1 / (np.array([e[2] for e in effects]) + tau2)
    random_log_or = float(np.sum(wr * y) / np.sum(wr))
    se = math.sqrt(1 / float(np.sum(wr)))
    return {
        "n": len(effects),
        "or": math.exp(random_log_or),
        "lo": math.exp(random_log_or - 1.96 * se),
        "hi": math.exp(random_log_or + 1.96 * se),
        "p": float(norm.sf(random_log_or / se)) if se else 1.0,
        "I2": max(0, (q - df) / q * 100) if q > 0 and df > 0 else 0,
        "positive": int(sum(e[1] > 0 for e in effects)),
        "details": [
            {
                "group": e[0],
                "log_or": e[1],
                "var": e[2],
                "a": e[3],
                "b": e[4],
                "c": e[5],
                "d": e[6],
            }
            for e in effects
        ],
    }


def pyseer(path, i):
    x = pd.read_csv(path, sep="\t", dtype={"variant": str})
    x.variant = x.variant.astype(str).str.upper()
    x["beta"] = pd.to_numeric(x["beta"], errors="coerce")
    x["lrt-pvalue"] = pd.to_numeric(x["lrt-pvalue"], errors="coerce")
    x["q"] = bh(x["lrt-pvalue"])
    return (
        x[["variant", "beta", "lrt-pvalue", "q"]]
        .drop_duplicates("variant")
        .rename(
            columns={
                "beta": f"whole_beta_{i}",
                "lrt-pvalue": f"whole_p_{i}",
                "q": f"whole_q_{i}",
            }
        )
    )


def main():
    a = args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    sel = pd.read_csv(
        a.selection, dtype={"canonical_sequence": str, "candidate_id": str}
    )
    sel.canonical_sequence = sel.canonical_sequence.str.upper()
    if sel.canonical_sequence.duplicated().any():
        raise RuntimeError("Duplicate canonical unitigs in discovery selection")

    meta = pd.read_csv(a.manifest, dtype={"assembly_ID": str})
    if meta.assembly_ID.duplicated().any():
        raise RuntimeError("Duplicate assembly IDs in sample manifest")
    val = meta[meta.split.eq("validation")].copy()
    vr = rtab(a.validation_rtab)
    ar = rtab(a.all_rtab)

    missing_validation = sorted(set(val.assembly_ID) - set(vr.columns))
    missing_all = sorted(set(meta.assembly_ID) - set(ar.columns))
    if missing_validation or missing_all:
        raise RuntimeError(
            f"Rtab/manifest mismatch: validation_missing={missing_validation[:10]} "
            f"all_missing={missing_all[:10]}"
        )

    whole = None
    for i, path in enumerate(a.whole_pyseer):
        current = pyseer(path, i)
        whole = current if whole is None else whole.merge(current, on="variant", how="outer")
    whole_beta_cols = [c for c in whole if c.startswith("whole_beta_")]
    whole_q_cols = [c for c in whole if c.startswith("whole_q_")]
    whole["whole_models"] = whole[whole_beta_cols].notna().sum(axis=1)
    whole["whole_beta_min"] = whole[whole_beta_cols].min(axis=1)
    whole["whole_q_max"] = whole[whole_q_cols].max(axis=1)

    rows = []
    detail = {}
    for _, selected in sel.iterrows():
        seq = selected.canonical_sequence
        if seq not in vr.index or seq not in ar.index:
            continue
        validation_table = cont(vr.loc[seq], val)
        all_table = cont(ar.loc[seq], meta)
        validation_or, validation_lo, validation_hi = orci(*validation_table)
        all_or, all_lo, all_hi = orci(*all_table)
        validation_p = float(
            fisher_exact(
                [
                    [validation_table[0], validation_table[1]],
                    [validation_table[2], validation_table[3]],
                ],
                alternative="greater",
            ).pvalue
        )
        source = random_meta(ar.loc[seq], meta, "source_group")
        country = random_meta(ar.loc[seq], meta, "ISO_country_code")
        detail[str(selected.candidate_id)] = {"source": source, "country": country}
        row = selected.to_dict()
        row.update(
            {
                "variant": seq,
                "validation_R_present": validation_table[0],
                "validation_S_present": validation_table[1],
                "validation_R_absent": validation_table[2],
                "validation_S_absent": validation_table[3],
                "validation_or": validation_or,
                "validation_ci_low": validation_lo,
                "validation_ci_high": validation_hi,
                "validation_p": validation_p,
                "all_or": all_or,
                "all_ci_low": all_lo,
                "all_ci_high": all_hi,
                "source_n": source.get("n", 0),
                "source_or": source.get("or"),
                "source_ci_low": source.get("lo"),
                "source_ci_high": source.get("hi"),
                "source_p": source.get("p"),
                "source_I2": source.get("I2"),
                "country_n": country.get("n", 0),
                "country_or": country.get("or"),
                "country_ci_low": country.get("lo"),
                "country_ci_high": country.get("hi"),
                "country_p": country.get("p"),
            }
        )
        rows.append(row)

    res = pd.DataFrame(rows)
    if res.empty:
        raise RuntimeError("No selected unitig overlapped validation and all matrices")

    res["validation_q"] = bh(res.validation_p)
    res = res.merge(whole, on="variant", how="left", validate="one_to_one")
    nmodels = len(a.whole_pyseer)
    res["validation_replication"] = (
        (res.validation_R_present + res.validation_S_present >= a.min_validation_present)
        & (res.validation_or > 1)
        & (res.validation_ci_low > 1)
        & (res.validation_q <= a.alpha)
    )
    res["whole_adjusted_stable"] = (
        (res.whole_models == nmodels)
        & (res.whole_beta_min > 0)
        & (res.whole_q_max <= a.alpha)
    )
    res["source_replication"] = (
        (res.source_n >= 3)
        & (pd.to_numeric(res.source_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(res.source_p, errors="coerce") <= a.alpha)
    )
    res["country_replication"] = (
        (res.country_n >= 3)
        & (pd.to_numeric(res.country_ci_low, errors="coerce") > 1)
        & (pd.to_numeric(res.country_p, errors="coerce") <= a.alpha)
    )
    res["strict_statistical_replication"] = (
        res.validation_replication
        & res.whole_adjusted_stable
        & res.source_replication
        & res.country_replication
    )
    res["score"] = (
        4 * res.validation_replication.astype(int)
        + 3 * res.whole_adjusted_stable.astype(int)
        + 4 * res.source_replication.astype(int)
        + 3 * res.country_replication.astype(int)
    )
    res = res.sort_values(
        ["strict_statistical_replication", "score", "validation_q", "whole_q_max"],
        ascending=[False, False, True, True],
    )

    res.to_csv(out / "ALL_UNITIG_REPLICATION_EVIDENCE.csv", index=False)
    strict = res[res.strict_statistical_replication].copy()
    strict.to_csv(out / "STRICT_STATISTICALLY_REPLICATED_UNITIGS.csv", index=False)
    (out / "UNITIG_META_ANALYSIS_DETAILS.json").write_text(
        json.dumps(detail, indent=2, ensure_ascii=False, default=str) + "\n"
    )
    summary = {
        "n_selected": int(len(sel)),
        "n_evaluated": int(len(res)),
        "n_validation_replicated": int(res.validation_replication.sum()),
        "n_whole_adjusted_stable": int(res.whole_adjusted_stable.sum()),
        "n_source_replicated": int(res.source_replication.sum()),
        "n_country_replicated": int(res.country_replication.sum()),
        "n_strict_statistically_replicated": int(
            res.strict_statistical_replication.sum()
        ),
        "strict_candidate_ids": strict.candidate_id.astype(str).tolist(),
        "status": (
            "PORTAL_COHORT_UNITIGS_REQUIRE_KNOWN_MECHANISM_INTERSECTION"
            if len(strict)
            else "NO_UNITIG_SURVIVED_COMPLETE_PORTAL_COHORT_GATE"
        ),
        "boundary": (
            "Statistical replication in the portal-residual cohort is not a novel resistance "
            "determinant. Candidates must survive the independent sequence-level known-mechanism "
            "filter, context/database/literature review, and biological validation."
        ),
    }
    (out / "UNITIG_REPLICATION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    report_columns = [
        c
        for c in [
            "candidate_id",
            "sequence_length",
            "p_max",
            "q_filtered_max",
            "validation_or",
            "validation_ci_low",
            "validation_ci_high",
            "validation_q",
            "beta_min",
            "whole_beta_min",
            "whole_q_max",
            "source_n",
            "source_or",
            "source_ci_low",
            "source_ci_high",
            "country_n",
            "country_or",
            "country_ci_low",
            "country_ci_high",
            "strict_statistical_replication",
        ]
        if c in res
    ]
    report = [
        "# Portal-residual whole-genome unitig replication audit",
        "",
        f"- Discovery-selected unitigs: **{len(sel):,}**",
        f"- Evaluated in untouched validation: **{len(res):,}**",
        f"- Held-out replicates: **{int(res.validation_replication.sum()):,}**",
        f"- Stable across adjusted models: **{int(res.whole_adjusted_stable.sum()):,}**",
        f"- Cross-source replicates: **{int(res.source_replication.sum()):,}**",
        f"- Cross-country replicates: **{int(res.country_replication.sum()):,}**",
        f"- Complete statistical gate: **{int(res.strict_statistical_replication.sum()):,}**",
        "",
        "## Claim boundary",
        "",
        summary["boundary"],
        "",
        "## Top evidence-ranked unitigs",
        "",
        res.head(30)[report_columns].to_markdown(index=False),
    ]
    (out / "UNITIG_REPLICATION_REPORT.md").write_text("\n".join(report) + "\n")
    hashes = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}"
        for path in sorted(out.rglob("*"))
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
