#!/usr/bin/env python3
"""Train a frozen precision filter from one completed Toy development pair.

The training labels are opened only after the development assembly was frozen. A
candidate is labelled safe when its appended contig has an accepted gold alignment and
is not flagged by the conservative cross-BINID chimera audit. The model is intended to
filter a *separately frozen* baseline-preserving assembly on a different participant
pair. It never receives sequence truth, taxonomic labels, or holdout performance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

FEATURES = [
    "length",
    "timepoints_supported",
    "assembler_count",
    "max_short_breadth",
    "max_long_breadth",
    "long_spanning_reads",
    "score",
    "best_backbone_identity",
    "best_backbone_shorter_coverage",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--decisions", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--chimera-audit", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=20260820)
    p.add_argument("--target-precision", type=float, default=0.80)
    p.add_argument("--minimum-threshold", type=float, default=0.50)
    return p.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    decisions_path = Path(a.decisions)
    metadata_path = Path(a.metadata)
    audit_path = Path(a.chimera_audit)

    decisions = pd.read_csv(decisions_path, sep="\t", dtype={"contig_id": str})
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"representative_id": str})
    audit = pd.read_csv(audit_path, sep="\t", dtype={"assembly_contig": str})

    selected = decisions[decisions["selected"].astype(str).str.lower().eq("true")].copy()
    metadata["output_id"] = "LANTERN_AUG_" + metadata["sequence_sha256"].str[:16].str.upper()
    selected = selected.merge(
        metadata[["representative_id", "sequence_sha256", "output_id"]],
        left_on="contig_id",
        right_on="representative_id",
        how="left",
        validate="one_to_one",
    )
    added_audit = audit[audit["assembly_contig"].str.startswith("LANTERN_AUG_")].copy()
    selected = selected.merge(
        added_audit,
        left_on="output_id",
        right_on="assembly_contig",
        how="left",
        validate="one_to_one",
        suffixes=("", "_gold"),
    )
    if selected["output_id"].isna().any() or selected["assembly_contig"].isna().any():
        raise SystemExit("selected candidate/output audit mapping is incomplete")

    selected["safe_label"] = (
        selected["has_accepted_alignment"].astype(str).str.lower().eq("true")
        & ~selected["cross_binid_chimera"].astype(str).str.lower().eq("true")
    )
    x = selected[FEATURES].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = selected["safe_label"].astype(int).to_numpy()
    if len(selected) < 100 or y.sum() < 25 or (len(y) - y.sum()) < 25:
        raise SystemExit("insufficient labelled candidates")

    primary = selected["primary_genome"].fillna("").astype(str)
    groups = np.array(
        [g if g else f"UNALIGNED_{oid}" for g, oid in zip(primary, selected["output_id"])],
        dtype=object,
    )
    model = HistGradientBoostingClassifier(
        max_iter=150,
        max_depth=3,
        learning_rate=0.04,
        l2_regularization=2.0,
        min_samples_leaf=20,
        random_state=a.seed,
    )
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=a.seed)
    oof = cross_val_predict(
        model,
        x,
        y,
        groups=groups,
        cv=splitter,
        method="predict_proba",
        n_jobs=-1,
    )[:, 1]
    precision, recall, thresholds = precision_recall_curve(y, oof)
    eligible = [
        i
        for i in range(len(thresholds))
        if precision[i] >= a.target_precision and thresholds[i] >= a.minimum_threshold
    ]
    if not eligible:
        raise SystemExit("no OOF threshold satisfies the frozen precision floor")
    index = max(eligible, key=lambda i: (recall[i], precision[i], thresholds[i]))
    threshold = float(thresholds[index])
    predicted = oof >= threshold
    observed_precision = float(y[predicted].mean()) if predicted.any() else 0.0
    observed_recall = float(y[predicted].sum() / y.sum())

    model.fit(x, y)
    model_path = out / "LANTERN_V5_PRECISION_FILTER.joblib"
    joblib.dump(model, model_path, compress=3)
    labelled = selected[
        [
            "contig_id",
            "output_id",
            "safe_label",
            "has_accepted_alignment",
            "cross_binid_chimera",
            "primary_genome",
            "primary_aligned_bp",
            *FEATURES,
        ]
    ].copy()
    labelled["oof_probability"] = oof
    labelled["oof_selected"] = predicted
    labelled.to_csv(out / "DEVELOPMENT_LABELS_AND_OOF.tsv", sep="\t", index=False)

    report = {
        "status": "MODEL_FROZEN_FOR_UNTOUCHED_APPLICATION",
        "version": "LANTERN-v5-precision-filter-20260820",
        "features": FEATURES,
        "model": {
            "class": "sklearn.ensemble.HistGradientBoostingClassifier",
            "max_iter": 150,
            "max_depth": 3,
            "learning_rate": 0.04,
            "l2_regularization": 2.0,
            "min_samples_leaf": 20,
            "random_state": a.seed,
        },
        "threshold": threshold,
        "threshold_selection": {
            "target_oof_precision": a.target_precision,
            "minimum_threshold": a.minimum_threshold,
            "observed_oof_precision": observed_precision,
            "observed_oof_recall": observed_recall,
        },
        "n_candidates": int(len(y)),
        "n_safe": int(y.sum()),
        "base_safe_fraction": float(y.mean()),
        "oof_roc_auc": float(roc_auc_score(y, oof)),
        "oof_average_precision": float(average_precision_score(y, oof)),
        "model_sha256": sha256(model_path),
        "input_sha256": {
            "decisions": sha256(decisions_path),
            "metadata": sha256(metadata_path),
            "chimera_audit": sha256(audit_path),
        },
        "claim_boundary": (
            "The model was trained on public Toy development gold after its assembly freeze. "
            "It is not evidence of holdout or CAMI challenge performance."
        ),
    }
    (out / "MODEL_FREEZE.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (out / "MODEL_FREEZE.md").write_text(
        "# LANTERN v5 precision-filter freeze\n\n"
        f"- Status: **{report['status']}**\n"
        f"- Labelled candidates: **{report['n_candidates']}**\n"
        f"- Safe candidates: **{report['n_safe']}**\n"
        f"- OOF ROC AUC: **{report['oof_roc_auc']:.4f}**\n"
        f"- OOF average precision: **{report['oof_average_precision']:.4f}**\n"
        f"- Frozen probability threshold: **{threshold:.8f}**\n"
        f"- OOF precision/recall at threshold: **{observed_precision:.4f}/{observed_recall:.4f}**\n\n"
        + report["claim_boundary"]
        + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
