#!/usr/bin/env python3
"""Train/apply a development-only precision filter for LANTERN whole extensions.

Gold labels are allowed only in ``train`` mode on the declared development data.  The
serialized model and numerical thresholds are subsequently frozen and may be applied to
an untouched holdout before its gold is opened.  Passing this script is not a CAMI rank,
a hidden-challenge result, or evidence of a novel organism.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterator

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

DNA = "ACGT"
K3 = [a + b + c for a in DNA for b in DNA for c in DNA]
FEATURE_VERSION = "LANTERN-v8-sequence-filter-1"


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    name = None
    seq: list[str] = []
    with open(path, "rt", encoding="ascii") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ValueError(f"empty FASTA id at line {line_no}")
                seq = []
            else:
                if name is None:
                    raise ValueError(f"sequence before header at line {line_no}")
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def write_record(fh, name: str, seq: str, width: int = 80) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def entropy(counter: Counter[str], total: int) -> float:
    if total <= 0:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counter.values() if n)


def longest_run(seq: str) -> int:
    if not seq:
        return 0
    best = current = 1
    previous = seq[0]
    for char in seq[1:]:
        if char == previous:
            current += 1
            best = max(best, current)
        else:
            previous = char
            current = 1
    return best


def seq_features(seq: str) -> dict[str, float]:
    n = len(seq)
    c1 = Counter(seq)
    c2 = Counter(seq[i:i + 2] for i in range(max(0, n - 1)))
    c3 = Counter(seq[i:i + 3] for i in range(max(0, n - 2)))
    denom3 = max(1, n - 2)
    run = longest_run(seq)
    row: dict[str, float] = {
        "length": float(n),
        "log_length": math.log1p(n),
        "gc_fraction": (c1.get("G", 0) + c1.get("C", 0)) / max(1, n),
        "n_fraction": c1.get("N", 0) / max(1, n),
        "mono_entropy": entropy(c1, n),
        "dinucleotide_entropy": entropy(c2, max(0, n - 1)),
        "max_homopolymer": float(run),
        "max_homopolymer_fraction": run / max(1, n),
        "distinct_5mer_fraction": len({seq[i:i + 5] for i in range(max(0, n - 4))}) / max(1, min(max(0, n - 4), 1024)),
        "top_3mer_fraction": max(c3.values()) / denom3 if c3 else 0.0,
    }
    for kmer in K3:
        row[f"k3_{kmer}"] = c3.get(kmer, 0) / denom3
    return row


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def feature_frame(assembly_fasta: str | Path, segments_tsv: str | Path) -> tuple[pd.DataFrame, dict[str, str]]:
    seqs = dict(read_fasta(assembly_fasta))
    segments = pd.read_csv(segments_tsv, sep="\t", dtype=str)
    if "output_id" not in segments:
        raise ValueError("extension segments lack output_id")
    rows: list[dict[str, object]] = []
    extension_sequences: dict[str, str] = {}
    for segment in segments.to_dict("records"):
        output_id = segment["output_id"]
        sequence = seqs.get(output_id)
        if sequence is None:
            raise ValueError(f"extension output absent from FASTA: {output_id}")
        extension_sequences[output_id] = sequence
        row: dict[str, object] = seq_features(sequence)
        row["output_id"] = output_id
        row["segment_is_whole"] = float(str(segment.get("segment_type", "")) == "whole")
        rows.append(row)
    frame = pd.DataFrame(rows).set_index("output_id", drop=False)
    return frame, extension_sequences


def approximate_chimera_change(selected: pd.DataFrame, baseline: dict[str, object]) -> float | None:
    base_bp = float(baseline["assembly_total_bp"])
    base_chimera_bp = float(baseline["cross_binid_chimeric_bp"])
    base_fraction = base_chimera_bp / base_bp if base_bp else 0.0
    selected_bp = float(selected["audit_length"].sum())
    selected_chimera_bp = float(selected.loc[selected["label_chimera"], "audit_length"].sum())
    new_fraction = (base_chimera_bp + selected_chimera_bp) / (base_bp + selected_bp) if base_bp + selected_bp else 0.0
    if base_fraction == 0:
        return 0.0 if new_fraction == 0 else None
    return (new_fraction - base_fraction) / base_fraction


def make_model(seed: int) -> ExtraTreesClassifier:
    return ExtraTreesClassifier(
        n_estimators=350,
        min_samples_leaf=8,
        max_features=0.70,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )


def train(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    features, _ = feature_frame(args.assembly_fasta, args.segments)
    audit = pd.read_csv(args.chimera_audit, sep="\t", dtype=str).rename(columns={"assembly_contig": "output_id"})
    required = {"output_id", "length", "cross_binid_chimera", "has_accepted_alignment", "primary_aligned_bp"}
    if not required.issubset(audit.columns):
        raise ValueError(f"chimera audit missing: {sorted(required - set(audit.columns))}")
    data = features.merge(audit, on="output_id", how="inner", validate="one_to_one")
    if len(data) != len(features):
        raise ValueError(f"audit/feature mismatch: {len(data)} vs {len(features)}")
    data["label_chimera"] = data["cross_binid_chimera"].map(bool_value)
    data["label_aligned"] = data["has_accepted_alignment"].map(bool_value)
    data["label_safe_aligned"] = data["label_aligned"] & ~data["label_chimera"]
    data["audit_length"] = pd.to_numeric(data["length_y"] if "length_y" in data else data["length"], errors="raise")
    data["primary_aligned_bp"] = pd.to_numeric(data["primary_aligned_bp"], errors="coerce").fillna(0.0)
    feature_columns = [c for c in features.columns if c != "output_id"]
    X = data[feature_columns].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    y_chimera = data["label_chimera"].astype(int)
    y_aligned = data["label_aligned"].astype(int)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.seed)
    chimera_model = make_model(args.seed)
    aligned_model = make_model(args.seed + 1)
    p_chimera = cross_val_predict(chimera_model, X, y_chimera, cv=cv, method="predict_proba", n_jobs=1)[:, 1]
    p_aligned = cross_val_predict(aligned_model, X, y_aligned, cv=cv, method="predict_proba", n_jobs=1)[:, 1]
    data["oof_p_chimera"] = p_chimera
    data["oof_p_aligned"] = p_aligned
    baseline = json.loads(Path(args.baseline_summary).read_text())
    grid: list[dict[str, object]] = []
    winner: dict[str, object] | None = None
    for chimera_max in np.arange(0.15, 0.91, 0.025):
        for aligned_min in np.arange(0.10, 0.91, 0.05):
            selected = data[(data["oof_p_chimera"] <= chimera_max) & (data["oof_p_aligned"] >= aligned_min)]
            if len(selected) < args.minimum_selected:
                continue
            relative_chimera = approximate_chimera_change(selected, baseline)
            useful_bp = float(selected.loc[selected["label_safe_aligned"], "primary_aligned_bp"].sum())
            aligned_bp = float(selected.loc[selected["label_aligned"], "primary_aligned_bp"].sum())
            row = {
                "chimera_probability_maximum": float(chimera_max),
                "aligned_probability_minimum": float(aligned_min),
                "selected_extensions": int(len(selected)),
                "selected_bp": int(selected["audit_length"].sum()),
                "selected_chimeric_extensions": int(selected["label_chimera"].sum()),
                "selected_chimeric_bp": int(selected.loc[selected["label_chimera"], "audit_length"].sum()),
                "selected_safe_aligned_extensions": int(selected["label_safe_aligned"].sum()),
                "selected_safe_primary_aligned_bp": useful_bp,
                "selected_primary_aligned_bp": aligned_bp,
                "approximate_relative_chimera_change": relative_chimera,
                "oof_safe_precision": float(selected["label_safe_aligned"].mean()),
            }
            grid.append(row)
            safe = relative_chimera is not None and relative_chimera <= args.development_chimera_margin
            if safe:
                key = (useful_bp, aligned_bp, -float(relative_chimera), -len(selected))
                if winner is None or key > winner["_key"]:
                    winner = {**row, "_key": key}
    if winner is None:
        raise RuntimeError("no OOF threshold satisfied the development chimera margin")
    winner.pop("_key", None)
    chimera_model.fit(X, y_chimera)
    aligned_model.fit(X, y_aligned)
    bundle = {
        "feature_version": FEATURE_VERSION,
        "feature_columns": feature_columns,
        "chimera_model": chimera_model,
        "aligned_model": aligned_model,
        "chimera_probability_maximum": winner["chimera_probability_maximum"],
        "aligned_probability_minimum": winner["aligned_probability_minimum"],
        "training_rows": len(data),
        "training_seed": args.seed,
    }
    joblib.dump(bundle, out / "LANTERN_V8_WHOLE_FILTER.joblib", compress=3)
    data.to_csv(out / "TRAINING_OOF_PREDICTIONS.tsv", sep="\t", index=False)
    pd.DataFrame(grid).to_csv(out / "THRESHOLD_GRID.tsv", sep="\t", index=False)
    metrics = {
        "feature_version": FEATURE_VERSION,
        "n_training_extensions": int(len(data)),
        "n_chimeric": int(y_chimera.sum()),
        "n_aligned": int(y_aligned.sum()),
        "n_safe_aligned": int(data["label_safe_aligned"].sum()),
        "chimera_roc_auc": float(roc_auc_score(y_chimera, p_chimera)),
        "chimera_average_precision": float(average_precision_score(y_chimera, p_chimera)),
        "aligned_roc_auc": float(roc_auc_score(y_aligned, p_aligned)),
        "aligned_average_precision": float(average_precision_score(y_aligned, p_aligned)),
        "selected_threshold": winner,
        "development_chimera_margin": args.development_chimera_margin,
        "claim_boundary": "Candidate labels and threshold selection use development gold only. Holdout accuracy is unknown until the serialized model is applied unchanged before gold access.",
    }
    (out / "MODEL_FREEZE.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    (out / "MODEL_SHA256.txt").write_text(f"{sha256(out / 'LANTERN_V8_WHOLE_FILTER.joblib')}  LANTERN_V8_WHOLE_FILTER.joblib\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))


def apply(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    bundle = joblib.load(args.model)
    if bundle.get("feature_version") != FEATURE_VERSION:
        raise ValueError("feature-version mismatch")
    features, extension_sequences = feature_frame(args.assembly_fasta, args.segments)
    columns = bundle["feature_columns"]
    missing = sorted(set(columns) - set(features.columns))
    if missing:
        raise ValueError(f"missing model features: {missing[:20]}")
    X = features[columns].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    features["predicted_chimera_probability"] = bundle["chimera_model"].predict_proba(X)[:, 1]
    features["predicted_aligned_probability"] = bundle["aligned_model"].predict_proba(X)[:, 1]
    features["selected"] = (
        (features["predicted_chimera_probability"] <= float(bundle["chimera_probability_maximum"]))
        & (features["predicted_aligned_probability"] >= float(bundle["aligned_probability_minimum"]))
    )
    selected_ids = list(features.index[features["selected"]])
    backbone = list(read_fasta(args.backbone_fasta))
    if not backbone:
        raise ValueError("empty backbone")
    seen_hashes = {hashlib.sha256(seq.encode()).hexdigest() for _, seq in backbone}
    written: list[str] = []
    with open(out / "LANTERN_V8_FILTERED_ASSEMBLY.fasta", "wt", encoding="ascii") as fh:
        for name, seq in backbone:
            write_record(fh, name, seq)
        for output_id in selected_ids:
            seq = extension_sequences[output_id]
            digest = hashlib.sha256(seq.encode()).hexdigest()
            if digest in seen_hashes:
                continue
            seen_hashes.add(digest)
            write_record(fh, output_id, seq)
            written.append(output_id)
    features["written"] = features.index.isin(written)
    features.to_csv(out / "V8_FILTER_DECISIONS.tsv", sep="\t", index=False)
    summary = {
        "feature_version": FEATURE_VERSION,
        "model_sha256": sha256(args.model),
        "chimera_probability_maximum": float(bundle["chimera_probability_maximum"]),
        "aligned_probability_minimum": float(bundle["aligned_probability_minimum"]),
        "candidate_extensions": int(len(features)),
        "selected_extensions_before_exact_dedup": int(len(selected_ids)),
        "written_extensions": int(len(written)),
        "written_extension_bp": int(sum(len(extension_sequences[x]) for x in written)),
        "backbone_records": int(len(backbone)),
        "backbone_bp": int(sum(len(seq) for _, seq in backbone)),
        "truth_accessed": False,
        "claim_boundary": "Selection uses the development-frozen model and pre-gold sequence features only.",
    }
    (out / "V8_FILTER_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    t = sub.add_parser("train")
    t.add_argument("--assembly-fasta", required=True)
    t.add_argument("--segments", required=True)
    t.add_argument("--chimera-audit", required=True)
    t.add_argument("--baseline-summary", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--seed", type=int, default=20260820)
    t.add_argument("--development-chimera-margin", type=float, default=0.08)
    t.add_argument("--minimum-selected", type=int, default=100)
    t.set_defaults(func=train)
    a = sub.add_parser("apply")
    a.add_argument("--model", required=True)
    a.add_argument("--assembly-fasta", required=True)
    a.add_argument("--segments", required=True)
    a.add_argument("--backbone-fasta", required=True)
    a.add_argument("--out", required=True)
    a.set_defaults(func=apply)
    return p


if __name__ == "__main__":
    ns = parser().parse_args()
    ns.func(ns)
