#!/usr/bin/env python3
"""Apply a frozen LANTERN-v5 candidate precision model without gold access.

All backbone records are preserved. Only augmentation records selected by the specified
LANTERN ablation and passing the frozen model threshold remain. Candidate output IDs are
reconstructed from the sequence hashes already recorded before gold evaluation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterator

import joblib
import pandas as pd

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
DNA_RE = re.compile(r"^[ACGTN]+$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-assembly", required=True)
    p.add_argument("--decisions", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--model-freeze", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name = None
    parts: list[str] = []
    with path.open("rt", encoding="ascii") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split(None, 1)[0]
                parts = []
            else:
                if name is None:
                    raise ValueError("sequence before first FASTA header")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def write_record(fh, name: str, seq: str) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i : i + 80] + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    assembly_path = Path(a.input_assembly)
    decisions_path = Path(a.decisions)
    metadata_path = Path(a.metadata)
    model_path = Path(a.model)
    freeze_path = Path(a.model_freeze)

    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("features") != FEATURES:
        raise SystemExit("frozen model feature list mismatch")
    if sha256(model_path) != freeze.get("model_sha256"):
        raise SystemExit("frozen model SHA-256 mismatch")
    threshold = float(freeze["threshold"])
    model = joblib.load(model_path)

    decisions = pd.read_csv(decisions_path, sep="\t", dtype={"contig_id": str})
    metadata = pd.read_csv(metadata_path, sep="\t", dtype={"representative_id": str})
    metadata["output_id"] = "LANTERN_AUG_" + metadata["sequence_sha256"].str[:16].str.upper()
    decisions = decisions.merge(
        metadata[["representative_id", "sequence_sha256", "output_id"]],
        left_on="contig_id",
        right_on="representative_id",
        how="left",
        validate="one_to_one",
    )
    if decisions["output_id"].isna().any():
        raise SystemExit("candidate metadata mapping incomplete")
    x = decisions[FEATURES].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    probabilities = model.predict_proba(x)[:, 1]
    decisions["v5_probability"] = probabilities
    decisions["v4_selected"] = decisions["selected"].astype(str).str.lower().eq("true")
    decisions["v5_selected"] = decisions["v4_selected"] & (decisions["v5_probability"] >= threshold)
    keep_outputs = set(decisions.loc[decisions["v5_selected"], "output_id"])

    records = list(read_fasta(assembly_path))
    if not records:
        raise SystemExit("empty input assembly")
    backbone = [(name, seq) for name, seq in records if not name.startswith("LANTERN_AUG_")]
    augmented = [(name, seq) for name, seq in records if name.startswith("LANTERN_AUG_")]
    present_aug = {name for name, _ in augmented}
    missing = sorted(keep_outputs - present_aug)
    if missing:
        raise SystemExit(f"selected candidate records absent from input full assembly: {missing[:10]}")

    output_path = out / "LANTERN_V5_ASSEMBLY.fasta"
    seen: set[str] = set()
    with output_path.open("wt", encoding="ascii") as fh:
        for name, seq in backbone + [(name, seq) for name, seq in augmented if name in keep_outputs]:
            if name in seen:
                raise SystemExit(f"duplicate output FASTA ID: {name}")
            if not seq or not DNA_RE.fullmatch(seq):
                raise SystemExit(f"invalid DNA sequence: {name}")
            seen.add(name)
            write_record(fh, name, seq)

    decisions.to_csv(out / "V5_FILTER_DECISIONS.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS",
        "method": "LANTERN-v5 frozen precision filter",
        "input_records": len(records),
        "backbone_records_preserved": len(backbone),
        "input_augmentation_records": len(augmented),
        "v4_selected_candidates": int(decisions["v4_selected"].sum()),
        "v5_selected_candidates": int(decisions["v5_selected"].sum()),
        "output_records": len(seen),
        "output_bp": sum(len(seq) for name, seq in backbone + [(n, s) for n, s in augmented if n in keep_outputs]),
        "threshold": threshold,
        "model_sha256": sha256(model_path),
        "input_sha256": {
            "assembly": sha256(assembly_path),
            "decisions": sha256(decisions_path),
            "metadata": sha256(metadata_path),
            "model_freeze": sha256(freeze_path),
        },
        "truth_input_accepted": False,
        "claim_boundary": "Truth-blind filtering only; performance is evaluated only after output freeze.",
    }
    (out / "V5_FILTER_SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out / "ASSEMBLY_SHA256.txt").write_text(f"{sha256(output_path)}  {output_path.name}\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
