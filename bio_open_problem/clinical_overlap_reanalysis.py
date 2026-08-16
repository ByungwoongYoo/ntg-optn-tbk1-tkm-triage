#!/usr/bin/env python3
"""Correct the DMS/clinical protein-overlap classification without reopening labels.

This script operates on the already frozen/evaluated clinical artifact. Clinical
proteins are marked as previously seen if any of the following hold:
1. exact target-sequence match to a DMS training target;
2. one target sequence is a >=50-aa exact substring of the other;
3. a conservative entry-name/gene-symbol heuristic matches.

The union is the primary conservative exclusion. Sequence-only exclusion is a
prespecified sensitivity analysis.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260817


def clean_sequence(value) -> str:
    return re.sub(r"[^A-Z]", "", str(value).upper())


def dms_entry_symbol(value) -> str:
    text = str(value).upper().strip()
    return re.sub(r"_HUMAN$", "", text)


def clinical_symbol(row: pd.Series) -> str:
    for column in ["EVE_model_path", "MSA_filename"]:
        text = str(row.get(column, "")).strip()
        match = re.search(r"refseq-([A-Za-z0-9]+)-", text, re.I)
        if match:
            return match.group(1).upper()
        match = re.search(r"([A-Za-z0-9]+)_HUMAN", text, re.I)
        if match:
            return match.group(1).upper()
    return ""


def build_overlap_map(dms_reference: pd.DataFrame, clinical_reference: pd.DataFrame):
    dms_sequences = [
        clean_sequence(value) for value in dms_reference["target_seq"]
        if len(clean_sequence(value)) >= 50
    ]
    exact_sequences = set(dms_sequences)
    dms_symbols = set(dms_reference["UniProt_ID"].map(dms_entry_symbol))
    rows = []
    for _, row in clinical_reference.iterrows():
        sequence = clean_sequence(row.get("target_seq", ""))
        symbol = clinical_symbol(row)
        exact = sequence in exact_sequences if sequence else False
        containment = False
        if sequence and not exact:
            for dms_sequence in dms_sequences:
                shorter, longer = (
                    (sequence, dms_sequence)
                    if len(sequence) <= len(dms_sequence)
                    else (dms_sequence, sequence)
                )
                if len(shorter) >= 50 and shorter in longer:
                    containment = True
                    break
        symbol_match = bool(symbol and symbol in dms_symbols)
        rows.append({
            "DMS_id": str(row["DMS_id"]),
            "clinical_symbol": symbol,
            "exact_sequence_match": exact,
            "sequence_containment_match": containment,
            "entry_symbol_match": symbol_match,
            "seen_sequence_only": exact or containment,
            "seen_conservative_union": exact or containment or symbol_match,
        })
    return pd.DataFrame(rows)


def bootstrap(values: np.ndarray, n_boot: int = 50000):
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(SEED)
    draws = np.empty(n_boot)
    for index in range(n_boot):
        draws[index] = np.mean(rng.choice(values, size=len(values), replace=True))
    return [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))]


def sign_flip(values: np.ndarray, n_perm: int = 200000):
    values = values[np.isfinite(values)]
    observed = float(np.mean(values))
    rng = np.random.default_rng(SEED)
    count = 0
    for _ in range(n_perm):
        permuted = float(np.mean(values * rng.choice([-1.0, 1.0], size=len(values))))
        count += int(permuted >= observed)
    return float((count + 1) / (n_perm + 1))


def summarize(frame: pd.DataFrame):
    difference = (
        frame["prediction_ridge100"]
        - frame["prediction_dms_selected_base"]
    ).to_numpy(float)
    return {
        "n_genes": int(len(frame)),
        "ridge_mean_auc": float(frame["prediction_ridge100"].mean()),
        "dms_selected_base_mean_auc": float(frame["prediction_dms_selected_base"].mean()),
        "mean_difference": float(np.mean(difference)),
        "median_difference": float(np.median(difference)),
        "paired_bootstrap_95ci": bootstrap(difference),
        "one_sided_sign_flip_p": sign_flip(difference),
        "fraction_genes_improved": float(np.mean(difference > 0)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-aucs", required=True)
    parser.add_argument("--dms-reference", required=True)
    parser.add_argument("--clinical-reference", required=True)
    parser.add_argument("--prediction-freeze", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    aucs = pd.read_csv(args.gene_aucs)
    dms_reference = pd.read_csv(args.dms_reference)
    clinical_reference = pd.read_csv(args.clinical_reference)
    freeze = json.loads(Path(args.prediction_freeze).read_text())
    overlap = build_overlap_map(dms_reference, clinical_reference)
    merged = aucs.drop(columns=["UniProt_ID"], errors="ignore").merge(
        overlap, on="DMS_id", how="left", validate="one_to_one"
    )
    if merged["seen_conservative_union"].isna().any():
        missing = merged.loc[
            merged["seen_conservative_union"].isna(), "DMS_id"
        ].tolist()[:20]
        raise RuntimeError(f"Clinical reference mapping failed for {missing}")

    primary = merged[~merged["seen_conservative_union"]].copy()
    sequence_sensitivity = merged[~merged["seen_sequence_only"]].copy()
    all_genes = merged.copy()
    result = {
        "correction_reason": (
            "The first transfer run treated a missing clinical UniProt_ID column as an empty "
            "identifier, so every clinical protein was incorrectly marked unseen. This "
            "reanalysis replaces that field with sequence and entry-symbol overlap checks."
        ),
        "prediction_freeze": freeze,
        "n_all_evaluable_genes": int(len(merged)),
        "n_seen_exact_sequence": int(merged["exact_sequence_match"].sum()),
        "n_seen_sequence_containment": int(merged["sequence_containment_match"].sum()),
        "n_seen_entry_symbol": int(merged["entry_symbol_match"].sum()),
        "n_seen_conservative_union": int(merged["seen_conservative_union"].sum()),
        "primary_unseen_conservative_union": summarize(primary),
        "sensitivity_unseen_sequence_only": summarize(sequence_sensitivity),
        "secondary_all_genes": summarize(all_genes),
    }
    primary_stats = result["primary_unseen_conservative_union"]
    result["decision"] = (
        "ROBUST_CLINICAL_TRANSFER_ADVANCE_AFTER_OVERLAP_CORRECTION"
        if primary_stats["paired_bootstrap_95ci"][0] > 0
        and primary_stats["one_sided_sign_flip_p"] < 0.05
        else "NO_ROBUST_CLINICAL_TRANSFER_ADVANCE_AFTER_OVERLAP_CORRECTION"
    )
    merged.to_csv(out / "clinical_gene_auc_with_overlap.csv", index=False)
    overlap.to_csv(out / "clinical_protein_overlap_map.csv", index=False)
    (out / "result.json").write_text(json.dumps(result, indent=2))
    report = f"""# Clinical-transfer overlap correction

The initial clinical-transfer analysis incorrectly treated all 717 evaluable genes as
unseen because the clinical reference file has no `UniProt_ID` column. This correction
uses exact sequence, >=50-aa sequence containment, and a conservative entry-symbol
match.

## Primary conservative unseen-protein result

- All evaluable genes: **{len(merged)}**
- Conservatively excluded as DMS-seen: **{result['n_seen_conservative_union']}**
- Remaining unseen genes: **{primary_stats['n_genes']}**
- DMS-trained ridge mean AUC: **{primary_stats['ridge_mean_auc']:.4f}**
- DMS-selected individual baseline mean AUC: **{primary_stats['dms_selected_base_mean_auc']:.4f}**
- Difference: **{primary_stats['mean_difference']:+.4f}**
- Paired bootstrap 95% CI: **[{primary_stats['paired_bootstrap_95ci'][0]:+.4f}, {primary_stats['paired_bootstrap_95ci'][1]:+.4f}]**
- One-sided sign-flip p: **{primary_stats['one_sided_sign_flip_p']:.6f}**
- Genes improved: **{100*primary_stats['fraction_genes_improved']:.1f}%**
- Decision: **`{result['decision']}`**

The original unlabeled prediction file and SHA-256 remain unchanged; only the
protein-overlap classification was corrected.
"""
    (out / "REPORT.md").write_text(report)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
