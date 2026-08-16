#!/usr/bin/env python3
"""Extract one canonical protein sequence per ProteinGym clinical CSV for homology clustering."""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import pandas as pd


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.strip().upper().encode()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical-score-zip", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    failures = []
    fasta_lines = []
    with zipfile.ZipFile(args.clinical_score_zip) as z:
        members = sorted(n for n in z.namelist() if n.lower().endswith(".csv"))
        for i, member in enumerate(members, 1):
            try:
                with z.open(member) as f:
                    d = pd.read_csv(f, usecols=lambda c: c in {"protein_sequence", "mutant", "DMS_bin_score"}, low_memory=False)
                if "protein_sequence" not in d.columns:
                    failures.append({"file": member, "reason": "missing protein_sequence"})
                    continue
                seqs = d["protein_sequence"].dropna().astype(str).str.strip().str.upper().unique()
                if len(seqs) != 1:
                    failures.append({"file": member, "reason": f"sequence_count={len(seqs)}"})
                    continue
                seq = seqs[0]
                if not seq or any(ch not in "ACDEFGHIKLMNPQRSTVWYUXOBZJ" for ch in seq):
                    failures.append({"file": member, "reason": "invalid sequence characters"})
                    continue
                protein_file = Path(member).name
                header = protein_file.removesuffix(".csv")
                cluster_seq = "".join(ch if ch in "ACDEFGHIKLMNPQRSTVWY" else "X" for ch in seq)
                h = seq_hash(seq)
                rows.append({
                    "protein_file": protein_file,
                    "header": header,
                    "sequence_hash": h,
                    "length": len(seq),
                    "zip_member": member,
                })
                fasta_lines.extend([f">{header}", cluster_seq])
            except Exception as exc:
                failures.append({"file": member, "reason": repr(exc)})
            if i % 500 == 0:
                print(f"prepare {i}/{len(members)} accepted={len(rows)}", flush=True)

    pd.DataFrame(rows).sort_values("protein_file").to_csv(out / "sequence_inventory.csv", index=False)
    (out / "clinical_sequences.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    result = {
        "members_scanned": len(members),
        "sequences_written": len(rows),
        "unique_sequence_hashes": len({r["sequence_hash"] for r in rows}),
        "failures": failures,
    }
    (out / "prepare_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
