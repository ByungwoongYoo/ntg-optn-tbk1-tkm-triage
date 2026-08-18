#!/usr/bin/env python3
"""Combine palm_annot and PALMdb similarity results into an audit table.

Novelty labels are screening labels, not taxonomic conclusions. In particular,
sequence divergence alone does not prove a new virus species or host infection.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterator


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return
    label = None
    parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if label is not None:
                    yield label, "".join(parts)
                label = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if label is not None:
            yield label, "".join(parts)


def parse_fev(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            cols = raw.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            d: dict[str, str] = {}
            for item in cols[1:]:
                if "=" in item:
                    key, value = item.split("=", 1)
                    d[key] = value
            out[cols[0]] = d
    return out


def as_float(value: str | None) -> float | None:
    try:
        return float(value) if value not in {None, ""} else None
    except ValueError:
        return None


def novelty_class(pident: float | None, qcov: float | None) -> str:
    if pident is None:
        return "no_PALMdb_hit"
    if qcov is not None and qcov < 40:
        return "weak_partial_hit"
    if pident < 35:
        return "very_divergent_lt35pct"
    if pident < 50:
        return "divergent_35_50pct"
    if pident < 70:
        return "moderate_50_70pct"
    if pident < 90:
        return "close_70_90pct"
    return "PALMdb_sOTU_level_ge90pct"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--accession", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--bioproject", default="")
    p.add_argument("--organism", default="")
    p.add_argument("--fev", type=Path, required=True)
    p.add_argument("--rdrp", type=Path, required=True)
    p.add_argument("--fullnt", type=Path, required=True)
    p.add_argument("--diamond", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--json", type=Path, required=True)
    p.add_argument("--markdown", type=Path, required=True)
    args = p.parse_args()

    fev = parse_fev(args.fev)
    aa = dict(read_fasta(args.rdrp) or [])
    nt = dict(read_fasta(args.fullnt) or [])

    best: dict[str, dict[str, str]] = {}
    if args.diamond.exists():
        with args.diamond.open(encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 12:
                    continue
                q = row[0]
                if q in best:
                    continue
                best[q] = {
                    "subject": row[1],
                    "pident": row[2],
                    "aln_len": row[3],
                    "qlen": row[4],
                    "slen": row[5],
                    "evalue": row[6],
                    "bitscore": row[7],
                    "qcov": row[8],
                    "scov": row[9],
                    "qstart": row[10],
                    "qend": row[11],
                }

    fieldnames = [
        "accession", "target_label", "bioproject", "organism", "query", "nt_contig",
        "nt_length", "aa_length", "palm_score", "pssm_score", "gene", "motif_order",
        "confidence", "group_guess", "palmprint_start", "palmprint_end", "motifs",
        "best_palmdb_sotu", "pident", "qcov", "scov", "evalue", "bitscore",
        "novelty_screen", "high_interest_screen",
    ]
    rows: list[dict[str, str]] = []
    for query, seq in aa.items():
        annotation = fev.get(query, {})
        hit = best.get(query, {})
        pident = as_float(hit.get("pident"))
        qcov = as_float(hit.get("qcov"))
        palm_score = as_float(
            annotation.get("score")
            or annotation.get("palm_score")
            or annotation.get("pa_score")
            or annotation.get("rdrp_score")
        )
        base = query.split("_frame=")[0]
        nt_seq = nt.get(base, "")
        novelty = novelty_class(pident, qcov)
        high_interest = (
            (palm_score is None or palm_score >= 75)
            and (pident is None or pident < 50)
            and (qcov is None or qcov >= 40)
            and len(seq) >= 90
        )
        rows.append(
            {
                "accession": args.accession,
                "target_label": args.label,
                "bioproject": args.bioproject,
                "organism": args.organism,
                "query": query,
                "nt_contig": base,
                "nt_length": str(len(nt_seq)),
                "aa_length": str(len(seq)),
                "palm_score": "" if palm_score is None else f"{palm_score:.3f}",
                "pssm_score": annotation.get("pssm_score", annotation.get("pssm_total_score", "")),
                "gene": annotation.get("gene", ""),
                "motif_order": annotation.get("order", ""),
                "confidence": annotation.get("confidence", ""),
                "group_guess": annotation.get("group", ""),
                "palmprint_start": annotation.get("pp_start", ""),
                "palmprint_end": annotation.get("pp_end", ""),
                "motifs": annotation.get("motifs", ""),
                "best_palmdb_sotu": hit.get("subject", ""),
                "pident": hit.get("pident", ""),
                "qcov": hit.get("qcov", ""),
                "scov": hit.get("scov", ""),
                "evalue": hit.get("evalue", ""),
                "bitscore": hit.get("bitscore", ""),
                "novelty_screen": novelty,
                "high_interest_screen": "true" if high_interest else "false",
            }
        )

    rows.sort(
        key=lambda r: (
            r["high_interest_screen"] != "true",
            float(r["pident"]) if r["pident"] else -1.0,
            -int(r["aa_length"] or 0),
        )
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "accession": args.accession,
        "target_label": args.label,
        "bioproject": args.bioproject,
        "organism": args.organism,
        "rdrp_candidates": len(rows),
        "high_interest_screen": sum(r["high_interest_screen"] == "true" for r in rows),
        "novelty_counts": {},
        "claim_boundary": (
            "These are computational RdRP/palmprint candidates. A candidate is not a confirmed "
            "new virus, does not prove infection of the nominal host, and requires independent "
            "sequence, contamination, replication, and taxonomic validation."
        ),
    }
    for row in rows:
        key = row["novelty_screen"]
        summary["novelty_counts"][key] = summary["novelty_counts"].get(key, 0) + 1
    args.json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        f"# RdRP screen: {args.accession}",
        "",
        f"- Target group: `{args.label}`",
        f"- BioProject: `{args.bioproject}`",
        f"- Nominal organism: `{args.organism}`",
        f"- RdRP candidates: **{len(rows)}**",
        f"- High-interest screening candidates: **{summary['high_interest_screen']}**",
        "",
        "## Claim boundary",
        "",
        summary["claim_boundary"],
        "",
        "## Candidate preview",
        "",
    ]
    if not rows:
        md.append("No candidate passed the configured palm_annot threshold in this accession.")
    else:
        md.extend([
            "| query | aa len | palm score | best PALMdb identity | qcov | screen label | high interest |",
            "|---|---:|---:|---:|---:|---|---|",
        ])
        for row in rows[:20]:
            md.append(
                f"| `{row['query']}` | {row['aa_length']} | {row['palm_score'] or 'NA'} | "
                f"{row['pident'] or 'none'} | {row['qcov'] or 'NA'} | {row['novelty_screen']} | "
                f"{row['high_interest_screen']} |"
            )
    args.markdown.write_text("\n".join(md) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
