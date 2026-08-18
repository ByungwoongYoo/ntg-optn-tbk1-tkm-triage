#!/usr/bin/env python3
"""Aggregate accession-level candidate tables and rank validation priorities."""
from __future__ import annotations

import csv
import glob
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    tables = [Path(x) for x in glob.glob("downloads/**/candidate_summary.tsv", recursive=True)]
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for table in tables:
        with table.open(encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = fields or list(reader.fieldnames or [])
            rows.extend(reader)

    def f(row: dict[str, str], key: str, default: float) -> float:
        try:
            return float(row.get(key, ""))
        except ValueError:
            return default

    rows.sort(
        key=lambda r: (
            r.get("high_interest_screen") != "true",
            f(r, "pident", -1.0),
            -f(r, "palm_score", 0.0),
            -f(r, "aa_length", 0.0),
        )
    )

    out = Path("aggregate")
    out.mkdir(exist_ok=True)
    combined = out / "all_candidates.tsv"
    with combined.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        if fields:
            writer.writeheader()
            writer.writerows(rows)
        else:
            handle.write("accession\tquery\n")

    summaries: list[dict] = []
    for fn in glob.glob("downloads/**/summary.json", recursive=True):
        try:
            summaries.append(json.loads(Path(fn).read_text(encoding="utf-8")))
        except Exception:
            pass
    accessions_scanned = len(summaries)
    high = [r for r in rows if r.get("high_interest_screen") == "true"]
    novelty = Counter(r.get("novelty_screen", "unknown") for r in rows)

    report_lines = [
        "# Public-sequence RNA-virus hunt — first exhaustive pass",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Executive result",
        "",
        f"- Accessions with completed accession-level reports: **{accessions_scanned}**",
        f"- RdRP candidates passing palm_annot threshold: **{len(rows)}**",
        f"- High-interest divergence screen candidates: **{len(high)}**",
        "",
        "A high-interest screening label is **not** a discovery claim. It only prioritizes "
        "sequences for read support, replicate occurrence, contamination checks, current-database "
        "search, ORF architecture, and phylogenetic validation.",
        "",
        "## Novelty-screen distribution",
        "",
    ]
    if novelty:
        for key, count in sorted(novelty.items()):
            report_lines.append(f"- `{key}`: {count}")
    else:
        report_lines.append("- No RdRP candidate was returned in the completed accessions.")

    report_lines.extend(["", "## Highest-priority candidates", ""])
    if high:
        report_lines.extend([
            "| rank | accession | nominal organism | query | palm score | PALMdb identity | qcov | screen |",
            "|---:|---|---|---|---:|---:|---:|---|",
        ])
        for idx, row in enumerate(high[:50], 1):
            report_lines.append(
                f"| {idx} | `{row.get('accession','')}` | {row.get('organism','')} | "
                f"`{row.get('query','')}` | {row.get('palm_score') or 'NA'} | "
                f"{row.get('pident') or 'none'} | {row.get('qcov') or 'NA'} | "
                f"{row.get('novelty_screen','')} |"
            )
    else:
        report_lines.append("No candidate met the first-pass high-interest rule.")

    report_lines.extend([
        "",
        "## Required validation before any public claim",
        "",
        "1. Recover the exact nucleotide contig and verify an uninterrupted RdRP ORF and conserved motifs.",
        "2. Map original reads to show breadth, depth, paired-read support, and absence of assembly-only joins.",
        "3. Search current GenBank/RefSeq/UniProt and PALMdb; the 2023 PALMdb screen alone is insufficient.",
        "4. Test independent samples or related runs and inspect batch/index-hopping patterns.",
        "5. Exclude host endogenous viral elements, retroelements, vectors, and laboratory contaminants.",
        "6. Build a reference alignment and phylogeny; use cautious taxonomic wording.",
        "7. Treat the SRA-listed organism as sample metadata, not proof of the true biological host.",
        "",
        "## Claim boundary",
        "",
        "This pass can identify candidate RNA-virus polymerase lineages in public sequence data. "
        "It cannot by itself establish a new ICTV taxon, active infection, pathogenicity, or the true host.",
    ])
    report = out / "FIRST_PASS_REPORT.md"
    report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest_rows = []
    for path in sorted(p for p in out.rglob("*") if p.is_file()):
        manifest_rows.append({"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = out / "MANIFEST.json"
    manifest.write_text(json.dumps(manifest_rows, indent=2) + "\n", encoding="utf-8")
    # Recompute manifest list including the manifest itself in a separate outer hash file.
    (out / "MANIFEST_SHA256.txt").write_text(f"{sha256(manifest)}  {manifest.name}\n", encoding="utf-8")
    print(report.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
