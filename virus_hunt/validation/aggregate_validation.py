#!/usr/bin/env python3
"""Aggregate strict replication, raw-read support, and current-database searches.

Decision labels are deliberately below the level of a formal virus discovery. The
strongest positive label is `supported_unreported_RdRP_lineage_candidate`.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(row: dict[str, str] | None, key: str, default: float = math.nan) -> float:
    if not row:
        return default
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return default


def best_by_query(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        q = row.get("qseqid", "")
        if q and q not in out:
            out[q] = row
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--replication", type=Path, required=True)
    p.add_argument("--blastp-viral", type=Path, required=True)
    p.add_argument("--blastp-all", type=Path, required=True)
    p.add_argument("--blastn-nt", type=Path, required=True)
    p.add_argument("--read-root", type=Path, required=True)
    p.add_argument("--phylogeny-status", type=Path, required=False)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, dict] = json.loads(args.manifest.read_text())
    replication = {r["lineage"]: r for r in read_tsv(args.replication)}

    blast_fields = [
        "qseqid", "saccver", "pident", "length", "qlen", "slen", "qcovs",
        "evalue", "bitscore", "staxids", "sscinames", "stitle",
    ]
    def load_blast(path: Path) -> list[dict[str, str]]:
        rows=[]
        if not path.exists(): return rows
        with path.open(encoding="utf-8", errors="replace") as handle:
            reader=csv.reader(handle,delimiter="\t")
            for x in reader:
                if len(x) >= len(blast_fields):
                    rows.append(dict(zip(blast_fields, x[:len(blast_fields)])))
        return rows

    viral_rows = load_blast(args.blastp_viral)
    all_rows = load_blast(args.blastp_all)
    nt_rows = load_blast(args.blastn_nt)
    viral_best = best_by_query(viral_rows)
    all_best = best_by_query(all_rows)
    nt_best = best_by_query(nt_rows)

    read_rows: list[dict[str, str]] = []
    for path in args.read_root.rglob("read_support.tsv"):
        read_rows.extend(read_tsv(path))
    reads_by_ref: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_rows:
        reads_by_ref[row.get("reference", "")].append(row)

    result_rows=[]
    for lineage, meta in manifest.items():
        rep = replication.get(lineage, {})
        rb = reads_by_ref.get(lineage, [])
        supported_runs = [
            r for r in rb
            if f(r, "breadth_1x", 0) >= 0.90
            and f(r, "mean_depth", 0) >= 5.0
            and f(r, "mapped_mapq20", 0) >= 20
        ]
        vb = viral_best.get(lineage)
        ab = all_best.get(lineage)
        nb = nt_best.get(lineage)
        protein_current_known = (
            (f(vb, "pident", 0) >= 90 and f(vb, "qcovs", 0) >= 80)
            or (f(ab, "pident", 0) >= 90 and f(ab, "qcovs", 0) >= 80)
        )
        nucleotide_current_known = f(nb, "pident", 0) >= 95 and f(nb, "qcovs", 0) >= 80
        unrestricted_overrides_viral = False
        if ab and vb:
            unrestricted_overrides_viral = (
                f(ab, "bitscore", 0) > 1.20 * f(vb, "bitscore", 0)
                and ab.get("saccver") != vb.get("saccver")
                and "virus" not in (ab.get("sscinames", "") + " " + ab.get("stitle", "")).lower()
            )
        elif ab and not vb:
            unrestricted_overrides_viral = "virus" not in (ab.get("sscinames", "") + " " + ab.get("stitle", "")).lower()

        motif_complete = all(meta.get(k) for k in ["motif_A", "motif_B", "motif_C"])
        assembly_replicated = rep.get("cross_run_replication") == "strong"
        raw_replicated = len({r.get("run") for r in supported_runs}) >= 2
        current_viral_support = vb is not None and f(vb, "evalue", 1) <= 1e-5

        if protein_current_known or nucleotide_current_known:
            decision = "currently_known_or_near_identical_sequence"
        elif unrestricted_overrides_viral:
            decision = "rejected_or_ambiguous_cellular_match"
        elif motif_complete and assembly_replicated and raw_replicated and current_viral_support:
            decision = "supported_unreported_RdRP_lineage_candidate"
        elif motif_complete and assembly_replicated and current_viral_support:
            decision = "assembly_replicated_RdRP_candidate_raw_support_incomplete"
        else:
            decision = "insufficient_validation"

        result_rows.append({
            "lineage": lineage,
            "decision": decision,
            "assembly_cross_run_replication": rep.get("cross_run_replication", "missing"),
            "assembly_distinct_runs": rep.get("distinct_runs", "0"),
            "raw_supported_run_count": str(len({r.get('run') for r in supported_runs})),
            "raw_supported_runs": ";".join(sorted({r.get('run','') for r in supported_runs})),
            "minimum_PALMdb_2023_identity": rep.get("minimum_palmdb_identity", str(meta.get("pident", ""))),
            "viral_hit_accession": (vb or {}).get("saccver", ""),
            "viral_hit_identity": (vb or {}).get("pident", ""),
            "viral_hit_qcov": (vb or {}).get("qcovs", ""),
            "viral_hit_taxon": (vb or {}).get("sscinames", ""),
            "viral_hit_title": (vb or {}).get("stitle", ""),
            "unrestricted_hit_accession": (ab or {}).get("saccver", ""),
            "unrestricted_hit_identity": (ab or {}).get("pident", ""),
            "unrestricted_hit_qcov": (ab or {}).get("qcovs", ""),
            "unrestricted_hit_taxon": (ab or {}).get("sscinames", ""),
            "unrestricted_hit_title": (ab or {}).get("stitle", ""),
            "nt_hit_accession": (nb or {}).get("saccver", ""),
            "nt_hit_identity": (nb or {}).get("pident", ""),
            "nt_hit_qcov": (nb or {}).get("qcovs", ""),
            "current_near_identical_protein": str(protein_current_known).lower(),
            "current_near_identical_nucleotide": str(nucleotide_current_known).lower(),
        })

    fields=list(result_rows[0]) if result_rows else ["lineage","decision"]
    with (args.out_dir/"FINAL_CANDIDATE_DECISIONS.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t");writer.writeheader();writer.writerows(result_rows)
    (args.out_dir/"FINAL_CANDIDATE_DECISIONS.json").write_text(json.dumps(result_rows,indent=2)+"\n")

    positive=[r for r in result_rows if r["decision"]=="supported_unreported_RdRP_lineage_candidate"]
    report=[
        "# Public-sequence Panax RNA-virus-lineage validation",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Decision",
        "",
        f"Strictly supported, currently unreported RdRP-lineage candidates: **{len(positive)}**.",
        "",
        "This wording does **not** assert a new virus species, active infection, pathogenicity, or Panax as the true host.",
        "",
        "## Per-lineage audit",
        "",
        "| Lineage | Decision | assembly runs | raw-read runs | current viral hit | protein identity/qcov | current nt identity/qcov |",
        "|---|---|---:|---:|---|---:|---:|",
    ]
    for r in result_rows:
        report.append(
            f"| `{r['lineage']}` | `{r['decision']}` | {r['assembly_distinct_runs']} | "
            f"{r['raw_supported_run_count']} | {r['viral_hit_taxon'] or r['viral_hit_accession'] or 'none'} | "
            f"{r['viral_hit_identity'] or 'NA'}/{r['viral_hit_qcov'] or 'NA'} | "
            f"{r['nt_hit_identity'] or 'NA'}/{r['nt_hit_qcov'] or 'NA'} |"
        )
    report += [
        "",
        "## Evidence standard used",
        "",
        "A positive screening decision required all of the following:",
        "",
        "1. complete A/B/C RdRP palm-motif evidence;",
        "2. close sequence recurrence across at least two independently archived runs;",
        "3. raw-read mapping across at least 90% of the representative contig at mean depth at least 5× in at least two runs;",
        "4. a significant current viral-protein database match;",
        "5. no current protein match at ≥90% identity over ≥80% of the query and no current nucleotide match at ≥95% identity over ≥80%;",
        "6. no substantially stronger unrestricted cellular-protein match.",
        "",
        "## Remaining biological uncertainty",
        "",
        "Even a positive decision remains a sequence-level candidate. Formal virus discovery and host assignment require additional taxonomic review, contamination analysis, and ideally independent biological sampling or targeted validation.",
    ]
    (args.out_dir/"FINAL_VALIDATION_REPORT.md").write_text("\n".join(report)+"\n")

    if read_rows:
        read_fields=list(read_rows[0])
        with (args.out_dir/"COMBINED_READ_SUPPORT.tsv").open("w",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=read_fields,delimiter="\t");writer.writeheader();writer.writerows(read_rows)

    manifest_rows=[]
    for path in sorted(p for p in args.out_dir.rglob("*") if p.is_file()):
        manifest_rows.append({"path":str(path.relative_to(args.out_dir)),"bytes":path.stat().st_size,"sha256":sha256(path)})
    (args.out_dir/"MANIFEST.json").write_text(json.dumps(manifest_rows,indent=2)+"\n")
    (args.out_dir/"MANIFEST_SHA256.txt").write_text(f"{sha256(args.out_dir/'MANIFEST.json')}  MANIFEST.json\n")
    print((args.out_dir/"FINAL_VALIDATION_REPORT.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
