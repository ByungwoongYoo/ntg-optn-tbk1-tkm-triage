#!/usr/bin/env python3
"""Conservative final adjudication of the Panax-associated RdRP candidates.

The strongest output is a *sequence-level divergent RdRP lineage candidate*.
The script never labels a new virus species, assigns Panax as the true host, or
infers infection, replication, pathogenicity, or treatment relevance.
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
from typing import Any

BLAST_FIELDS = [
    "qseqid", "saccver", "pident", "length", "qlen", "slen", "qcovs",
    "evalue", "bitscore", "staxids", "sscinames", "stitle",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_blast(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open(encoding="utf-8", errors="replace") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) >= len(BLAST_FIELDS):
                rows.append(dict(zip(BLAST_FIELDS, values[: len(BLAST_FIELDS)])))
    return rows


def num(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def best_by_query(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        q = row.get("qseqid", "")
        if not q:
            continue
        if q not in out or num(row.get("bitscore"), -1) > num(out[q].get("bitscore"), -1):
            out[q] = row
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def search_success(status: dict[str, Any], label: str) -> bool:
    try:
        return bool(status["searches"][label]["success"])
    except Exception:
        return False


def current_near_identical_protein(*hits: dict[str, str] | None) -> bool:
    for hit in hits:
        if hit and num(hit.get("pident"), 0) >= 90 and num(hit.get("qcovs"), 0) >= 80:
            return True
    return False


def current_near_identical_nt(hit: dict[str, str] | None) -> bool:
    return bool(hit and num(hit.get("pident"), 0) >= 95 and num(hit.get("qcovs"), 0) >= 80)


def nonviral_override(nonviral: dict[str, str] | None, viral: dict[str, str] | None) -> bool:
    if not nonviral:
        return False
    nbit = num(nonviral.get("bitscore"), 0)
    vbit = num((viral or {}).get("bitscore"), 0)
    npid = num(nonviral.get("pident"), 0)
    nqcov = num(nonviral.get("qcovs"), 0)
    ne = num(nonviral.get("evalue"), 1)
    if vbit > 0 and nbit > 1.20 * vbit and nqcov >= 40:
        return True
    if npid >= 70 and nqcov >= 70 and ne <= 1e-20:
        return True
    return False


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--members", type=Path, required=True)
    p.add_argument("--replication", type=Path, required=True)
    p.add_argument("--independence", type=Path, required=True)
    p.add_argument("--blast-status", type=Path, required=True)
    p.add_argument("--blastp-viral", type=Path, required=True)
    p.add_argument("--blastp-nonviral", type=Path, required=True)
    p.add_argument("--blastp-nr", type=Path, required=True)
    p.add_argument("--blastn-nt", type=Path, required=True)
    p.add_argument("--read-root", type=Path, required=True)
    p.add_argument("--phylogeny-status", type=Path, required=False)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, dict[str, Any]] = json.loads(args.manifest.read_text(encoding="utf-8"))
    members = read_tsv(args.members)
    member_by_lineage: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in members:
        member_by_lineage[row.get("lineage", "")].append(row)
    replication = {r.get("lineage", ""): r for r in read_tsv(args.replication)}
    independence = json.loads(args.independence.read_text(encoding="utf-8")) if args.independence.exists() else {}
    blast_status = json.loads(args.blast_status.read_text(encoding="utf-8")) if args.blast_status.exists() else {}

    viral_best = best_by_query(load_blast(args.blastp_viral))
    nonviral_best = best_by_query(load_blast(args.blastp_nonviral))
    nr_best = best_by_query(load_blast(args.blastp_nr))
    nt_best = best_by_query(load_blast(args.blastn_nt))
    phylogeny = {r.get("query", ""): r for r in read_tsv(args.phylogeny_status)} if args.phylogeny_status else {}

    read_rows: list[dict[str, str]] = []
    for path in sorted(args.read_root.rglob("read_support.tsv")):
        read_rows.extend(read_tsv(path))
    reads_by_ref: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_rows:
        reads_by_ref[row.get("reference", "")].append(row)

    required_searches = ["blastp_viral", "blastp_nonviral", "blastp_nr", "blastn_nt"]
    database_audit_complete = all(search_success(blast_status, x) for x in required_searches)
    independence_grade = independence.get("independence_grade", "metadata_missing")
    archived_sample_independence = independence_grade == "distinct_archived_samples_and_experiments"

    result_rows: list[dict[str, Any]] = []
    for lineage, meta in sorted(manifest.items()):
        rep = replication.get(lineage, {})
        lineage_reads = reads_by_ref.get(lineage, [])
        supported = [
            r for r in lineage_reads
            if num(r.get("breadth_1x"), 0) >= 0.90
            and num(r.get("mean_depth"), 0) >= 5.0
            and num(r.get("mapped_mapq20"), 0) >= 20
        ]
        supported_runs = sorted({r.get("run", "") for r in supported if r.get("run")})
        all_mapped_runs = sorted({r.get("run", "") for r in lineage_reads if r.get("run")})

        vb = viral_best.get(lineage)
        nb = nonviral_best.get(lineage)
        ab = nr_best.get(lineage)
        ntb = nt_best.get(lineage)
        phy = phylogeny.get(lineage, {})

        motif_complete = all(meta.get(k) for k in ("motif_A", "motif_B", "motif_C"))
        palm_hmm_support = bool(meta.get("hmm_rdrp_plus")) and bool(meta.get("nearest_palmdb"))
        assembly_replicated = rep.get("cross_run_replication") == "strong"
        raw_replicated = len(supported_runs) >= 2
        near_protein = current_near_identical_protein(vb, nb, ab)
        near_nt = current_near_identical_nt(ntb)
        cellular_override = nonviral_override(nb, vb)
        viral_remote_support = bool(vb and num(vb.get("evalue"), 1) <= 1e-5 and num(vb.get("qcovs"), 0) >= 25)
        phylogeny_built = phy.get("mafft") == "exit_0" and phy.get("fasttree") == "exit_0"

        if not database_audit_complete:
            decision = "database_audit_incomplete"
        elif near_protein or near_nt:
            decision = "current_database_near_identity_detected"
        elif cellular_override:
            decision = "ambiguous_or_nonviral_preferred_match"
        elif not (motif_complete and palm_hmm_support and assembly_replicated):
            decision = "insufficient_RdRP_or_replication_evidence"
        elif not raw_replicated:
            decision = "assembly_replicated_RdRP_candidate_raw_support_incomplete"
        elif archived_sample_independence:
            decision = "validated_divergent_RdRP_lineage_candidate"
        else:
            decision = "validated_within_project_RdRP_candidate_sample_independence_limited"

        member_rows = member_by_lineage.get(lineage, [])
        palm_identities = [num(x.get("pident"), math.nan) for x in member_rows]
        palm_identities = [x for x in palm_identities if not math.isnan(x)]

        result_rows.append({
            "lineage": lineage,
            "decision": decision,
            "classification_hint": meta.get("classification_hint", ""),
            "motif_A": meta.get("motif_A", ""),
            "motif_B": meta.get("motif_B", ""),
            "motif_C": meta.get("motif_C", ""),
            "pssm_score": meta.get("pssm_score", ""),
            "hmm_model": meta.get("hmm_rdrp_plus", ""),
            "nearest_PALMdb_2023": meta.get("nearest_palmdb", ""),
            "PALMdb_identity_min": min(palm_identities) if palm_identities else meta.get("pident", ""),
            "PALMdb_identity_max": max(palm_identities) if palm_identities else meta.get("pident", ""),
            "assembly_cross_run_replication": rep.get("cross_run_replication", "missing"),
            "assembly_distinct_runs": rep.get("distinct_runs", "0"),
            "raw_runs_evaluated": ";".join(all_mapped_runs),
            "raw_supported_run_count": len(supported_runs),
            "raw_supported_runs": ";".join(supported_runs),
            "independence_grade": independence_grade,
            "database_audit_complete": str(database_audit_complete).lower(),
            "viral_remote_support": str(viral_remote_support).lower(),
            "viral_hit_accession": (vb or {}).get("saccver", ""),
            "viral_hit_identity": (vb or {}).get("pident", ""),
            "viral_hit_qcov": (vb or {}).get("qcovs", ""),
            "viral_hit_evalue": (vb or {}).get("evalue", ""),
            "viral_hit_taxon": (vb or {}).get("sscinames", ""),
            "viral_hit_title": (vb or {}).get("stitle", ""),
            "nonviral_hit_accession": (nb or {}).get("saccver", ""),
            "nonviral_hit_identity": (nb or {}).get("pident", ""),
            "nonviral_hit_qcov": (nb or {}).get("qcovs", ""),
            "nonviral_hit_bitscore": (nb or {}).get("bitscore", ""),
            "nonviral_hit_taxon": (nb or {}).get("sscinames", ""),
            "nonviral_hit_title": (nb or {}).get("stitle", ""),
            "nr_top_accession": (ab or {}).get("saccver", ""),
            "nr_top_taxon": (ab or {}).get("sscinames", ""),
            "nt_hit_accession": (ntb or {}).get("saccver", ""),
            "nt_hit_identity": (ntb or {}).get("pident", ""),
            "nt_hit_qcov": (ntb or {}).get("qcovs", ""),
            "current_near_identical_protein": str(near_protein).lower(),
            "current_near_identical_nucleotide": str(near_nt).lower(),
            "nonviral_override": str(cellular_override).lower(),
            "phylogeny_built": str(phylogeny_built).lower(),
            "phylogeny_homologs_retrieved": phy.get("homologs_retrieved", ""),
        })

    fields = list(result_rows[0]) if result_rows else ["lineage", "decision"]
    with (args.out_dir / "FINAL_CANDIDATE_DECISIONS.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(result_rows)
    (args.out_dir / "FINAL_CANDIDATE_DECISIONS.json").write_text(
        json.dumps(result_rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    strict = [r for r in result_rows if r["decision"] == "validated_divergent_RdRP_lineage_candidate"]
    limited = [r for r in result_rows if r["decision"] == "validated_within_project_RdRP_candidate_sample_independence_limited"]
    known = [r for r in result_rows if r["decision"] == "current_database_near_identity_detected"]
    ambiguous = [r for r in result_rows if r["decision"] == "ambiguous_or_nonviral_preferred_match"]

    if strict:
        overall = "STRICT_SEQUENCE_LEVEL_CANDIDATES_SUPPORTED"
    elif limited:
        overall = "SEQUENCE_LEVEL_CANDIDATES_SUPPORTED_BUT_SAMPLE_INDEPENDENCE_LIMITED"
    else:
        overall = "NO_CANDIDATE_MET_THE_FULL_SEQUENCE_LEVEL_THRESHOLD"

    report = [
        "# Final Panax-associated RdRP validation report",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Overall result",
        "",
        f"**`{overall}`**",
        "",
        f"- Strict divergent RdRP lineage candidates: **{len(strict)}**",
        f"- Supported but archived-sample independence limited: **{len(limited)}**",
        f"- Current near-identity detections: **{len(known)}**",
        f"- Nonviral/ambiguous preferred matches: **{len(ambiguous)}**",
        f"- ENA independence grade: **`{independence_grade}`**",
        f"- Required current-database searches complete: **{database_audit_complete}**",
        "",
        "The strongest permitted interpretation is a sequence-level divergent RNA-dependent RNA polymerase (RdRP) lineage candidate in public Panax-associated transcriptomic data. This report does **not** establish a new virus species, active viral replication, infection of Panax, pathogenicity, transmissibility, or any medical effect.",
        "",
        "## Per-lineage evidence",
        "",
        "| Lineage | Decision | Assembly runs | Raw-supported runs | PALMdb identity range | Current viral hit | Protein id/qcov | NT id/qcov | Tree |",
        "|---|---|---:|---:|---:|---|---:|---:|---|",
    ]
    for row in result_rows:
        palmr = f"{row['PALMdb_identity_min']}-{row['PALMdb_identity_max']}%"
        report.append(
            f"| `{row['lineage']}` | `{row['decision']}` | {row['assembly_distinct_runs']} | "
            f"{row['raw_supported_run_count']} | {palmr} | "
            f"{row['viral_hit_taxon'] or row['viral_hit_accession'] or 'none returned'} | "
            f"{row['viral_hit_identity'] or 'NA'}/{row['viral_hit_qcov'] or 'NA'} | "
            f"{row['nt_hit_identity'] or 'NA'}/{row['nt_hit_qcov'] or 'NA'} | "
            f"{row['phylogeny_built']} |"
        )

    report += [
        "",
        "## Strict evidence gates",
        "",
        "A strict candidate required: complete A/B/C palm motifs; concordant PSSM/HMM and PALMdb support; strong cross-run assembly recurrence; raw-read coverage of at least 90% at mean depth at least 5× with at least 20 MAPQ≥20 reads in at least two runs; successful current viral, nonviral, unrestricted-protein, and nucleotide database searches; no protein hit at ≥90% identity over ≥80% of the query; no nucleotide hit at ≥95% identity over ≥80% of the query; and no substantially stronger nonviral explanation.",
        "",
        "## Independence boundary",
        "",
        f"The ENA audit grade is `{independence_grade}`. Distinct archive accessions do not by themselves prove independent plants, independent infections, or the true host.",
        "",
        "## What remains necessary for a formal virus discovery",
        "",
        "Formal taxonomic or host claims would require expert viral-taxonomy review, complete or substantially extended genome architecture, stronger contamination/index-hopping exclusion, and ideally independent biological sampling with strand-aware or targeted validation.",
    ]
    (args.out_dir / "FINAL_VALIDATION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    boundaries = [
        "# Claim boundary",
        "",
        "Permitted if supported by the final table:",
        "- divergent RdRP sequence-lineage candidate detected in public Panax-associated RNA-seq data;",
        "- candidate recurs across specified archived runs;",
        "- candidate is supported by raw-read mapping and current database audit.",
        "",
        "Not permitted from these data alone:",
        "- a formally named new virus species or genus;",
        "- Panax notoginseng is the true host;",
        "- active infection or viral replication;",
        "- disease causation, agricultural risk, human infection, or medical relevance;",
        "- absence from every unpublished or unindexed dataset.",
    ]
    (args.out_dir / "CLAIM_BOUNDARY.md").write_text("\n".join(boundaries) + "\n", encoding="utf-8")

    if read_rows:
        read_fields = list(read_rows[0])
        with (args.out_dir / "COMBINED_READ_SUPPORT.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=read_fields, delimiter="\t")
            writer.writeheader(); writer.writerows(read_rows)

    # Write manifest after all principal outputs. The manifest excludes itself to
    # avoid self-reference and is accompanied by its own hash.
    manifest_rows = []
    for path in sorted(p for p in args.out_dir.rglob("*") if p.is_file() and p.name not in {"MANIFEST.json", "MANIFEST_SHA256.txt"}):
        manifest_rows.append({
            "path": str(path.relative_to(args.out_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    (args.out_dir / "MANIFEST.json").write_text(json.dumps(manifest_rows, indent=2) + "\n", encoding="utf-8")
    (args.out_dir / "MANIFEST_SHA256.txt").write_text(
        f"{sha256(args.out_dir / 'MANIFEST.json')}  MANIFEST.json\n", encoding="utf-8"
    )
    print((args.out_dir / "FINAL_VALIDATION_REPORT.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
