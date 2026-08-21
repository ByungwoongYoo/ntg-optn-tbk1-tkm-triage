#!/usr/bin/env python3
"""Build and audit local current-database evidence for Panax A1/A2/B.

The two subcommands deliberately separate reference extraction from domain
parsing: curated proteins must be recovered and title-validated before Pfam is
run on the combined candidate/reference set.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
CONTROLS = ("PNX_Duplo_A_control", "PNX_Duplo_B_control")
BLAST_FIELDS = (
    "qseqid", "saccver", "pident", "length", "qlen", "slen", "qstart",
    "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "staxids",
    "sscinames", "stitle", "sseq",
)
UNIVEC_FIELDS = (
    "qseqid", "saccver", "pident", "length", "qlen", "qstart", "qend",
    "sstart", "send", "evalue", "score", "bitscore", "qseq", "sseq",
)


def sha_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_accession(header: str) -> str:
    first = header.split()[0]
    tokens = [first, *first.split("|")]
    pattern = re.compile(r"^(?:[A-Z]{1,4}|WP)_\d+(?:\.\d+)?$")
    for token in tokens:
        if token.startswith("PNX_"):
            return token
        if pattern.fullmatch(token):
            return token
    return first


def iter_fasta(path: Path) -> Iterable[tuple[str, str, str]]:
    """Yield FASTA records without loading a potentially large file at once."""
    name: str | None = None
    header: str | None = None
    chunks: list[str] = []
    yielded = 0
    with path.open(errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and header is not None:
                    sequence = "".join(chunks)
                    if not sequence:
                        raise SystemExit(f"empty FASTA record in {path}: {name}")
                    yield name, header, sequence
                    yielded += 1
                header = line[1:]
                name = canonical_accession(header)
                chunks = []
            elif name is None:
                raise SystemExit(f"sequence before first FASTA header in {path}")
            else:
                chunks.append(line)
    if name is not None:
        sequence = "".join(chunks)
        if not sequence:
            raise SystemExit(f"empty FASTA record in {path}: {name}")
        if header is None:
            raise SystemExit(f"missing FASTA header in {path}")
        yield name, header, sequence
        yielded += 1
    if not yielded:
        raise SystemExit(f"empty or malformed FASTA: {path}")


def read_fasta(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    seqs: dict[str, str] = {}
    headers: dict[str, str] = {}
    for name, header, sequence in iter_fasta(path):
        if name in headers:
            raise SystemExit(f"duplicate FASTA identifier in {path}: {name}")
        headers[name] = header
        seqs[name] = sequence
    return seqs, headers


def write_fasta(path: Path, records: Iterable[tuple[str, str, str]]) -> None:
    with path.open("w") as handle:
        for name, description, sequence in records:
            handle.write(f">{name} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else ["status"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_title(value: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.casefold()).split())


def extract_references(args: argparse.Namespace) -> int:
    manifest = read_tsv(args.manifest)
    required = {
        "accession", "context_group", "role", "expected_title_fragment",
        "expected_product_fragment", "expected_aa_length", "selection_basis", "selection_source",
    }
    if not manifest or not required.issubset(manifest[0]):
        raise SystemExit(f"curated manifest lacks required columns: {sorted(required)}")
    accessions = [row["accession"] for row in manifest]
    if len(accessions) != len(set(accessions)):
        raise SystemExit("duplicate curated reference accession")
    if sum(row["role"] == "rooting_sensitivity_reference" for row in manifest) != 2:
        raise SystemExit("curated reference manifest must contain exactly two rooting-sensitivity references")

    found_seq: dict[str, str] = {}
    found_header: dict[str, str] = {}
    wanted = set(accessions)
    for fasta in args.protein_fasta:
        for accession, header, sequence in iter_fasta(fasta):
            if accession not in wanted:
                continue
            if accession in found_seq:
                raise SystemExit(f"curated accession appears in multiple RefSeq files: {accession}")
            found_seq[accession] = sequence.upper()
            found_header[accession] = header
    missing = sorted(set(accessions) - set(found_seq))
    if missing:
        raise SystemExit(f"curated RefSeq accessions absent from current release: {missing}")

    args.out.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str, str]] = []
    provenance: list[dict[str, object]] = []
    for row in manifest:
        accession = row["accession"]
        sequence = found_seq[accession]
        header = found_header[accession]
        fragment = normalize_title(row["expected_title_fragment"])
        if fragment not in normalize_title(header):
            raise SystemExit(
                f"current RefSeq title mismatch for {accession}: expected fragment "
                f"{row['expected_title_fragment']!r}; observed {header!r}"
            )
        product_fragment = normalize_title(row["expected_product_fragment"])
        if product_fragment not in normalize_title(header):
            raise SystemExit(
                f"current RefSeq product mismatch for {accession}: expected fragment "
                f"{row['expected_product_fragment']!r}; observed {header!r}"
            )
        if not re.fullmatch(r"[ABCDEFGHIKJLMNPQRSTVWXYZ*]+", sequence):
            raise SystemExit(f"invalid amino-acid character in curated reference {accession}")
        sequence = sequence.rstrip("*")
        if "*" in sequence or len(sequence) < 300 or len(sequence) != int(row["expected_aa_length"]):
            raise SystemExit(f"curated reference is not a usable polyprotein/protein: {accession}")
        records.append((accession, f"context_group={row['context_group']}; role={row['role']}", sequence))
        provenance.append({
            **row, "observed_refseq_header": header, "protein_length": len(sequence),
            "sequence_sha256": sha_text(sequence),
        })
    write_fasta(args.out / "CURATED_REFERENCE_FULL.faa", records)
    write_tsv(args.out / "CURATED_REFERENCE_PROVENANCE.tsv", provenance)
    (args.out / "CURATED_REFERENCE_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return 0


def parse_blast(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(errors="replace", newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if not values:
                continue
            if len(values) != len(fields):
                raise SystemExit(f"malformed BLAST row in {path}: {len(values)} != {len(fields)}")
            rows.append(dict(zip(fields, values)))
    return rows


def parse_domtbl(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for raw in path.read_text(errors="replace").splitlines():
        if not raw or raw.startswith("#"):
            continue
        values = raw.split(maxsplit=22)
        if len(values) < 22:
            raise SystemExit(f"malformed HMMER domtblout row: {raw[:120]}")
        rows.append({
            "target_name": values[0], "target_accession": values[1], "model_length": int(values[2]),
            "query": values[3], "query_length": int(values[5]), "full_evalue": float(values[6]),
            "full_score": float(values[7]), "i_evalue": float(values[12]), "domain_score": float(values[13]),
            "hmm_from": int(values[15]), "hmm_to": int(values[16]), "ali_from": int(values[17]),
            "ali_to": int(values[18]), "env_from": int(values[19]), "env_to": int(values[20]),
            "accuracy": float(values[21]), "description": values[22] if len(values) > 22 else "",
        })
    return rows


def masked_fraction(original: str, masked: str) -> float:
    if len(original) != len(masked):
        raise SystemExit(f"masked sequence length mismatch: {len(original)} != {len(masked)}")
    changed = sum(
        (m.islower() or m.upper() in {"N", "X"}) and o.upper() not in {"N", "X"}
        for o, m in zip(original, masked)
    )
    return changed / len(original) if original else 1.0


def summarize_hits(rows: list[dict[str, str]], expected: set[str], label: str) -> dict[str, dict[str, object]]:
    unexpected = sorted({row["qseqid"] for row in rows} - expected)
    if unexpected:
        raise SystemExit(f"unexpected {label} query IDs: {unexpected}")
    out: dict[str, dict[str, object]] = {}
    for query in sorted(expected):
        hits = [row for row in rows if row["qseqid"] == query]
        hits.sort(key=lambda row: float(row["bitscore"]), reverse=True)
        top = hits[0] if hits else None
        out[query] = {
            "hit_count": len(hits),
            "near_identical_90pct_qcov80_count": sum(
                float(row["pident"]) >= 90 and float(row["qcovs"]) >= 80 for row in hits
            ),
            "near_identical_95pct_qcov80_count": sum(
                float(row["pident"]) >= 95 and float(row["qcovs"]) >= 80 for row in hits
            ),
            "top_hit": None if top is None else {key: top[key] for key in BLAST_FIELDS if key != "sseq"},
        }
    return out


def finalize(args: argparse.Namespace) -> int:
    args.out.mkdir(parents=True, exist_ok=True)
    candidate_orfs, _ = read_fasta(args.query_root / "panax_three_partial_orfs.faa")
    candidate_nt, _ = read_fasta(args.query_root / "panax_three_contigs.fna")
    control_orfs, _ = read_fasta(args.query_root / "panax_pipeline_controls_orfs.faa")
    if tuple(candidate_orfs) != CANDIDATES or set(control_orfs) != set(CONTROLS):
        raise SystemExit("exact candidate/control query sets were not preserved")
    references, _ = read_fasta(args.reference_full)
    reference_meta = {row["accession"]: row for row in read_tsv(args.reference_provenance)}
    if set(references) != set(reference_meta):
        raise SystemExit("curated reference FASTA/provenance mismatch")

    dom_rows = [row for row in parse_domtbl(args.domtbl)
                if str(row["target_accession"]).split(".")[0] == "PF00680"]
    by_query: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in dom_rows:
        by_query[str(row["query"])].append(row)
    expected_domain_ids = set(CANDIDATES) | set(references)
    unexpected_domains = sorted(set(by_query) - expected_domain_ids)
    if unexpected_domains:
        raise SystemExit(f"unexpected PF00680 query IDs: {unexpected_domains}")

    combined = {**candidate_orfs, **references}
    domain_gate: list[dict[str, object]] = []
    core_records: list[tuple[str, str, str]] = []
    domain_failures: list[str] = []
    for query in [*CANDIDATES, *references]:
        hits = sorted(by_query.get(query, []), key=lambda row: float(row["domain_score"]), reverse=True)
        if not hits:
            domain_failures.append(f"{query}:no_PF00680_hit")
            continue
        best = hits[0]
        sequence = combined[query]
        if int(best["query_length"]) != len(sequence):
            raise SystemExit(f"HMMER query length mismatch for {query}")
        start, end = int(best["ali_from"]), int(best["ali_to"])
        core = sequence[start - 1:end]
        model_coverage = (int(best["hmm_to"]) - int(best["hmm_from"]) + 1) / int(best["model_length"])
        minimum_coverage = 0.50 if query in CANDIDATES else 0.40
        passed = bool(
            float(best["i_evalue"]) <= 1e-3 and float(best["domain_score"]) >= 20
            and len(core) >= 100 and model_coverage >= minimum_coverage
        )
        if not passed:
            domain_failures.append(f"{query}:weak_or_partial_PF00680")
        meta = {
            "sequence_id": query, "role": "candidate" if query in CANDIDATES else reference_meta[query]["role"],
            "context_group": "Panax_candidate" if query in CANDIDATES else reference_meta[query]["context_group"],
            "PF00680_version": best["target_accession"], "full_sequence_length": len(sequence),
            "domain_ali_from": start, "domain_ali_to": end, "core_length": len(core),
            "hmm_from": best["hmm_from"], "hmm_to": best["hmm_to"],
            "model_length": best["model_length"], "model_coverage": f"{model_coverage:.6f}",
            "i_evalue": f"{float(best['i_evalue']):.6g}", "domain_score": f"{float(best['domain_score']):.3f}",
            "accuracy": f"{float(best['accuracy']):.4f}", "domain_gate_pass": str(passed).lower(),
            "full_sequence_sha256": sha_text(sequence), "core_sha256": sha_text(core),
        }
        domain_gate.append(meta)
        core_records.append((query, f"role={meta['role']}; context_group={meta['context_group']}; PF00680={start}-{end}", core))
    write_tsv(args.out / "DOMAIN_GATE.tsv", domain_gate)
    write_fasta(args.out / "RDRP_CORES.faa", core_records)

    protein_rows = parse_blast(args.refseq_blastp, BLAST_FIELDS)
    nucleotide_rows = parse_blast(args.refseq_blastn, BLAST_FIELDS)
    search_expected = set(CANDIDATES) | set(CONTROLS)
    protein_summary = summarize_hits(protein_rows, search_expected, "RefSeq protein")
    nucleotide_summary = summarize_hits(nucleotide_rows, search_expected, "RefSeq nucleotide")
    control_failures = [
        f"{control}:no_local_refseq_protein_hit" for control in CONTROLS
        if not protein_summary[control]["hit_count"]
    ] + [
        f"{control}:no_local_refseq_nucleotide_hit" for control in CONTROLS
        if not nucleotide_summary[control]["hit_count"]
    ]
    local_rows: list[dict[str, object]] = []
    for query in [*CANDIDATES, *CONTROLS]:
        p, n = protein_summary[query], nucleotide_summary[query]
        ptop, ntop = p["top_hit"], n["top_hit"]
        local_rows.append({
            "query": query, "role": "candidate" if query in CANDIDATES else "Durnavirales_like_pipeline_control",
            "protein_hit_count": p["hit_count"], "protein_top_accession": "" if ptop is None else ptop["saccver"],
            "protein_top_identity": "" if ptop is None else ptop["pident"],
            "protein_top_qcov": "" if ptop is None else ptop["qcovs"],
            "protein_top_evalue": "" if ptop is None else ptop["evalue"],
            "protein_top_title": "" if ptop is None else ptop["stitle"],
            "protein_near_identical_90pct_qcov80_count": p["near_identical_90pct_qcov80_count"],
            "nucleotide_hit_count": n["hit_count"], "nucleotide_top_accession": "" if ntop is None else ntop["saccver"],
            "nucleotide_top_identity": "" if ntop is None else ntop["pident"],
            "nucleotide_top_qcov": "" if ntop is None else ntop["qcovs"],
            "nucleotide_top_evalue": "" if ntop is None else ntop["evalue"],
            "nucleotide_near_identical_95pct_qcov80_count": n["near_identical_95pct_qcov80_count"],
        })
    write_tsv(args.out / "REFSEQ_LOCAL_SUMMARY.tsv", local_rows)

    univec_rows = parse_blast(args.univec, UNIVEC_FIELDS)
    unexpected_vector = sorted({row["qseqid"] for row in univec_rows} - set(CANDIDATES))
    if unexpected_vector:
        raise SystemExit(f"unexpected UniVec query IDs: {unexpected_vector}")
    original_dna = candidate_nt
    dust, _ = read_fasta(args.dust_masked)
    seg, _ = read_fasta(args.seg_masked)
    if set(dust) != set(CANDIDATES) or set(seg) != set(CANDIDATES):
        raise SystemExit("DUST/SEG masked query sets are incomplete")
    contamination_rows: list[dict[str, object]] = []
    contamination_failures: list[str] = []
    for query in CANDIDATES:
        hits = [row for row in univec_rows if row["qseqid"] == query]
        # NCBI VecScreen categories use raw BLAST score and different terminal
        # and internal thresholds. A hit is terminal when it touches the first
        # or last 25 query bases; adjacent-hit propagation is applied below.
        terminal = []
        for row in hits:
            lo, hi = sorted((int(row["qstart"]), int(row["qend"])))
            terminal.append(lo <= 25 or hi >= int(row["qlen"]) - 24)
        changed = True
        while changed:
            changed = False
            for i, row_i in enumerate(hits):
                if terminal[i]:
                    continue
                lo_i, hi_i = sorted((int(row_i["qstart"]), int(row_i["qend"])))
                for j, row_j in enumerate(hits):
                    if i == j or not terminal[j]:
                        continue
                    lo_j, hi_j = sorted((int(row_j["qstart"]), int(row_j["qend"])))
                    interval_gap = max(0, max(lo_i, lo_j) - min(hi_i, hi_j) - 1)
                    if interval_gap <= 25:
                        terminal[i] = True; changed = True; break
        categories = []
        for row, is_terminal in zip(hits, terminal):
            score = int(float(row["score"]))
            if (is_terminal and score >= 24) or (not is_terminal and score >= 30):
                category = "strong"
            elif (is_terminal and 19 <= score <= 23) or (not is_terminal and 25 <= score <= 29):
                category = "moderate"
            elif (is_terminal and 16 <= score <= 18) or (not is_terminal and 23 <= score <= 24):
                category = "weak"
            else:
                category = "below_weak"
            categories.append(category)
        strong = [row for row, category in zip(hits, categories) if category == "strong"]
        moderate = [row for row, category in zip(hits, categories) if category == "moderate"]
        weak = [row for row, category in zip(hits, categories) if category == "weak"]
        dust_fraction = masked_fraction(original_dna[query].upper(), dust[query])
        seg_fraction = masked_fraction(candidate_orfs[query].upper(), seg[query])
        passed = not strong and not moderate and dust_fraction < 0.50 and seg_fraction < 0.50
        if not passed:
            contamination_failures.append(f"{query}:vector_or_low_complexity_gate")
        contamination_rows.append({
            "candidate": query, "univec_raw_hit_count": len(hits), "univec_moderate_hit_count": len(moderate),
            "univec_strong_hit_count": len(strong), "univec_weak_hit_count": len(weak),
            "dust_masked_fraction": f"{dust_fraction:.6f}",
            "seg_masked_fraction": f"{seg_fraction:.6f}", "contamination_gate_pass": str(passed).lower(),
            "classification_note": "official NCBI VecScreen search preset and terminal/internal raw-score categories; moderate or strong produces PENDING/FAIL; raw alignments retained",
        })
    write_tsv(args.out / "CONTAMINATION_GATE.tsv", contamination_rows)

    duplicate_cores: list[list[str]] = []
    by_hash: dict[str, list[str]] = defaultdict(list)
    for name, _, sequence in core_records:
        by_hash[sha_text(sequence)].append(name)
    duplicate_cores = [names for names in by_hash.values() if len(names) > 1]
    technical_failures = domain_failures + control_failures + contamination_failures
    candidate_status = {}
    domain_by_id = {str(row["sequence_id"]): row for row in domain_gate}
    contamination_by_id = {str(row["candidate"]): row for row in contamination_rows}
    for candidate in CANDIDATES:
        p, n = protein_summary[candidate], nucleotide_summary[candidate]
        candidate_status[candidate] = {
            "PF00680_gate_pass": domain_by_id.get(candidate, {}).get("domain_gate_pass") == "true",
            "contamination_gate_pass": contamination_by_id[candidate]["contamination_gate_pass"] == "true",
            "current_refseq_viral_protein_hit_count": p["hit_count"],
            "current_refseq_viral_nucleotide_hit_count": n["hit_count"],
            "current_refseq_near_identical_protein": bool(p["near_identical_90pct_qcov80_count"]),
            "current_refseq_near_identical_nucleotide": bool(n["near_identical_95pct_qcov80_count"]),
        }
        candidate_status[candidate]["local_sequence_gate_pass"] = bool(
            candidate_status[candidate]["PF00680_gate_pass"]
            and candidate_status[candidate]["contamination_gate_pass"]
            and p["hit_count"]
        )
    technical_complete = not technical_failures and len(core_records) == len(expected_domain_ids)
    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "technical_complete": technical_complete,
        "failures": technical_failures,
        "candidate_status": candidate_status,
        "candidate_count": len(CANDIDATES), "curated_reference_count": len(references),
        "rooting_sensitivity_reference_count": sum(row["role"] == "rooting_sensitivity_reference" for row in reference_meta.values()),
        "exact_duplicate_core_groups": duplicate_cores,
        "database_provenance_sha256": sha_file(args.database_provenance),
        "interpretation_boundary": "Local matches and PF00680 support sequence-level Picornavirales-like candidates; they do not establish a new taxon, true host, active replication, or disease causation.",
    }
    (args.out / "LOCAL_EVIDENCE_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    report = [
        "# Panax A1/A2/B local sequence evidence", "",
        f"Technical gate: **{'PASS' if technical_complete else 'FAIL'}**", "",
        "All candidates were rebuilt from hash-locked source contigs. Current RefSeq viral files were selected from the official release catalog and MD5-checked; current Pfam-A was MD5-checked before local hmmscan; UniVec_Core was searched locally with the NCBI VecScreen BLAST preset.", "",
        "| Candidate | PF00680 | current RefSeq viral protein hits | near-identical protein | strong UniVec | local gate |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for candidate in CANDIDATES:
        c = candidate_status[candidate]
        report.append(
            f"| `{candidate}` | {c['PF00680_gate_pass']} | {c['current_refseq_viral_protein_hit_count']} | "
            f"{c['current_refseq_near_identical_protein']} | {contamination_by_id[candidate]['univec_strong_hit_count']} | "
            f"{c['local_sequence_gate_pass']} |"
        )
    report += ["", "This is a sequence audit of partial candidates, not a formal species or host assignment."]
    (args.out / "LOCAL_EVIDENCE_REPORT.md").write_text("\n".join(report) + "\n")
    for path in (args.reference_provenance, args.database_provenance):
        target = args.out / path.name
        if target.resolve() != path.resolve():
            target.write_bytes(path.read_bytes())
    sums = []
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            sums.append(f"{sha_file(path)}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n")
    print(json.dumps(status, indent=2))
    return 0 if technical_complete else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract-references")
    extract.add_argument("--manifest", type=Path, required=True)
    extract.add_argument("--protein-fasta", type=Path, action="append", required=True)
    extract.add_argument("--out", type=Path, required=True)
    extract.set_defaults(function=extract_references)
    final = sub.add_parser("finalize")
    final.add_argument("--query-root", type=Path, required=True)
    final.add_argument("--reference-full", type=Path, required=True)
    final.add_argument("--reference-provenance", type=Path, required=True)
    final.add_argument("--domtbl", type=Path, required=True)
    final.add_argument("--refseq-blastp", type=Path, required=True)
    final.add_argument("--refseq-blastn", type=Path, required=True)
    final.add_argument("--univec", type=Path, required=True)
    final.add_argument("--dust-masked", type=Path, required=True)
    final.add_argument("--seg-masked", type=Path, required=True)
    final.add_argument("--database-provenance", type=Path, required=True)
    final.add_argument("--out", type=Path, required=True)
    final.set_defaults(function=finalize)
    args = parser.parse_args()
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
