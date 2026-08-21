#!/usr/bin/env python3
"""Aggregate exact A1/A2/B sequence-audit artifacts and fail closed."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
CONTROLS = ("PNX_Duplo_A_control", "PNX_Duplo_B_control")
MODE_CONTROLS = {
    "protein_viral": CONTROLS,
    "nt_viral": CONTROLS,
    "nt_megablast": CONTROLS,
    "protein_cellular": ("PNX_Panax_L2_control",),
    "nt_cellular": ("PNX_Panax_cpDNA_control",),
    "nt_panax": ("PNX_Panax_cpDNA_control",),
}
REMOTE_MODES = (
    "protein_viral", "protein_cellular", "protein_tsa", "protein_environmental",
    "nt_viral", "nt_cellular", "nt_megablast", "nt_panax", "nt_tsa",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(root: Path) -> list[str]:
    manifest = root / "SHA256SUMS.txt"
    failures: list[str] = []
    if not manifest.is_file():
        return [f"missing_checksum_manifest:{root.name}"]
    seen = set()
    for number, line in enumerate(manifest.read_text(errors="replace").splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            failures.append(f"malformed_checksum_row:{root.name}:{number}")
            continue
        expected, raw_relative = fields[0], fields[1].lstrip(" *")
        relative_path = Path(raw_relative)
        if relative_path.is_absolute():
            failures.append(f"unsafe_checksum_path:{root.name}:{raw_relative}")
            continue
        path = (root / relative_path).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError:
            failures.append(f"unsafe_checksum_path:{root.name}:{raw_relative}")
            continue
        normalized = path.relative_to(root.resolve()).as_posix()
        if normalized == "SHA256SUMS.txt":
            failures.append(f"checksum_manifest_lists_itself:{root.name}")
            continue
        if normalized in seen:
            failures.append(f"duplicate_checksum_path:{root.name}:{normalized}")
        seen.add(normalized)
        if not path.is_file():
            failures.append(f"missing_checksummed_file:{root.name}:{normalized}")
        elif sha(path) != expected:
            failures.append(f"checksum_mismatch:{root.name}:{normalized}")
    if not seen:
        failures.append(f"empty_checksum_manifest:{root.name}")
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.resolve() != manifest.resolve()
    }
    for relative in sorted(actual - seen):
        failures.append(f"unchecksummed_file:{root.name}:{relative}")
    return failures


def load_json(path: Path, failures: list[str], label: str, object_only: bool = True):
    if not path.is_file():
        failures.append(f"missing_json:{label}")
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception as exc:
        failures.append(f"invalid_json:{label}:{type(exc).__name__}")
        return {}
    if object_only and not isinstance(value, dict):
        failures.append(f"nonobject_json:{label}")
        return {}
    return value


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fields = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collected", type=Path, required=True)
    parser.add_argument("--upstream-results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")

    failures: list[str] = []
    expected_artifacts = [
        "panax-query-preflight", "panax-local-evidence", "panax-rdrp-phylogeny",
        *[f"panax-remote-{mode}" for mode in REMOTE_MODES],
    ]
    observed_artifacts = sorted(path.name for path in args.collected.iterdir() if path.is_dir()) if args.collected.is_dir() else []
    missing_artifacts = sorted(set(expected_artifacts) - set(observed_artifacts))
    unexpected_artifacts = sorted(set(observed_artifacts) - set(expected_artifacts))
    failures.extend(f"missing_artifact:{name}" for name in missing_artifacts)
    failures.extend(f"unexpected_artifact:{name}" for name in unexpected_artifacts)
    for name in expected_artifacts:
        root = args.collected / name
        if root.is_dir():
            failures.extend(verify_manifest(root))

    upstream = load_json(args.upstream_results, failures, "UPSTREAM_RESULTS")
    required_jobs = {"preflight", "local_evidence", "remote_search", "phylogeny"}
    if set(upstream) != required_jobs:
        failures.append(f"upstream_result_key_mismatch:{sorted(upstream)}")
    for job in sorted(required_jobs):
        if upstream.get(job) != "success":
            failures.append(f"upstream_job_not_success:{job}:{upstream.get(job,'missing')}")

    preflight = args.collected / "panax-query-preflight"
    query_manifest = load_json(preflight / "QUERY_MANIFEST.json", failures, "QUERY_MANIFEST", object_only=False)
    observed_candidates = [row.get("candidate") for row in query_manifest] if isinstance(query_manifest, list) else []
    if tuple(observed_candidates) != CANDIDATES:
        failures.append(f"candidate_query_set_mismatch:{observed_candidates}")
    control_manifest_path = preflight / "PIPELINE_CONTROL_MANIFEST.json"
    try:
        control_manifest = json.loads(control_manifest_path.read_text()) if control_manifest_path.is_file() else []
    except Exception:
        control_manifest = []
    if tuple(row.get("control") for row in control_manifest) != CONTROLS:
        failures.append("pipeline_control_set_mismatch")

    local = load_json(
        args.collected / "panax-local-evidence" / "LOCAL_EVIDENCE_STATUS.json",
        failures, "LOCAL_EVIDENCE_STATUS",
    )
    if not local.get("technical_complete"):
        failures.append("local_evidence_not_complete")
    tree = load_json(
        args.collected / "panax-rdrp-phylogeny" / "TREE_QC.json",
        failures, "TREE_QC",
    )
    if not tree.get("technical_complete"):
        failures.append("phylogeny_not_complete")

    remote_rows: list[dict[str, object]] = []
    remote_status: dict[str, dict] = {}
    remote_complete: dict[str, bool] = {}
    for mode in REMOTE_MODES:
        status = load_json(
            args.collected / f"panax-remote-{mode}" / "SEARCH_STATUS.json",
            failures, f"SEARCH_STATUS:{mode}",
        )
        remote_status[mode] = status
        expected_controls = set(MODE_CONTROLS.get(mode, ()))
        expected_ids = set(CANDIDATES) | expected_controls
        observed_ids = set(status.get("query_ids", []))
        observed_controls = set(status.get("validation_control_ids", []))
        control_results = status.get("validation_control_results", {})
        controls_valid = bool(
            set(control_results) == expected_controls
            and all(control_results.get(control, {}).get("validated") for control in expected_controls)
        ) if expected_controls else control_results in ({}, None)
        per_query = status.get("per_query", {})
        observed_per_query = set(per_query) if isinstance(per_query, dict) else set()
        complete = bool(
            status.get("technical_complete") and status.get("mode") == mode
            and observed_ids == expected_ids and observed_per_query == expected_ids
            and observed_controls == expected_controls
            and controls_valid
        )
        remote_complete[mode] = complete
        if not complete:
            failures.append(f"remote_search_not_complete:{mode}")
        for query in sorted(expected_ids):
            details = status.get("per_query", {}).get(query, {})
            remote_rows.append({
                "mode": mode, "database": status.get("database", ""), "query": query,
                "technical_complete": str(complete).lower(), "hit_count": details.get("hit_count", ""),
                "near_identical_qcov80_pident90_count": details.get("near_identical_qcov80_pident90_count", ""),
                "near_identical_qcov80_pident95_count": details.get("near_identical_qcov80_pident95_count", ""),
                "top_accession": (details.get("top_hit") or {}).get("saccver", ""),
                "top_identity": (details.get("top_hit") or {}).get("pident", ""),
                "top_qcov": (details.get("top_hit") or {}).get("qcovs", ""),
                "top_evalue": (details.get("top_hit") or {}).get("evalue", ""),
                "top_title": (details.get("top_hit") or {}).get("stitle", ""),
            })
        if expected_controls and complete:
            for control in sorted(expected_controls):
                if not status.get("per_query", {}).get(control, {}).get("hit_count"):
                    failures.append(f"pipeline_control_no_hit:{mode}:{control}")
    write_tsv(args.out / "REMOTE_COMPLETENESS.tsv", remote_rows)

    local_candidates = local.get("candidate_status", {})
    gate_rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        lc = local_candidates.get(candidate, {})
        protein_near = bool(lc.get("current_refseq_near_identical_protein"))
        nucleotide_near = bool(lc.get("current_refseq_near_identical_nucleotide"))
        remote_all_complete = True
        for mode, status in remote_status.items():
            details = status.get("per_query", {}).get(candidate, {})
            if not remote_complete.get(mode, False):
                remote_all_complete = False
            if mode.startswith("protein_") and details.get("near_identical_qcov80_pident90_count", 0):
                protein_near = True
            if mode.startswith("nt_") and details.get("near_identical_qcov80_pident95_count", 0):
                nucleotide_near = True
        local_pass = bool(local.get("technical_complete") and lc.get("local_sequence_gate_pass"))
        phylogeny_complete = bool(tree.get("technical_complete"))
        if not local_pass:
            decision = "local_sequence_gate_failed_or_pending"
        elif not phylogeny_complete:
            decision = "technical_phylogeny_incomplete"
        elif protein_near or nucleotide_near:
            decision = "currently_known_or_near_identical_sequence"
        elif remote_all_complete:
            decision = "divergent_Picornavirales_like_partial_sequence_candidate"
        else:
            decision = "technical_search_incomplete"
        gate_rows.append({
            "candidate": candidate, "exact_query_preflight": str(candidate in observed_candidates).lower(),
            "PF00680_gate": str(bool(lc.get("PF00680_gate_pass"))).lower(),
            "contamination_gate": str(bool(lc.get("contamination_gate_pass"))).lower(),
            "current_refseq_viral_protein_hits": lc.get("current_refseq_viral_protein_hit_count", ""),
            "near_identical_protein_detected": str(protein_near).lower(),
            "near_identical_nucleotide_detected": str(nucleotide_near).lower(),
            "all_required_remote_searches_complete": str(remote_all_complete).lower(),
            "homologous_core_phylogeny_complete": str(phylogeny_complete).lower(),
            "sequence_level_decision": decision,
        })
    write_tsv(args.out / "SEQUENCE_GATE_MATRIX.tsv", gate_rows)
    (args.out / "SEQUENCE_GATE_MATRIX.json").write_text(json.dumps(gate_rows, indent=2) + "\n")

    technical_complete = not failures
    completeness = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "technical_complete": technical_complete, "failures": failures,
        "expected_artifacts": expected_artifacts, "observed_artifacts": observed_artifacts,
        "upstream_results": upstream,
        "database_coverage_caveat": "Standard-task nr/nt coverage is split into explicit viral and cellular-organism Entrez partitions because the unfiltered remote service returned structurally invalid zero-statistic archives. Records outside those taxonomic roots or without usable taxonomy can fall outside the partitions. NCBI nt also excludes bulk WGS and some project-based TSA/environmental sequence; therefore absence of a near-identical hit is not an exhaustive GenBank novelty proof.",
        "claim_boundary": "A passing result retains a hash-locked partial Picornavirales-like sequence candidate for read-support analysis. It does not establish a formal new species, true Panax host, active replication, root-rot association, causality, pathogenicity, transmission, or agricultural/medical effect.",
    }
    (args.out / "TECHNICAL_COMPLETENESS.json").write_text(json.dumps(completeness, indent=2) + "\n")
    (args.out / "CLAIM_BOUNDARY.md").write_text(
        "# Claim boundary\n\n"
        "A passing sequence gate supports only this statement: hash-locked partial RNA-sequence candidates with Picornavirales-like PF00680/RdRP evidence were recovered from Panax notoginseng-associated root RNA-seq data. It does not establish formal virus species, the true biological host, active replication, root-rot association or causation, pathogenicity, transmission, or agricultural/medical effects. No-hit and sequence divergence are not taxonomic novelty proofs.\n\n"
        "Standard-task nr/nt coverage is partitioned into explicit viral and cellular-organism Entrez searches because the unfiltered remote service returned structurally invalid zero-statistic archives. Records outside those taxonomic roots or lacking usable taxonomy can fall outside the partitions. The NCBI nucleotide collection also does not provide one universal remote alias covering all bulk WGS, TSA, and environmental project sequence; these coverage gaps are retained explicitly rather than hidden.\n"
    )
    report = [
        "# Panax A1/A2/B final sequence gate", "",
        f"Technical completeness: **{'PASS' if technical_complete else 'FAIL/PENDING'}**", "",
        "| Candidate | Final sequence-level decision |",
        "|---|---|",
        *[f"| `{row['candidate']}` | `{row['sequence_level_decision']}` |" for row in gate_rows],
        "", "These are partial sequence-level candidates. The gate does not assign a formal taxon or biological host.",
    ]
    if failures:
        report += ["", "## Missing or failed evidence", "", *[f"- `{failure}`" for failure in failures]]
    (args.out / "FINAL_SEQUENCE_AUDIT.md").write_text("\n".join(report) + "\n")
    sums = []
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            sums.append(f"{sha(path)}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n")
    print(json.dumps(completeness, indent=2))
    return 0 if technical_complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
