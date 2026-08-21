#!/usr/bin/env python3
"""Aggregate exactly six fail-closed Panax de novo recovery audits.

The candidate-level workflow-defined rule is one passing source-run recovery
plus at least one passing non-source-run recovery on a single contig.  This is
a technical reproducibility gate for predefined partial sequences, not a
taxonomic or biological-host claim.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from audit_panax_denovo import (
    CLAIM_BOUNDARY,
    RECOVERY_FIELDS,
    load_structure,
    load_thresholds,
    read_fasta,
    sha256_path,
    write_checksums,
    write_tsv,
)


EXPECTED_RUNS = tuple(f"DRR8539{number:02d}" for number in range(7, 13))
GATE_FIELDS = (
    "query", "source_run", "source_recovery_runs",
    "non_source_recovery_runs", "source_recovery_count",
    "non_source_recovery_count", "candidate_support_rule",
    "threshold_scope", "gate_status",
)


class AggregateError(ValueError):
    """Input failure that makes the aggregate technically incomplete."""


@dataclass(frozen=True)
class AggregateOutputs:
    recovery_rows: list[dict[str, object]]
    gate_rows: list[dict[str, object]]
    technical_status: str
    failures: list[str]
    thresholds: dict[str, object]
    observed_runs: list[str]


def load_json(path: Path, label: str) -> dict[str, object]:
    if not path.is_file():
        raise AggregateError(f"missing_{label}:{path}")
    try:
        value = json.loads(path.read_text())
    except Exception as exc:  # pragma: no cover
        raise AggregateError(f"invalid_{label}_json:{type(exc).__name__}") from exc
    if not isinstance(value, dict):
        raise AggregateError(f"{label}_must_be_object")
    return value


def finite_float(value: object, label: str, minimum: float, maximum: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AggregateError(f"invalid_numeric_field:{label}:{value!r}") from exc
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise AggregateError(f"numeric_field_out_of_range:{label}:{number}")
    return number


def finite_integer(value: object, label: str, minimum: int, maximum: int) -> int:
    number = finite_float(value, label, float(minimum), float(maximum))
    if not number.is_integer():
        raise AggregateError(f"noninteger_field:{label}:{number}")
    return int(number)


def expected_thresholds(
    run: str,
    source_run: str,
    thresholds: dict[str, object],
) -> tuple[float, float, str]:
    if run == source_run:
        return (
            float(thresholds["source_query_coverage"]),
            float(thresholds["source_identity"]),
            "source",
        )
    return (
        float(thresholds["non_source_query_coverage"]),
        float(thresholds["non_source_identity"]),
        "non_source",
    )


def safe_artifact_path(
    base: Path, relative: object, run: str, label: str,
    containment_root: Path | None = None,
) -> Path:
    if not isinstance(relative, str) or not relative.strip():
        raise AggregateError(f"per_run_provenance_path:{run}:{label}")
    candidate = Path(relative)
    if candidate.is_absolute():
        raise AggregateError(f"per_run_absolute_path:{run}:{label}")
    resolved = (base / candidate).resolve()
    allowed_root = (containment_root or base).resolve()
    try:
        resolved.relative_to(allowed_root)
    except ValueError as exc:
        raise AggregateError(f"per_run_path_escape:{run}:{label}") from exc
    return resolved


def validate_retained_assembly(
    value: dict[str, object], run: str, run_root: Path
) -> None:
    de_novo_root = run_root / "de_novo"
    retained = safe_artifact_path(
        de_novo_root, value["retained_assembly_file"], run,
        "retained_assembly_file", run_root,
    )
    manifest = safe_artifact_path(
        de_novo_root, value["retained_assembly_manifest"], run,
        "retained_assembly_manifest", run_root,
    )
    if not retained.is_file() or not manifest.is_file():
        raise AggregateError(f"missing_retained_assembly_evidence:{run}")
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = [
            "uncompressed_sha256", "compressed_sha256",
            "uncompressed_bytes", "compressed_bytes",
        ]
        if reader.fieldnames != expected:
            raise AggregateError(f"retained_assembly_manifest_schema:{run}")
        rows = list(reader)
    if len(rows) != 1:
        raise AggregateError(f"retained_assembly_manifest_rows:{run}:{len(rows)}")
    row = rows[0]
    compressed_sha = sha256_path(retained)
    if row["compressed_sha256"] != compressed_sha:
        raise AggregateError(f"retained_assembly_compressed_sha256:{run}")
    try:
        declared_compressed_bytes = int(row["compressed_bytes"])
        declared_uncompressed_bytes = int(row["uncompressed_bytes"])
    except (TypeError, ValueError) as exc:
        raise AggregateError(f"retained_assembly_manifest_integer:{run}") from exc
    if declared_compressed_bytes != retained.stat().st_size:
        raise AggregateError(f"retained_assembly_compressed_size:{run}")
    digest = hashlib.sha256()
    uncompressed_bytes = 0
    try:
        with gzip.open(retained, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
                uncompressed_bytes += len(block)
    except (OSError, EOFError) as exc:
        raise AggregateError(f"retained_assembly_invalid_gzip:{run}") from exc
    if uncompressed_bytes <= 0 or uncompressed_bytes != declared_uncompressed_bytes:
        raise AggregateError(f"retained_assembly_uncompressed_size:{run}")
    uncompressed_sha = digest.hexdigest()
    if (
        row["uncompressed_sha256"] != uncompressed_sha
        or value["assembly_sha256"] != uncompressed_sha
    ):
        raise AggregateError(f"retained_assembly_uncompressed_sha256:{run}")


def validate_provenance_summary(
    value: object, run: str, run_root: Path,
    expected_query_sha: str, expected_structure_sha: str,
) -> None:
    if not isinstance(value, dict):
        raise AggregateError(f"missing_per_run_provenance:{run}")
    required = {
        "run", "assembly_sha256", "assembly_file", "input_scope",
        "candidate_baiting", "mapping_seeded", "reference_guided",
        "target_read_selection", "assembler", "assembler_version",
        "assembly_exit_code", "search_exit_code", "fastq_sha256_manifest",
        "fastq_sha256_manifest_sha256", "assembly_method_manifest",
        "assembly_method_manifest_sha256", "retained_assembly_file",
        "retained_assembly_manifest", "retained_assembly_compression",
        "candidate_query_sha256", "candidate_structure_sha256",
    }
    if not required.issubset(value):
        raise AggregateError(
            f"per_run_provenance_fields:{run}:{sorted(required-set(value))}"
        )
    if value["run"] != run or value["input_scope"] != "complete_paired_fastq":
        raise AggregateError(f"per_run_provenance_identity:{run}")
    for field in (
        "candidate_baiting", "mapping_seeded", "reference_guided",
        "target_read_selection",
    ):
        if value[field] is not False:
            raise AggregateError(f"per_run_provenance_not_seed_free:{run}:{field}")
    if value["assembler"] != "MEGAHIT" or not str(value["assembler_version"]).strip():
        raise AggregateError(f"per_run_provenance_assembler:{run}")
    if str(value["assembly_exit_code"]) != "0" or str(value["search_exit_code"]) != "0":
        raise AggregateError(f"per_run_provenance_exit_code:{run}")
    for field in (
        "assembly_sha256", "fastq_sha256_manifest_sha256",
        "assembly_method_manifest_sha256", "candidate_query_sha256",
        "candidate_structure_sha256",
    ):
        digest = value[field]
        if not isinstance(digest, str) or len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        ):
            raise AggregateError(f"per_run_provenance_sha256:{run}:{field}")
    for field in (
        "assembly_file", "fastq_sha256_manifest", "assembly_method_manifest",
        "retained_assembly_file", "retained_assembly_manifest",
    ):
        if not isinstance(value[field], str) or not value[field].strip():
            raise AggregateError(f"per_run_provenance_path:{run}:{field}")
    expected_paths = {
        "assembly_file": "megahit/final.contigs.fa",
        "fastq_sha256_manifest": "../FASTQ_SHA256.txt",
        "assembly_method_manifest": "ASSEMBLY_METHOD.tsv",
        "retained_assembly_file": "FULL_ASSEMBLY.fna.gz",
        "retained_assembly_manifest": "FULL_ASSEMBLY_MANIFEST.tsv",
    }
    for field, expected in expected_paths.items():
        if value[field] != expected:
            raise AggregateError(f"per_run_unexpected_path:{run}:{field}")
    if value["retained_assembly_compression"] != "gzip_nondeterministic_metadata_disabled":
        raise AggregateError(f"per_run_retained_compression:{run}")
    if value["candidate_query_sha256"] != expected_query_sha:
        raise AggregateError(f"per_run_query_sha256:{run}")
    if value["candidate_structure_sha256"] != expected_structure_sha:
        raise AggregateError(f"per_run_structure_sha256:{run}")
    external_provenance = load_json(
        run_root / "de_novo" / "ASSEMBLY_PROVENANCE.json",
        "assembly_provenance",
    )
    if external_provenance != value:
        raise AggregateError(f"embedded_provenance_mismatch:{run}")
    for path_field, hash_field in (
        ("fastq_sha256_manifest", "fastq_sha256_manifest_sha256"),
        ("assembly_method_manifest", "assembly_method_manifest_sha256"),
    ):
        path = safe_artifact_path(
            run_root / "de_novo", value[path_field], run, path_field,
            run_root,
        )
        if not path.is_file() or sha256_path(path) != value[hash_field]:
            raise AggregateError(f"per_run_manifest_sha256:{run}:{path_field}")
    validate_retained_assembly(value, run, run_root)


def validate_recovery_table(
    path: Path,
    audit_path: Path,
    expected_queries: dict[str, int],
    source_runs: dict[str, str],
    thresholds: dict[str, object],
    threshold_sha: str,
    query_sha: str,
    structure_sha: str,
) -> tuple[str, str, list[str], list[dict[str, object]]]:
    audit = load_json(audit_path, "de_novo_audit")
    technical_status = audit.get("technical_status")
    if technical_status not in {"pass", "technical_incomplete"}:
        raise AggregateError(
            f"invalid_per_run_technical_status:{audit_path.parent.name}:"
            f"{technical_status}"
        )
    if technical_status == "pass":
        if audit.get("technical_complete") is not True:
            raise AggregateError(f"per_run_complete_flag_mismatch:{audit_path.parent.name}")
        if audit.get("mapping_seed_free_provenance_validated") is not True:
            raise AggregateError(
                f"per_run_mapping_seed_free_not_validated:{audit_path.parent.name}"
            )
    elif audit.get("technical_complete") is not False:
        raise AggregateError(f"per_run_incomplete_flag_mismatch:{audit_path.parent.name}")
    raw_failures = audit.get("failures", [])
    if not isinstance(raw_failures, list) or any(
        not isinstance(failure, str) for failure in raw_failures
    ):
        raise AggregateError(f"invalid_per_run_failures:{audit_path.parent.name}")
    if technical_status == "pass" and raw_failures:
        raise AggregateError(f"complete_run_has_failures:{audit_path.parent.name}")
    run = audit.get("run")
    if not isinstance(run, str) or not run:
        raise AggregateError(f"invalid_per_run_audit_run:{audit_path}")
    if technical_status == "pass":
        run_root = (
            audit_path.parent.parent
            if audit_path.parent.name == "de_novo_audit"
            else audit_path.parent
        )
        validate_provenance_summary(
            audit.get("provenance"), run, run_root, query_sha, structure_sha
        )
    audit_thresholds = audit.get("workflow_defined_thresholds")
    if not isinstance(audit_thresholds, dict) or audit_thresholds.get(
        "source_file_sha256"
    ) != threshold_sha:
        raise AggregateError(f"per_run_threshold_hash_mismatch:{run}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(RECOVERY_FIELDS):
            raise AggregateError(
                f"recovery_table_fields:{run}:{reader.fieldnames}:"
                f"required={list(RECOVERY_FIELDS)}"
            )
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise AggregateError(f"malformed_recovery_row:{run}")
    if len(rows) != len(expected_queries):
        raise AggregateError(
            f"recovery_row_count:{run}:{len(rows)}!={len(expected_queries)}"
        )
    observed_queries = [row.get("query", "") for row in rows]
    if (
        len(set(observed_queries)) != len(observed_queries)
        or set(observed_queries) != set(expected_queries)
    ):
        raise AggregateError(
            f"recovery_query_set:{run}:observed={sorted(observed_queries)}:"
            f"expected={sorted(expected_queries)}"
        )

    validated: list[dict[str, object]] = []
    for row in rows:
        query = row["query"]
        if row.get("run") != run:
            raise AggregateError(f"recovery_run_mismatch:{run}:{query}:{row.get('run')}")
        if row.get("source_run") != source_runs[query]:
            raise AggregateError(f"recovery_source_run_mismatch:{run}:{query}")
        coverage_threshold, identity_threshold, role = expected_thresholds(
            run, source_runs[query], thresholds
        )
        if row.get("run_role") != role:
            raise AggregateError(f"recovery_role_mismatch:{run}:{query}:{row.get('run_role')}")
        declared_coverage_threshold = finite_float(
            row.get("query_coverage_threshold"),
            f"{run}:{query}:query_coverage_threshold", 0.0, 1.0,
        )
        declared_identity_threshold = finite_float(
            row.get("identity_threshold"),
            f"{run}:{query}:identity_threshold", 0.0, 1.0,
        )
        if not math.isclose(declared_coverage_threshold, coverage_threshold, abs_tol=1e-9):
            raise AggregateError(f"recovery_coverage_threshold_mismatch:{run}:{query}")
        if not math.isclose(declared_identity_threshold, identity_threshold, abs_tol=1e-9):
            raise AggregateError(f"recovery_identity_threshold_mismatch:{run}:{query}")
        declared_gap_threshold = finite_float(
            row.get("maximum_internal_query_gap_nt"),
            f"{run}:{query}:maximum_internal_query_gap_nt", 0.0, 1_000_000.0,
        )
        maximum_internal_gap = float(thresholds["maximum_internal_query_gap_nt"])
        if not math.isclose(declared_gap_threshold, maximum_internal_gap, abs_tol=1e-9):
            raise AggregateError(f"recovery_gap_threshold_mismatch:{run}:{query}")
        declared_subject_gap_threshold = finite_float(
            row.get("maximum_internal_subject_gap_nt"),
            f"{run}:{query}:maximum_internal_subject_gap_nt", 0.0, 1_000_000.0,
        )
        maximum_internal_subject_gap = float(
            thresholds["maximum_internal_subject_gap_nt"]
        )
        if not math.isclose(
            declared_subject_gap_threshold,
            maximum_internal_subject_gap,
            abs_tol=1e-9,
        ):
            raise AggregateError(
                f"recovery_subject_gap_threshold_mismatch:{run}:{query}"
            )
        coverage = finite_float(
            row.get("best_single_contig_query_coverage"),
            f"{run}:{query}:query_coverage", 0.0, 1.0,
        )
        query_length = finite_integer(
            row.get("query_length_nt"), f"{run}:{query}:query_length",
            1, 10_000_000,
        )
        if query_length != expected_queries[query]:
            raise AggregateError(f"recovery_query_length_mismatch:{run}:{query}")
        has_contig = bool(row.get("best_single_contig"))
        identity_text = row.get("coordinate_weighted_identity", "")
        identity = (
            finite_float(identity_text, f"{run}:{query}:identity", 0.0, 1.0)
            if identity_text != "" else 0.0
        )
        if row.get("collinear_hsp_geometry") not in {"true", "false"}:
            raise AggregateError(f"invalid_collinearity_flag:{run}:{query}")
        collinear = row["collinear_hsp_geometry"] == "true"
        gap_text = row.get("max_internal_query_gap_nt", "")
        observed_gap = (
            finite_float(
                gap_text, f"{run}:{query}:max_internal_query_gap_nt",
                0.0, 1_000_000.0,
            )
            if gap_text != "" else 0.0
        )
        subject_gap_text = row.get("max_internal_subject_gap_nt", "")
        observed_subject_gap = (
            finite_float(
                subject_gap_text,
                f"{run}:{query}:max_internal_subject_gap_nt",
                0.0, 1_000_000.0,
            )
            if subject_gap_text != "" else 0.0
        )
        covered_query_nt = finite_integer(
            row.get("covered_query_nt"), f"{run}:{query}:covered_query_nt",
            0, query_length,
        )
        hsp_count = finite_integer(
            row.get("hsp_count"), f"{run}:{query}:hsp_count", 0, 1_000_000,
        )
        orientation = row.get("relative_orientation", "")
        if has_contig:
            finite_integer(
                row.get("best_contig_length_nt"),
                f"{run}:{query}:best_contig_length", 1, 1_000_000_000,
            )
            if identity_text == "" or gap_text == "" or subject_gap_text == "":
                raise AggregateError(f"missing_best_contig_metric:{run}:{query}")
            if covered_query_nt <= 0 or hsp_count <= 0:
                raise AggregateError(f"empty_best_contig_support:{run}:{query}")
            if not math.isclose(
                coverage, covered_query_nt / query_length, abs_tol=1e-6,
            ):
                raise AggregateError(f"recovery_coverage_count_mismatch:{run}:{query}")
            if orientation not in {"plus", "minus", "mixed"}:
                raise AggregateError(f"invalid_recovery_orientation:{run}:{query}")
            if collinear and orientation not in {"plus", "minus"}:
                raise AggregateError(f"collinear_orientation_mismatch:{run}:{query}")
            finite_float(
                row.get("best_evalue"), f"{run}:{query}:best_evalue",
                0.0, 1e300,
            )
            finite_float(
                row.get("total_bitscore"), f"{run}:{query}:total_bitscore",
                0.000001, 1e300,
            )
        else:
            if (
                coverage != 0.0 or covered_query_nt != 0 or identity_text != ""
                or gap_text != "" or subject_gap_text != ""
                or row.get("best_contig_length_nt", "") != ""
                or orientation != "" or collinear or hsp_count != 0
                or row.get("best_evalue", "") != ""
                or row.get("total_bitscore", "") != ""
            ):
                raise AggregateError(f"nonempty_metric_without_contig:{run}:{query}")
        recomputed = bool(
            has_contig
            and collinear
            and coverage >= coverage_threshold
            and identity >= identity_threshold
            and observed_gap <= maximum_internal_gap
            and observed_subject_gap <= maximum_internal_subject_gap
        )
        declared_gate = row.get("gate_status")
        if declared_gate not in {"pass", "fail", "technical_incomplete"}:
            raise AggregateError(f"invalid_recovery_gate_status:{run}:{query}:{declared_gate}")
        if technical_status == "technical_incomplete":
            if declared_gate != "technical_incomplete":
                raise AggregateError(
                    f"incomplete_run_row_not_incomplete:{run}:{query}:{declared_gate}"
                )
        elif declared_gate == "technical_incomplete" or (
            (declared_gate == "pass") != recomputed
        ):
            raise AggregateError(f"recovery_gate_recompute_mismatch:{run}:{query}")
        expected_recovery_status = (
            f"{role}_technical_incomplete"
            if declared_gate == "technical_incomplete"
            else f"{role}_{'recovered' if declared_gate == 'pass' else 'not_recovered'}"
        )
        if row.get("recovery_status") != expected_recovery_status:
            raise AggregateError(f"recovery_status_mismatch:{run}:{query}")
        validated.append(row)
    recovered = sum(row["gate_status"] == "pass" for row in validated)
    expected_run_gate = (
        "technical_incomplete" if technical_status != "pass"
        else ("pass" if recovered else "fail")
    )
    if (
        audit.get("recovery_gate_status") != expected_run_gate
        or audit.get("recovered_query_count") != recovered
        or audit.get("query_count") != len(validated)
    ):
        raise AggregateError(f"per_run_audit_summary_mismatch:{run}")
    return run, str(technical_status), list(raw_failures), validated


def aggregate_audits(
    collected: Path,
    query_path: Path,
    structure_path: Path,
    threshold_path: Path,
) -> AggregateOutputs:
    if not collected.is_dir():
        raise AggregateError(f"missing_collected_directory:{collected}")
    queries = read_fasta(query_path, "query")
    thresholds = load_thresholds(threshold_path)
    structure = load_structure(structure_path, queries)
    source_runs = {query: row["source_run"] for query, row in structure.items()}
    query_lengths = {query: len(record.sequence) for query, record in queries.items()}
    threshold_sha = sha256_path(threshold_path)
    query_sha = sha256_path(query_path)
    structure_sha = sha256_path(structure_path)
    tables = sorted(collected.rglob("DE_NOVO_RECOVERY.tsv"))
    if len(tables) != len(EXPECTED_RUNS):
        raise AggregateError(
            f"per_run_recovery_table_count:{len(tables)}!={len(EXPECTED_RUNS)}"
        )

    rows_by_run: dict[str, list[dict[str, object]]] = {}
    statuses_by_run: dict[str, str] = {}
    declared_failures: list[str] = []
    for table in tables:
        run, technical_status, run_failures, rows = validate_recovery_table(
            table,
            table.parent / "DE_NOVO_AUDIT.json",
            query_lengths,
            source_runs,
            thresholds,
            threshold_sha,
            query_sha,
            structure_sha,
        )
        if run in rows_by_run:
            raise AggregateError(f"duplicate_per_run_recovery:{run}")
        rows_by_run[run] = rows
        statuses_by_run[run] = technical_status
        if technical_status == "technical_incomplete":
            if run_failures:
                declared_failures.extend(
                    f"per_run_technical_incomplete:{run}:{failure}"
                    for failure in run_failures
                )
            else:
                declared_failures.append(f"per_run_technical_incomplete:{run}")
    if set(rows_by_run) != set(EXPECTED_RUNS):
        raise AggregateError(
            f"per_run_set_mismatch:observed={sorted(rows_by_run)}:"
            f"expected={list(EXPECTED_RUNS)}"
        )

    all_rows = [
        row for run in EXPECTED_RUNS
        for row in sorted(rows_by_run[run], key=lambda item: str(item["query"]))
    ]
    gate_rows: list[dict[str, object]] = []
    for query in queries:
        passing = [
            row for row in all_rows
            if row["query"] == query and row["gate_status"] == "pass"
        ]
        source_passing = sorted(
            str(row["run"]) for row in passing if row["run"] == source_runs[query]
        )
        non_source_passing = sorted(
            str(row["run"]) for row in passing if row["run"] != source_runs[query]
        )
        passed = bool(source_passing and non_source_passing)
        gate_status = (
            "technical_incomplete"
            if any(status != "pass" for status in statuses_by_run.values())
            else ("pass" if passed else "fail")
        )
        gate_rows.append({
            "query": query,
            "source_run": source_runs[query],
            "source_recovery_runs": ",".join(source_passing),
            "non_source_recovery_runs": ",".join(non_source_passing),
            "source_recovery_count": len(source_passing),
            "non_source_recovery_count": len(non_source_passing),
            "candidate_support_rule": thresholds["candidate_support_rule"],
            "threshold_scope": "workflow_defined_technical_recovery_rule",
            "gate_status": gate_status,
        })
    return AggregateOutputs(
        recovery_rows=all_rows,
        gate_rows=gate_rows,
        technical_status=(
            "technical_incomplete"
            if any(status != "pass" for status in statuses_by_run.values())
            else "pass"
        ),
        failures=declared_failures,
        thresholds=thresholds,
        observed_runs=list(EXPECTED_RUNS),
    )


def write_outputs(out: Path, result: AggregateOutputs) -> dict[str, object]:
    write_tsv(out / "ALL_DE_NOVO_RECOVERY.tsv", RECOVERY_FIELDS, result.recovery_rows)
    write_tsv(out / "DE_NOVO_CANDIDATE_GATE.tsv", GATE_FIELDS, result.gate_rows)
    if result.technical_status != "pass":
        overall = "technical_incomplete"
    elif result.gate_rows and all(row["gate_status"] == "pass" for row in result.gate_rows):
        overall = "pass"
    else:
        overall = "fail"
    audit = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "technical_status": result.technical_status,
        "technical_complete": result.technical_status == "pass",
        "failures": result.failures,
        "overall_gate_status": overall,
        "expected_runs": list(EXPECTED_RUNS),
        "per_run_recovery_table_count": len(result.observed_runs),
        "observed_runs": result.observed_runs,
        "workflow_defined_thresholds": result.thresholds,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (out / "DE_NOVO_AGGREGATE_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    report = [
        "# Six-run mapping-seed-free de novo recovery gate", "",
        f"Technical status: **{audit['technical_status']}**", "",
        f"Overall gate status: **{audit['overall_gate_status']}**", "",
        "| Query | Source run | Passing source | Passing non-source | Gate |",
        "|---|---|---|---|---|",
    ]
    for row in result.gate_rows:
        report.append(
            f"| `{row['query']}` | `{row['source_run']}` | "
            f"`{row['source_recovery_runs']}` | "
            f"`{row['non_source_recovery_runs']}` | `{row['gate_status']}` |"
        )
    report.extend(["", CLAIM_BOUNDARY])
    if result.failures:
        report.extend(["", "## Technical failures", ""])
        report.extend(f"- `{failure}`" for failure in result.failures)
    (out / "DE_NOVO_AGGREGATE_REPORT.md").write_text("\n".join(report) + "\n")
    write_checksums(out)
    return audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collected", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")
    try:
        result = aggregate_audits(
            collected=args.collected,
            query_path=args.queries,
            structure_path=args.structure,
            threshold_path=args.thresholds,
        )
    except (AggregateError, ValueError) as exc:
        result = AggregateOutputs(
            recovery_rows=[],
            gate_rows=[],
            technical_status="technical_incomplete",
            failures=[str(exc)],
            thresholds={},
            observed_runs=[],
        )
    audit = write_outputs(args.out, result)
    print(json.dumps(audit, indent=2))
    return 0 if audit["overall_gate_status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
