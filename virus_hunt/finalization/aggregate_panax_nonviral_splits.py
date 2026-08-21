#!/usr/bin/env python3
"""Validate and aggregate split Panax nonviral BLAST requests.

The NCBI remote nr and nt complement backends reproducibly returned
structurally empty archives for combined A1+A2+B requests.  This module accepts
only the narrow mode-specific replacement contracts: exactly three requests,
each containing one immutable candidate plus every immutable positive control
required by that mode.  Every split remains an independently checksummed raw
evidence bundle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
CONTROL = "PNX_Panax_L2_control"
CONTROL_SPEC = {
    "expected_accession": "YP_009121238.1",
    "min_query_coverage": 99.0,
    "min_identity": 99.0,
}
STRATEGY = "protein_nonviral_candidate_control_splits_v1"
MODE_CONFIGS = {
    "protein_nonviral": {
        "program": "blastp", "database": "nr", "suffix": "faa",
        "controls": {CONTROL: CONTROL_SPEC},
        "root_strategy": STRATEGY,
        "child_strategy": "protein_nonviral_candidate_control_split_v1",
        "extra": (
            "-entrez_query", "all[filter] NOT txid10239[ORGN]",
            "-seg", "yes", "-comp_based_stats", "2",
        ),
    },
    "nt_nonviral": {
        "program": "blastn", "database": "nt", "suffix": "fna",
        "controls": {
            "PNX_Panax_cpDNA_control": {
                "expected_accession": "NC_026447.1",
                "min_query_coverage": 99.0,
                "min_identity": 99.0,
            },
            "PNX_NonPanax_mtDNA_control": {
                "expected_accession": "NC_012920.1",
                "min_query_coverage": 99.0,
                "min_identity": 99.0,
            },
        },
        "root_strategy": "nt_nonviral_candidate_controls_splits_v1",
        "child_strategy": "nt_nonviral_candidate_controls_split_v1",
        "extra": (
            "-task", "blastn", "-entrez_query",
            "all[filter] NOT txid10239[ORGN]", "-dust", "yes",
            "-soft_masking", "true",
        ),
    },
}
HIT_FIELDS = (
    "qseqid", "saccver", "sallacc", "sallseqid", "pident", "length",
    "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "qcovs", "staxids", "sscinames", "stitle", "qseq", "sseq",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fasta_records(path: Path) -> dict[str, dict[str, Any]]:
    raw = path.read_bytes()
    text = raw.decode("utf-8", errors="strict")
    records: dict[str, dict[str, Any]] = {}
    current: str | None = None
    body: list[str] = []
    chunks: list[str] = []

    def finish() -> None:
        nonlocal current, body, chunks
        if current is None:
            return
        sequence = "".join(line.strip() for line in body).upper()
        if not sequence:
            raise ValueError(f"empty FASTA sequence: {current}")
        payload = "".join(chunks).encode()
        if not payload.endswith(b"\n"):
            payload += b"\n"
        records[current] = {"sequence": sequence, "payload": payload}

    for line in text.splitlines(keepends=True):
        if line.startswith(">"):
            finish()
            current = line[1:].split()[0]
            if not current or current in records:
                raise ValueError("missing or duplicate FASTA record ID")
            body = []
            chunks = [line]
        elif current is None:
            if line.strip():
                raise ValueError("sequence before first FASTA header")
        else:
            body.append(line)
            chunks.append(line)
    finish()
    if not records:
        raise ValueError("FASTA has no records")
    return records


def _query_spec(record_id: str, record: dict[str, Any]) -> dict[str, Any]:
    sequence = record["sequence"]
    return {
        "id": record_id,
        "length": len(sequence),
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
    }


def _verify_sha_manifest(root: Path) -> list[str]:
    failures: list[str] = []
    manifest = root / "SHA256SUMS.txt"
    try:
        lines = manifest.read_text(errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        return [f"split_checksum_manifest_unreadable:{root.name}:{type(exc).__name__}"]
    listed: set[str] = set()
    for line in lines:
        if "  " not in line:
            failures.append(f"split_checksum_line_malformed:{root.name}")
            continue
        digest, relative = line.split("  ", 1)
        relative = relative.removeprefix("./")
        if (
            len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
            or not relative
            or Path(relative).is_absolute()
            or ".." in Path(relative).parts
            or relative == "SHA256SUMS.txt"
            or relative in listed
        ):
            failures.append(f"split_checksum_entry_invalid:{root.name}:{relative}")
            continue
        listed.add(relative)
        target = root / relative
        try:
            observed = sha256(target)
        except OSError:
            failures.append(f"split_checksummed_file_missing:{root.name}:{relative}")
            continue
        if observed != digest:
            failures.append(f"split_checksum_mismatch:{root.name}:{relative}")
    observed_files = {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    }
    if listed != observed_files:
        failures.append(
            f"split_checksum_coverage_mismatch:{root.name}:"
            f"missing={sorted(observed_files - listed)}:extra={sorted(listed - observed_files)}"
        )
    return failures


def _load_json(path: Path, failure: str, failures: list[str]) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(errors="strict"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        failures.append(failure)
        return {}
    if not isinstance(payload, dict):
        failures.append(failure)
        return {}
    return payload


def _validate_attempt_provenance(
    split_root: Path,
    candidate: str,
    status: dict[str, Any],
    failures: list[str],
) -> tuple[list[dict[str, str]], str]:
    """Bind the final archive to one canonical successful preserved attempt."""
    try:
        with (split_root / "REMOTE_ATTEMPTS.tsv").open(
            newline="", errors="strict"
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError):
        failures.append(f"split_attempt_history_unreadable:{candidate}")
        return [], ""
    if not rows:
        failures.append(f"split_attempt_history_empty:{candidate}")
        return [], ""
    if status.get("attempt_history") != rows:
        failures.append(f"split_attempt_status_mismatch:{candidate}")
    if status.get("attempt_count") != len(rows):
        failures.append(f"split_attempt_count_mismatch:{candidate}")

    try:
        attempt_count_text = (
            split_root / "ATTEMPT_COUNT.txt"
        ).read_text(errors="strict").strip()
        success_attempt_text = (
            split_root / "SUCCESS_ATTEMPT.txt"
        ).read_text(errors="strict").strip()
        attempt_count = int(attempt_count_text)
        success_attempt = int(success_attempt_text)
        search_success = (
            (split_root / "SEARCH_SUCCESS.txt").read_text(errors="strict") == "1\n"
        )
        termination_reason = (
            (split_root / "TERMINATION_REASON.txt")
            .read_text(errors="strict").strip()
        )
    except (OSError, UnicodeError, ValueError):
        failures.append(f"split_attempt_marker_invalid:{candidate}")
        attempt_count_text = ""
        success_attempt_text = ""
        attempt_count = -1
        success_attempt = -1
        search_success = False
        termination_reason = ""
    if (
        attempt_count_text != str(attempt_count)
        or success_attempt_text != str(success_attempt)
        or attempt_count != len(rows)
        or success_attempt <= 0
        or success_attempt != len(rows)
        or not search_success
        or termination_reason != "success"
        or status.get("termination_reason") != "success"
    ):
        failures.append(f"split_success_marker_mismatch:{candidate}")

    final_archive = split_root / "RESULTS.asn"
    if not final_archive.is_file() or final_archive.stat().st_size <= 0:
        failures.append(f"split_archive_missing:{candidate}")
        final_size = -1
        final_digest = ""
    else:
        final_size = final_archive.stat().st_size
        final_digest = sha256(final_archive)

    observed_attempts: list[int] = []
    success_rows: list[dict[str, str]] = []
    success_like_rows: list[dict[str, str]] = []
    for row in rows:
        try:
            attempt = int(row["attempt"])
            archive_bytes = int(row["result_archive_bytes"])
            archive_digest = row["result_archive_sha256"]
        except (KeyError, TypeError, ValueError):
            failures.append(f"split_attempt_archive_metadata_invalid:{candidate}")
            continue
        if row["attempt"] != str(attempt):
            failures.append(f"split_attempt_id_noncanonical:{candidate}")
        observed_attempts.append(attempt)
        if row["attempt"] == success_attempt_text:
            success_rows.append(row)
        if (
            all(
                row.get(key) == "0"
                for key in (
                    "blast_rc", "json_formatter_rc", "tsv_formatter_rc",
                    "validator_rc", "retryable",
                )
            )
            and row.get("failure_stage") == "none"
            and row.get("failure_class") == "none"
        ):
            success_like_rows.append(row)
        attempt_archive = (
            split_root / "ATTEMPT_ARCHIVES" / f"attempt{attempt}.asn"
        )
        if archive_bytes < 0:
            failures.append(f"split_attempt_archive_metadata_invalid:{candidate}")
        elif archive_bytes > 0:
            if (
                len(archive_digest) != 64
                or any(char not in "0123456789abcdef" for char in archive_digest)
                or not attempt_archive.is_file()
                or attempt_archive.stat().st_size != archive_bytes
                or sha256(attempt_archive) != archive_digest
            ):
                failures.append(
                    f"split_attempt_archive_mismatch:{candidate}:attempt{attempt}"
                )
        elif archive_digest != "NA" or attempt_archive.exists():
            failures.append(
                f"split_unexpected_empty_attempt_archive:{candidate}:attempt{attempt}"
            )
    if observed_attempts != list(range(1, len(rows) + 1)):
        failures.append(f"split_attempt_sequence_mismatch:{candidate}")
    if len(success_rows) != 1:
        failures.append(f"split_success_attempt_missing:{candidate}")
    if success_like_rows != success_rows:
        failures.append(f"split_success_attempt_ambiguous:{candidate}")
    if len(success_rows) == 1 and success_like_rows == success_rows:
        row = success_rows[0]
        try:
            success_size = int(row["result_archive_bytes"])
        except (KeyError, TypeError, ValueError):
            success_size = -1
        if (
            any(
                row.get(key) != "0"
                for key in (
                    "blast_rc", "json_formatter_rc", "tsv_formatter_rc",
                    "validator_rc", "retryable",
                )
            )
            or row.get("failure_stage") != "none"
            or row.get("failure_class") != "none"
            or success_size != final_size
            or row.get("result_archive_sha256") != final_digest
        ):
            failures.append(f"split_success_archive_mismatch:{candidate}")
    return rows, final_digest


def _expected_argv(
    query_argument: str, mode: str = "protein_nonviral",
) -> list[str]:
    config = MODE_CONFIGS[mode]
    split_root = str(Path(query_argument).parent)
    return [
        config["program"], "-remote", "-query", query_argument,
        "-db", config["database"],
        "-evalue", "1e-5", "-max_target_seqs", "100", "-max_hsps", "1",
        "-outfmt", "11", "-out", f"{split_root}/RESULTS.asn",
        *config["extra"],
    ]


def _report_signature(
    result_path: Path,
    candidate: str,
    lengths: dict[str, int],
    failures: list[str],
) -> str:
    payload = _load_json(
        result_path, f"split_result_json_unreadable:{candidate}", failures
    )
    reports = payload.get("BlastOutput2")
    if not isinstance(reports, list):
        failures.append(f"split_result_reports_malformed:{candidate}")
        return ""
    observed: dict[str, int] = {}
    signatures: list[dict[str, Any]] = []
    for item in reports:
        try:
            report = item["report"]
            search = report["results"]["search"]
            qid = str(search["query_title"]).split()[0]
            qlen = int(search["query_len"])
            stat = search["stat"]
            db_num = int(stat["db_num"])
            db_len = int(stat["db_len"])
            statistic = [float(stat[key]) for key in ("kappa", "lambda", "entropy")]
        except (KeyError, TypeError, ValueError, IndexError):
            failures.append(f"split_result_report_invalid:{candidate}")
            continue
        if (
            qid in observed
            or qid not in lengths
            or qlen != lengths[qid]
            or db_num <= 0
            or db_len <= 0
            or any(not math.isfinite(value) or value <= 0 for value in statistic)
        ):
            failures.append(f"split_result_report_contract_mismatch:{candidate}:{qid}")
        observed[qid] = qlen
        signatures.append({
            "program": report.get("program"),
            "version": report.get("version"),
            "reference": report.get("reference"),
            "search_target": report.get("search_target"),
            "params": report.get("params"),
            "db_num": db_num,
            "db_len": db_len,
        })
    if observed != lengths:
        failures.append(f"split_result_query_set_mismatch:{candidate}")
    canonical = {
        json.dumps(signature, sort_keys=True, separators=(",", ":"))
        for signature in signatures
    }
    if len(canonical) != 1:
        failures.append(f"split_result_signature_mismatch_within_request:{candidate}")
        return ""
    return next(iter(canonical))


def _read_hits(
    path: Path, candidate: str, control_ids: set[str], failures: list[str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    try:
        with path.open(errors="strict") as handle:
            iterator = csv.reader(handle, delimiter="\t")
            for line_number, values in enumerate(iterator, 1):
                if not values:
                    continue
                if len(values) != len(HIT_FIELDS):
                    failures.append(f"split_hit_field_count:{candidate}:{line_number}")
                    continue
                row = dict(zip(HIT_FIELDS, values))
                if row["qseqid"] not in {candidate, *control_ids}:
                    failures.append(
                        f"split_hit_unexpected_query:{candidate}:{row['qseqid']}"
                    )
                    continue
                rows.append(row)
    except (OSError, UnicodeError):
        failures.append(f"split_hits_unreadable:{candidate}")
    return rows


def _control_accession_present(
    row: dict[str, str], spec: dict[str, Any],
) -> bool:
    tokens = {row["saccver"]}
    tokens.update(token for token in row["sallacc"].split(";") if token)
    for identifier in row["sallseqid"].split(";"):
        tokens.update(token for token in identifier.split("|") if token)
    return spec["expected_accession"] in tokens


def _summarize_hits(
    rows: list[dict[str, str]],
    expected_lengths: dict[str, int],
    failures: list[str],
    label: str,
) -> dict[str, dict[str, Any]]:
    """Recompute the exact runner ``per_query`` projection from raw TSV rows."""
    summaries: dict[str, dict[str, Any]] = {}
    valid_rows: list[dict[str, str]] = []
    for row in rows:
        query_id = row["qseqid"]
        try:
            qlen = int(row["qlen"])
            bitscore = float(row["bitscore"])
            pident = float(row["pident"])
            qcovs = float(row["qcovs"])
        except (TypeError, ValueError):
            failures.append(f"split_hit_summary_numeric_invalid:{label}:{query_id}")
            continue
        if (
            query_id not in expected_lengths
            or qlen != expected_lengths[query_id]
            or any(not math.isfinite(value) for value in (bitscore, pident, qcovs))
        ):
            failures.append(f"split_hit_summary_contract_invalid:{label}:{query_id}")
            continue
        valid_rows.append(row)
    for query_id, query_length in expected_lengths.items():
        hits = sorted(
            (row for row in valid_rows if row["qseqid"] == query_id),
            key=lambda row: float(row["bitscore"]), reverse=True,
        )
        distinct: dict[str, dict[str, str]] = {}
        for hit in hits:
            distinct.setdefault(hit["saccver"], hit)
        hits = list(distinct.values())
        top = hits[0] if hits else None
        summaries[query_id] = {
            "query_length": query_length,
            "hit_count": len(hits),
            "near_identical_qcov80_pident90_count": sum(
                float(row["qcovs"]) >= 80 and float(row["pident"]) >= 90
                for row in hits
            ),
            "near_identical_qcov80_pident95_count": sum(
                float(row["qcovs"]) >= 80 and float(row["pident"]) >= 95
                for row in hits
            ),
            "top_hit": None if top is None else {
                key: top[key]
                for key in HIT_FIELDS
                if key not in {"qseq", "sseq", "sscinames"}
            },
        }
    return summaries


def _recompute_control_results(
    rows: list[dict[str, str]],
    control_specs: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    for control_id, spec in control_specs.items():
        matches = []
        for row in rows:
            if row["qseqid"] != control_id:
                continue
            try:
                valid = (
                    _control_accession_present(row, spec)
                    and float(row["pident"]) >= spec["min_identity"]
                    and float(row["qcovs"]) >= spec["min_query_coverage"]
                )
            except ValueError:
                valid = False
            if valid:
                matches.append(row)
        accession = spec["expected_accession"]
        results[control_id] = {
            **spec,
            "validated_accessions": sorted(
                {row["saccver"] for row in matches}
                | ({accession} if matches else set())
            )[:10],
            "validated": bool(matches),
        }
    return results


def _blank_query(length: int) -> dict[str, Any]:
    return {
        "query_length": length,
        "hit_count": 0,
        "near_identical_qcov80_pident90_count": 0,
        "near_identical_qcov80_pident95_count": 0,
        "top_hit": None,
    }


def aggregate(
    out: Path,
    candidate_fasta: Path,
    control_fasta: Path | list[Path],
    query_prefix: str,
    mode: str = "protein_nonviral",
) -> dict[str, Any]:
    config = MODE_CONFIGS[mode]
    control_specs: dict[str, dict[str, Any]] = config["controls"]
    control_ids = sorted(control_specs)
    control_fastas = (
        [control_fasta] if isinstance(control_fasta, Path) else control_fasta
    )
    failures: list[str] = []
    candidates = _fasta_records(candidate_fasta)
    controls: dict[str, dict[str, Any]] = {}
    for path in control_fastas:
        for control_id, record in _fasta_records(path).items():
            if control_id in controls:
                raise ValueError(f"duplicate control FASTA record: {control_id}")
            controls[control_id] = record
    if set(candidates) != set(CANDIDATES):
        raise ValueError(f"candidate FASTA set mismatch: {sorted(candidates)}")
    if set(controls) != set(control_specs):
        raise ValueError(f"control FASTA set mismatch: {sorted(controls)}")
    query_suffix = config["suffix"]
    root_query = out / f"SEARCH_QUERIES.{query_suffix}"
    expected_root_bytes = candidate_fasta.read_bytes()
    if not expected_root_bytes.endswith(b"\n"):
        expected_root_bytes += b"\n"
    for path in control_fastas:
        payload = path.read_bytes()
        expected_root_bytes += payload
        if not expected_root_bytes.endswith(b"\n"):
            expected_root_bytes += b"\n"
    if root_query.read_bytes() != expected_root_bytes:
        failures.append("aggregate_query_bytes_mismatch")

    split_summaries: dict[str, dict[str, Any]] = {}
    split_signatures: dict[str, str] = {}
    candidate_details: dict[str, dict[str, Any]] = {}
    control_details: dict[str, list[dict[str, Any]]] = {
        control_id: [] for control_id in control_ids
    }
    control_results: dict[str, list[dict[str, Any]]] = {
        control_id: [] for control_id in control_ids
    }
    aggregate_rows: list[dict[str, str]] = []
    attempt_rows: list[dict[str, str]] = []
    archive_rows: list[tuple[str, str]] = []

    for candidate in CANDIDATES:
        split_root = out / "SPLITS" / candidate
        if not split_root.is_dir():
            failures.append(f"split_missing:{candidate}")
            continue
        failures.extend(_verify_sha_manifest(split_root))
        expected_query_bytes = candidates[candidate]["payload"] + b"".join(
            controls[control_id]["payload"]
            for control_id in controls
        )
        query_argument = (
            f"{query_prefix}/SPLITS/{candidate}/SEARCH_QUERIES.{query_suffix}"
        )
        query_path = split_root / f"SEARCH_QUERIES.{query_suffix}"
        try:
            observed_query_bytes = query_path.read_bytes()
        except OSError:
            observed_query_bytes = b""
            failures.append(f"split_query_missing:{candidate}")
        if observed_query_bytes != expected_query_bytes:
            failures.append(f"split_query_bytes_mismatch:{candidate}")
        query_hash = hashlib.sha256(expected_query_bytes).hexdigest()
        candidate_spec = _query_spec(candidate, candidates[candidate])
        query_specs = [candidate_spec] + [
            _query_spec(control_id, controls[control_id])
            for control_id in controls
        ]
        expected_lengths = {row["id"]: row["length"] for row in query_specs}
        expected_payload = {
            "query_file": query_argument,
            "query_file_sha256": query_hash,
            "candidate_ids": [candidate],
            "validation_control_ids": control_ids,
            "validation_controls": [
                {"id": control_id, **control_specs[control_id]}
                for control_id in control_ids
            ],
            "queries": query_specs,
            "split_candidate_id": candidate,
        }
        observed_expected = _load_json(
            split_root / "EXPECTED_QUERIES.json",
            f"split_expected_queries_unreadable:{candidate}", failures,
        )
        if observed_expected != expected_payload:
            failures.append(f"split_expected_queries_mismatch:{candidate}")
        try:
            command = shlex.split((split_root / "COMMAND.txt").read_text(errors="strict"))
        except (OSError, UnicodeError, ValueError):
            command = []
        if command != _expected_argv(query_argument, mode):
            failures.append(f"split_command_mismatch:{candidate}")
        try:
            if (split_root / "QUERY_SHA256.txt").read_text(errors="strict") != (
                f"{query_hash}  {query_argument}\n"
            ):
                failures.append(f"split_query_sha_mismatch:{candidate}")
        except (OSError, UnicodeError):
            failures.append(f"split_query_sha_missing:{candidate}")

        status = _load_json(
            split_root / "SEARCH_STATUS.json",
            f"split_status_unreadable:{candidate}", failures,
        )
        expected_status = {
            "mode": mode,
            "database": config["database"],
            "query_file": query_argument,
            "query_sha256": query_hash,
            "query_count": 1 + len(control_ids),
            "query_ids": [candidate, *list(controls)],
            "validation_control_ids": control_ids,
            "expected_query_lengths": expected_lengths,
            "request_strategy": config["child_strategy"],
            "split_candidate_id": candidate,
        }
        for key, expected_value in expected_status.items():
            if status.get(key) != expected_value:
                failures.append(f"split_status_mismatch:{candidate}:{key}")
        per_query = status.get("per_query")
        if not isinstance(per_query, dict) or set(per_query) != {
            candidate, *control_ids
        }:
            failures.append(f"split_per_query_mismatch:{candidate}")
            per_query = {}
        control_result_map = status.get("validation_control_results")
        if not isinstance(control_result_map, dict):
            control_result_map = {}
        if (
            status.get("technical_complete") is not True
            or status.get("command_completed_successfully") is not True
            or status.get("result_archive_valid") is not True
        ):
            failures.append(f"split_not_technically_complete:{candidate}")
        for control_id in control_ids:
            control_result = control_result_map.get(control_id, {})
            if (
                not isinstance(control_result, dict)
                or control_result.get("validated") is not True
                or any(
                    control_result.get(key) != value
                    for key, value in control_specs[control_id].items()
                )
            ):
                failures.append(
                    f"split_control_status_failed:{candidate}:{control_id}"
                )

        result_signature = _report_signature(
            split_root / "RESULTS.json", candidate, expected_lengths, failures
        )
        split_signatures[candidate] = result_signature
        rows = _read_hits(
            split_root / "HITS.tsv", candidate, set(control_ids), failures
        )
        recomputed_per_query = _summarize_hits(
            rows, expected_lengths, failures, candidate
        )
        if per_query != recomputed_per_query:
            failures.append(f"split_per_query_summary_mismatch:{candidate}")
        recomputed_control_results = _recompute_control_results(
            rows, control_specs
        )
        if control_result_map != recomputed_control_results:
            failures.append(f"split_control_result_summary_mismatch:{candidate}")
        valid_control_rows = []
        validated_control_ids: set[str] = set()
        for row in rows:
            control_id = row["qseqid"]
            if control_id not in control_specs:
                continue
            spec = control_specs[control_id]
            try:
                valid = (
                    _control_accession_present(row, spec)
                    and float(row["pident"]) >= spec["min_identity"]
                    and float(row["qcovs"]) >= spec["min_query_coverage"]
                )
            except ValueError:
                valid = False
            if valid:
                valid_control_rows.append(row)
                validated_control_ids.add(control_id)
        for control_id in control_ids:
            if control_id not in validated_control_ids:
                failures.append(
                    f"split_control_hit_failed:{candidate}:{control_id}"
                )
        aggregate_rows.extend(row for row in rows if row["qseqid"] == candidate)
        aggregate_rows.extend(valid_control_rows)
        candidate_details[candidate] = recomputed_per_query.get(candidate, {})
        for control_id in control_ids:
            control_details[control_id].append(
                recomputed_per_query.get(control_id, {})
            )
            control_results[control_id].append(
                recomputed_control_results.get(control_id, {})
            )

        child_attempt_rows, archive_hash = _validate_attempt_provenance(
            split_root, candidate, status, failures
        )
        attempt_rows.extend(
            {"split_candidate_id": candidate, **row}
            for row in child_attempt_rows
        )
        if archive_hash:
            archive_rows.append((archive_hash, f"SPLITS/{candidate}/RESULTS.asn"))
        split_summaries[candidate] = {
            "relative_root": f"SPLITS/{candidate}",
            "query_sha256": query_hash,
            "result_archive_sha256": archive_hash,
            "search_status_sha256": (
                sha256(split_root / "SEARCH_STATUS.json")
                if (split_root / "SEARCH_STATUS.json").is_file() else ""
            ),
            "database_signature_sha256": (
                hashlib.sha256(result_signature.encode()).hexdigest()
                if result_signature else ""
            ),
            "technical_complete": status.get("technical_complete") is True,
            "attempt_count": status.get("attempt_count"),
        }

    nonempty_signatures = {value for value in split_signatures.values() if value}
    if len(split_signatures) != len(CANDIDATES) or len(nonempty_signatures) != 1:
        failures.append("cross_split_database_or_request_signature_mismatch")
    for control_id in control_ids:
        details = control_details[control_id]
        results = control_results[control_id]
        if len(details) != len(CANDIDATES) or any(
            item != details[0] for item in details[1:]
        ):
            failures.append(f"cross_split_control_summary_mismatch:{control_id}")
        if len(results) != len(CANDIDATES) or any(
            item != results[0] for item in results[1:]
        ):
            failures.append(f"cross_split_control_result_mismatch:{control_id}")

    # Preserve raw candidate rows and one copy of each identical control row.
    unique_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    for row in aggregate_rows:
        key = tuple(row[field] for field in HIT_FIELDS)
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)
    technical_complete = not failures
    with (out / "HITS.tsv").open("w") as handle:
        if technical_complete:
            for row in unique_rows:
                handle.write("\t".join(row[field] for field in HIT_FIELDS) + "\n")

    root_records = {**candidates, **controls}
    lengths = {name: len(record["sequence"]) for name, record in root_records.items()}
    if technical_complete:
        per_query = {candidate: candidate_details[candidate] for candidate in CANDIDATES}
        for control_id in control_ids:
            per_query[control_id] = control_details[control_id][0]
        aggregate_control_results = {
            control_id: control_results[control_id][0]
            for control_id in control_ids
        }
    else:
        per_query = {name: _blank_query(length) for name, length in lengths.items()}
        aggregate_control_results = {
            control_id: {
                **control_specs[control_id],
                "validated_accessions": [], "validated": False,
            }
            for control_id in control_ids
        }
    query_file = f"{query_prefix}/SEARCH_QUERIES.{query_suffix}"
    root_query_hash = hashlib.sha256(expected_root_bytes).hexdigest()
    expected_root = {
        "query_file": query_file,
        "query_file_sha256": root_query_hash,
        "candidate_ids": sorted(CANDIDATES),
        "validation_control_ids": control_ids,
        "validation_controls": [
            {"id": control_id, **control_specs[control_id]}
            for control_id in control_ids
        ],
        "queries": [_query_spec(candidate, candidates[candidate]) for candidate in CANDIDATES]
        + [_query_spec(control_id, controls[control_id]) for control_id in controls],
        "request_strategy": config["root_strategy"],
        "split_candidate_ids": list(CANDIDATES),
    }
    (out / "EXPECTED_QUERIES.json").write_text(json.dumps(expected_root, indent=2) + "\n")
    (out / "QUERY_SHA256.txt").write_text(f"{root_query_hash}  {query_file}\n")
    (out / "RESULT_ARCHIVES.sha256").write_text(
        "".join(f"{digest}  {relative}\n" for digest, relative in archive_rows)
    )
    split_manifest = {
        "strategy": config["root_strategy"],
        "query_prefix": query_prefix,
        "candidate_ids": list(CANDIDATES),
        "validation_control_ids": control_ids,
        "splits": split_summaries,
        "validation_failures": failures,
    }
    (out / "SPLIT_REQUESTS.json").write_text(json.dumps(split_manifest, indent=2) + "\n")

    attempt_fields = ["split_candidate_id"]
    for row in attempt_rows:
        for key in row:
            if key not in attempt_fields:
                attempt_fields.append(key)
    with (out / "REMOTE_ATTEMPTS.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=attempt_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(attempt_rows)

    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mode": mode,
        "database": config["database"],
        "request_strategy": config["root_strategy"],
        "split_candidate_ids": list(CANDIDATES),
        "split_results": split_summaries,
        "query_file": query_file,
        "query_sha256": root_query_hash,
        "query_count": len(root_records),
        "query_ids": [*CANDIDATES, *list(controls)],
        "validation_control_ids": control_ids,
        "validation_control_results": aggregate_control_results,
        "expected_query_lengths": lengths,
        "command_completed_successfully": technical_complete,
        "result_archive_valid": technical_complete,
        "attempt_count": len(attempt_rows),
        "attempt_history": attempt_rows,
        "termination_reason": "success" if technical_complete else "split_validation_failed",
        "technical_complete": technical_complete,
        "result_row_count": len(unique_rows) if technical_complete else 0,
        "unvalidated_diagnostic_row_count": 0 if technical_complete else len(aggregate_rows),
        "per_query": per_query,
        "split_validation_failures": failures,
        "annotation_validation": (
            "candidate annotations are JSON/TSV-bound within each independently "
            "validated split; every repeated control and the database signature "
            "agree across all three requests"
            if technical_complete else
            "no hit or control annotation is asserted because at least one split "
            "or cross-split binding failed"
        ),
        "interpretation_boundary": (
            "An empty hit table is evidence only when technical_complete is true; "
            "no-hit or divergence does not establish a new taxon."
        ),
    }
    (out / "SEARCH_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    (out / "SEARCH_SUCCESS.txt").write_text("1\n" if technical_complete else "0\n")
    (out / "ATTEMPT_COUNT.txt").write_text(f"{len(attempt_rows)}\n")
    (out / "SUCCESS_ATTEMPT.txt").write_text("1\n" if technical_complete else "0\n")
    (out / "TERMINATION_REASON.txt").write_text(status["termination_reason"] + "\n")
    return status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--control-fasta", type=Path, action="append", required=True)
    parser.add_argument("--query-prefix", required=True)
    parser.add_argument("--mode", choices=sorted(MODE_CONFIGS), default="protein_nonviral")
    args = parser.parse_args()
    status = aggregate(
        args.out, args.candidate_fasta, args.control_fasta, args.query_prefix,
        args.mode,
    )
    print(json.dumps(status, indent=2))
    return 0 if status["technical_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
