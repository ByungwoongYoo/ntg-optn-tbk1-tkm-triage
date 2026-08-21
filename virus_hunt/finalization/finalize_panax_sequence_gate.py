#!/usr/bin/env python3
"""Aggregate exact A1/A2/B sequence-audit artifacts and fail closed."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shlex
from datetime import datetime, timezone
from pathlib import Path


_VALIDATOR_PATH = Path(__file__).resolve().with_name(
    "validate_panax_remote_archive.py"
)
_VALIDATOR_SPEC = importlib.util.spec_from_file_location(
    "panax_remote_archive_validator", _VALIDATOR_PATH
)
if _VALIDATOR_SPEC is None or _VALIDATOR_SPEC.loader is None:
    raise RuntimeError(f"cannot load remote archive validator: {_VALIDATOR_PATH}")
_VALIDATOR_MODULE = importlib.util.module_from_spec(_VALIDATOR_SPEC)
_VALIDATOR_SPEC.loader.exec_module(_VALIDATOR_MODULE)
validate_remote_archive = _VALIDATOR_MODULE.validate


CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
CONTROLS = ("PNX_Duplo_A_control", "PNX_Duplo_B_control")
CURRENT_PANEL_FIELDS = (
    "accession", "context_group", "expected_title", "expected_length",
    "sequence_sha256", "expected_queries", "distinct_rank",
)
MODE_CONTROL_SPECS = {
    "protein_viral": {
        "PNX_Duplo_A_control": {}, "PNX_Duplo_B_control": {},
    },
    "nt_viral": {
        "PNX_Duplo_A_control": {}, "PNX_Duplo_B_control": {},
    },
    "nt_megablast": {
        "PNX_Duplo_A_control": {}, "PNX_Duplo_B_control": {},
    },
    "protein_nonviral": {
        "PNX_Panax_L2_control": {
            "expected_accession": "YP_009121238.1",
            "min_query_coverage": 99.0,
            "min_identity": 99.0,
        },
    },
    "nt_nonviral": {
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
    "nt_panax": {
        "PNX_Panax_cpDNA_control": {
            "expected_accession": "NC_026447.1",
            "min_query_coverage": 99.0,
            "min_identity": 99.0,
        },
    },
    "protein_tsa": {},
    "protein_environmental": {},
    "nt_tsa": {},
}
REMOTE_MODES = (
    "protein_viral", "protein_nonviral", "protein_tsa", "protein_environmental",
    "nt_viral", "nt_nonviral", "nt_megablast", "nt_panax", "nt_tsa",
)
REMOTE_REQUEST_CONTRACT = {
    "protein_viral": {
        "program": "blastp", "database": "nr",
        "query_argument": "queries/panax_candidates_plus_controls_orfs.faa",
        "query_parts": (("preflight", "panax_candidates_plus_controls_orfs.faa"),),
        "extra": ("-entrez_query", "txid10239[ORGN]", "-seg", "yes",
                  "-comp_based_stats", "2"),
    },
    "protein_nonviral": {
        "program": "blastp", "database": "nr",
        "query_argument": "remote-protein_nonviral/SEARCH_QUERIES.faa",
        "query_parts": (
            ("preflight", "panax_three_partial_orfs.faa"),
            ("source", "remote_partition_controls.faa"),
        ),
        "extra": ("-entrez_query", "all[filter] NOT txid10239[ORGN]",
                  "-seg", "yes", "-comp_based_stats", "2"),
    },
    "protein_tsa": {
        "program": "blastp", "database": "tsa_nr",
        "query_argument": "queries/panax_three_partial_orfs.faa",
        "query_parts": (("preflight", "panax_three_partial_orfs.faa"),),
        "extra": ("-seg", "yes", "-comp_based_stats", "2"),
    },
    "protein_environmental": {
        "program": "blastp", "database": "env_nr",
        "query_argument": "queries/panax_three_partial_orfs.faa",
        "query_parts": (("preflight", "panax_three_partial_orfs.faa"),),
        "extra": ("-seg", "yes", "-comp_based_stats", "2"),
    },
    "nt_viral": {
        "program": "blastn", "database": "nt",
        "query_argument": "queries/panax_candidates_plus_controls_contigs.fna",
        "query_parts": (("preflight", "panax_candidates_plus_controls_contigs.fna"),),
        "extra": ("-task", "blastn", "-entrez_query", "txid10239[ORGN]",
                  "-dust", "yes", "-soft_masking", "true"),
    },
    "nt_nonviral": {
        "program": "blastn", "database": "nt",
        "query_argument": "remote-nt_nonviral/SEARCH_QUERIES.fna",
        "query_parts": (
            ("preflight", "panax_three_contigs.fna"),
            ("source", "remote_partition_controls.fna"),
            ("source", "remote_nonpanax_control.fna"),
        ),
        "extra": ("-task", "blastn", "-entrez_query",
                  "all[filter] NOT txid10239[ORGN]", "-dust", "yes",
                  "-soft_masking", "true"),
    },
    "nt_megablast": {
        "program": "blastn", "database": "nt",
        "query_argument": "queries/panax_candidates_plus_controls_contigs.fna",
        "query_parts": (("preflight", "panax_candidates_plus_controls_contigs.fna"),),
        "extra": ("-task", "megablast", "-dust", "yes",
                  "-soft_masking", "true"),
    },
    "nt_panax": {
        "program": "blastn", "database": "nt",
        "query_argument": "remote-nt_panax/SEARCH_QUERIES.fna",
        "query_parts": (
            ("preflight", "panax_three_contigs.fna"),
            ("source", "remote_partition_controls.fna"),
        ),
        "extra": ("-task", "blastn", "-entrez_query", "txid44586[ORGN]",
                  "-dust", "yes", "-soft_masking", "true"),
    },
    "nt_tsa": {
        "program": "blastn", "database": "tsa_nt",
        "query_argument": "queries/panax_three_contigs.fna",
        "query_parts": (("preflight", "panax_three_contigs.fna"),),
        "extra": ("-task", "blastn", "-dust", "yes",
                  "-soft_masking", "true"),
    },
}
REMOTE_HIT_FIELDS = (
    "qseqid", "saccver", "sallacc", "sallseqid", "pident", "length",
    "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "qcovs", "staxids", "sscinames", "stitle", "qseq", "sseq",
)
PROTEIN_NONVIRAL_SPLIT_STRATEGY = (
    "protein_nonviral_candidate_control_splits_v1"
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_panel_sha(rows: list[dict[str, str]]) -> str:
    payload = [
        {field: row[field] for field in CURRENT_PANEL_FIELDS}
        for row in sorted(rows, key=lambda item: item["accession"])
    ]
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def expected_remote_argv(mode: str) -> list[str]:
    request = REMOTE_REQUEST_CONTRACT[mode]
    return [
        request["program"], "-remote", "-query", request["query_argument"],
        "-db", request["database"], "-evalue", "1e-5",
        "-max_target_seqs", "100", "-max_hsps", "1",
        "-outfmt", "11", "-out", f"remote-{mode}/RESULTS.asn",
        *request["extra"],
    ]


def _read_fasta_contract(payload: bytes) -> list[dict[str, object]]:
    text = payload.decode("utf-8", errors="strict")
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            if not name:
                raise ValueError("missing FASTA record ID")
            chunks = []
        elif name is None:
            raise ValueError("sequence before first FASTA header")
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    names = [name for name, _ in records]
    if not records or len(names) != len(set(names)) or any(not sequence for _, sequence in records):
        raise ValueError("empty or duplicated FASTA record")
    return [
        {
            "id": name,
            "length": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
        for name, sequence in records
    ]


def _fasta_record_payloads(payload: bytes) -> dict[str, bytes]:
    """Return exact newline-terminated FASTA record byte slices."""
    text = payload.decode("utf-8", errors="strict")
    records: dict[str, bytes] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in text.splitlines(keepends=True):
        if line.startswith(">"):
            if name is not None:
                record = "".join(chunks).encode()
                records[name] = record if record.endswith(b"\n") else record + b"\n"
            name = line[1:].split()[0]
            if not name or name in records:
                raise ValueError("missing or duplicate FASTA record ID")
            chunks = [line]
        elif name is None:
            if line.strip():
                raise ValueError("sequence before first FASTA header")
        else:
            chunks.append(line)
    if name is not None:
        record = "".join(chunks).encode()
        records[name] = record if record.endswith(b"\n") else record + b"\n"
    if not records:
        raise ValueError("FASTA has no records")
    return records


def _protein_nonviral_split_argv(query_argument: str) -> list[str]:
    split_root = str(Path(query_argument).parent)
    return [
        "blastp", "-remote", "-query", query_argument, "-db", "nr",
        "-evalue", "1e-5", "-max_target_seqs", "100", "-max_hsps", "1",
        "-outfmt", "11", "-out", f"{split_root}/RESULTS.asn",
        "-entrez_query", "all[filter] NOT txid10239[ORGN]",
        "-seg", "yes", "-comp_based_stats", "2",
    ]


def _validate_split_attempt_archives(
    split_root: Path, mode: str, candidate: str, child: dict,
) -> list[str]:
    failures: list[str] = []
    attempts_path = split_root / "REMOTE_ATTEMPTS.tsv"
    try:
        with attempts_path.open(newline="", errors="strict") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError):
        return [f"remote_split_attempt_history_unreadable:{mode}:{candidate}"]
    if not rows:
        return [f"remote_split_attempt_history_empty:{mode}:{candidate}"]
    label = f"{mode}:{candidate}"
    if child.get("attempt_history") != rows:
        failures.append(f"remote_split_attempt_status_mismatch:{label}")
    if child.get("attempt_count") != len(rows):
        failures.append(f"remote_split_attempt_count_mismatch:{label}")
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
            (split_root / "TERMINATION_REASON.txt").read_text(errors="strict").strip()
        )
    except (OSError, UnicodeError, ValueError):
        failures.append(f"remote_split_attempt_marker_invalid:{label}")
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
        or child.get("termination_reason") != "success"
    ):
        failures.append(f"remote_split_success_marker_mismatch:{label}")
    observed_attempts: list[int] = []
    success_like_rows: list[dict[str, str]] = []
    for row in rows:
        try:
            attempt = int(row["attempt"])
            archive_bytes = int(row["result_archive_bytes"])
            archive_digest = row["result_archive_sha256"]
        except (KeyError, TypeError, ValueError):
            failures.append(
                f"remote_split_attempt_archive_metadata_invalid:{mode}:{candidate}"
            )
            continue
        observed_attempts.append(attempt)
        if row["attempt"] != str(attempt):
            failures.append(f"remote_split_attempt_id_noncanonical:{label}")
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
            failures.append(
                f"remote_split_attempt_archive_metadata_invalid:{mode}:{candidate}"
            )
        elif archive_bytes > 0:
            if (
                len(archive_digest) != 64
                or any(char not in "0123456789abcdef" for char in archive_digest)
                or not attempt_archive.is_file()
                or attempt_archive.stat().st_size != archive_bytes
                or sha(attempt_archive) != archive_digest
            ):
                failures.append(
                    f"remote_split_attempt_archive_mismatch:{mode}:"
                    f"{candidate}:attempt{attempt}"
                )
        elif archive_digest != "NA" or attempt_archive.exists():
            failures.append(
                f"remote_split_unexpected_empty_attempt_archive:{mode}:"
                f"{candidate}:attempt{attempt}"
            )
    if observed_attempts != list(range(1, len(rows) + 1)):
        failures.append(f"remote_split_attempt_sequence_mismatch:{label}")

    final_archive = split_root / "RESULTS.asn"
    if not final_archive.is_file() or final_archive.stat().st_size <= 0:
        failures.append(f"remote_split_final_archive_invalid:{label}")
        final_size = -1
        final_digest = ""
    else:
        final_size = final_archive.stat().st_size
        final_digest = sha(final_archive)
    success_rows = [
        row for row in rows if row.get("attempt") == success_attempt_text
    ]
    if len(success_rows) != 1:
        failures.append(f"remote_split_success_attempt_missing:{label}")
    if success_like_rows != success_rows:
        failures.append(f"remote_split_success_attempt_ambiguous:{label}")
    if len(success_rows) == 1 and success_like_rows == success_rows:
        success_row = success_rows[0]
        try:
            success_size = int(success_row["result_archive_bytes"])
        except (KeyError, TypeError, ValueError):
            success_size = -1
        if (
            any(
                success_row.get(key) != "0"
                for key in (
                    "blast_rc", "json_formatter_rc", "tsv_formatter_rc",
                    "validator_rc", "retryable",
                )
            )
            or success_row.get("failure_stage") != "none"
            or success_row.get("failure_class") != "none"
            or success_size != final_size
            or success_row.get("result_archive_sha256") != final_digest
        ):
            failures.append(f"remote_split_success_archive_mismatch:{label}")
    return failures


def _summarize_split_hits(
    rows: list[dict[str, str]],
    expected_lengths: dict[str, int],
    failures: list[str],
    label: str,
) -> dict[str, dict[str, object]]:
    valid_rows: list[dict[str, str]] = []
    for row in rows:
        query_id = row["qseqid"]
        try:
            qlen = int(row["qlen"])
            bitscore = float(row["bitscore"])
            pident = float(row["pident"])
            qcovs = float(row["qcovs"])
        except (TypeError, ValueError):
            failures.append(f"remote_split_hit_summary_numeric_invalid:{label}:{query_id}")
            continue
        if (
            query_id not in expected_lengths
            or qlen != expected_lengths[query_id]
            or any(not math.isfinite(value) for value in (bitscore, pident, qcovs))
        ):
            failures.append(f"remote_split_hit_summary_contract_invalid:{label}:{query_id}")
            continue
        valid_rows.append(row)
    summaries: dict[str, dict[str, object]] = {}
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
                for key in REMOTE_HIT_FIELDS
                if key not in {"qseq", "sseq", "sscinames"}
            },
        }
    return summaries


def _split_control_results_from_hits(
    rows: list[dict[str, str]],
    control_specs: dict[str, dict[str, object]],
) -> dict[str, dict[str, object]]:
    results: dict[str, dict[str, object]] = {}
    for control_id, spec in control_specs.items():
        matches: list[dict[str, str]] = []
        accession = str(spec["expected_accession"])
        for row in rows:
            if row["qseqid"] != control_id:
                continue
            tokens = {row["saccver"]}
            tokens.update(token for token in row["sallacc"].split(";") if token)
            for identifier in row["sallseqid"].split(";"):
                tokens.update(token for token in identifier.split("|") if token)
            try:
                valid = bool(
                    accession in tokens
                    and float(row["pident"]) >= float(spec["min_identity"])
                    and float(row["qcovs"]) >= float(spec["min_query_coverage"])
                )
            except ValueError:
                valid = False
            if valid:
                matches.append(row)
        results[control_id] = {
            **spec,
            "validated_accessions": sorted(
                {row["saccver"] for row in matches}
                | ({accession} if matches else set())
            )[:10],
            "validated": bool(matches),
        }
    return results


def _validate_and_hash_split_result(
    split_root: Path,
    expected_path: Path,
    mode: str,
    candidate: str,
    failures: list[str],
) -> str:
    try:
        structural, control, _ = validate_remote_archive(
            split_root / "RESULTS.json",
            expected_path,
            mode,
            split_root / "HITS.tsv",
            query_path_override=(
                split_root / (
                    "SEARCH_QUERIES.faa"
                    if mode == "protein_nonviral"
                    else "SEARCH_QUERIES.fna"
                )
            ),
        )
    except (
        OSError, UnicodeError, ValueError, KeyError, TypeError,
        json.JSONDecodeError,
    ) as exc:
        failures.append(
            f"remote_split_archive_validator_error:{mode}:{candidate}:"
            f"{type(exc).__name__}"
        )
        return ""
    if structural:
        failures.append(
            f"remote_split_archive_structural_failure:{mode}:{candidate}"
        )
    if control:
        failures.append(f"remote_split_archive_control_failure:{mode}:{candidate}")
    try:
        payload = json.loads((split_root / "RESULTS.json").read_text(errors="strict"))
        signatures = []
        for item in payload["BlastOutput2"]:
            report = item["report"]
            stat = report["results"]["search"]["stat"]
            signatures.append({
                "program": report["program"],
                "version": report["version"],
                "reference": report["reference"],
                "search_target": report["search_target"],
                "params": report["params"],
                "db_num": int(stat["db_num"]),
                "db_len": int(stat["db_len"]),
            })
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        failures.append(f"remote_split_result_signature_invalid:{mode}:{candidate}")
        return ""
    canonical = {
        json.dumps(signature, sort_keys=True, separators=(",", ":"))
        for signature in signatures
    }
    if len(canonical) != 1:
        failures.append(f"remote_split_result_signature_mismatch:{mode}:{candidate}")
        return ""
    return hashlib.sha256(next(iter(canonical)).encode()).hexdigest()


def validate_protein_nonviral_split_contract(
    collected: Path, status: dict,
) -> list[str]:
    """Bind the three split requests to immutable sources and each other."""
    mode = "protein_nonviral"
    failures: list[str] = []
    remote_root = collected / "panax-remote-protein_nonviral"
    candidate_path = collected / "panax-query-preflight" / "panax_three_partial_orfs.faa"
    control_path = Path(__file__).resolve().parent / "remote_partition_controls.faa"
    try:
        candidate_bytes = candidate_path.read_bytes()
        control_bytes = control_path.read_bytes()
        candidate_records = _fasta_record_payloads(candidate_bytes)
        control_records = _fasta_record_payloads(control_bytes)
    except (OSError, UnicodeError, ValueError) as exc:
        return [f"remote_split_query_source_invalid:{type(exc).__name__}"]
    if set(candidate_records) != set(CANDIDATES) or set(control_records) != {
        "PNX_Panax_L2_control"
    }:
        return ["remote_split_query_source_set_mismatch"]
    control = "PNX_Panax_L2_control"
    control_specs = MODE_CONTROL_SPECS[mode]
    query_argument = "remote-protein_nonviral/SEARCH_QUERIES.faa"
    root_bytes = candidate_bytes
    if not root_bytes.endswith(b"\n"):
        root_bytes += b"\n"
    root_bytes += control_bytes
    root_sha = hashlib.sha256(root_bytes).hexdigest()
    try:
        if (remote_root / "SEARCH_QUERIES.faa").read_bytes() != root_bytes:
            failures.append("remote_constructed_query_mismatch:protein_nonviral")
    except OSError:
        failures.append("remote_constructed_query_missing:protein_nonviral")
    try:
        if (remote_root / "CANDIDATE_QUERIES.faa").read_bytes() != candidate_bytes:
            failures.append("remote_split_candidate_copy_mismatch:protein_nonviral")
    except OSError:
        failures.append("remote_split_candidate_copy_missing:protein_nonviral")
    try:
        if (remote_root / "remote_partition_controls.faa").read_bytes() != control_bytes:
            failures.append("remote_split_control_copy_mismatch:protein_nonviral")
    except OSError:
        failures.append("remote_split_control_copy_missing:protein_nonviral")

    root_queries = _read_fasta_contract(root_bytes)
    root_lengths = {row["id"]: row["length"] for row in root_queries}
    root_expected = {
        "query_file": query_argument,
        "query_file_sha256": root_sha,
        "candidate_ids": sorted(CANDIDATES),
        "validation_control_ids": [control],
        "validation_controls": [{"id": control, **control_specs[control]}],
        "queries": root_queries,
        "request_strategy": PROTEIN_NONVIRAL_SPLIT_STRATEGY,
        "split_candidate_ids": list(CANDIDATES),
    }
    if load_json(
        remote_root / "EXPECTED_QUERIES.json", failures,
        "EXPECTED_QUERIES:protein_nonviral",
    ) != root_expected:
        failures.append("remote_expected_query_contract_mismatch:protein_nonviral")
    try:
        if (remote_root / "QUERY_SHA256.txt").read_text(errors="strict") != (
            f"{root_sha}  {query_argument}\n"
        ):
            failures.append("remote_query_sha_manifest_mismatch:protein_nonviral")
    except (OSError, UnicodeError):
        failures.append("remote_query_sha_manifest_missing:protein_nonviral")

    status_contract = {
        "mode": mode,
        "database": "nr",
        "request_strategy": PROTEIN_NONVIRAL_SPLIT_STRATEGY,
        "split_candidate_ids": list(CANDIDATES),
        "query_file": query_argument,
        "query_sha256": root_sha,
        "query_count": len(root_queries),
        "query_ids": [row["id"] for row in root_queries],
        "validation_control_ids": [control],
        "expected_query_lengths": root_lengths,
        "technical_complete": True,
        "command_completed_successfully": True,
        "result_archive_valid": True,
    }
    for key, expected_value in status_contract.items():
        if status.get(key) != expected_value:
            failures.append(f"remote_status_contract_mismatch:{mode}:{key}")
    if status.get("split_validation_failures") != []:
        failures.append("remote_split_validation_failures_present:protein_nonviral")

    root_per_query = status.get("per_query")
    root_control_results = status.get("validation_control_results")
    if not isinstance(root_per_query, dict) or set(root_per_query) != {
        *CANDIDATES, control
    }:
        failures.append("remote_split_root_per_query_mismatch:protein_nonviral")
        root_per_query = {}
    if not isinstance(root_control_results, dict) or set(root_control_results) != {
        control
    }:
        failures.append("remote_split_root_control_result_mismatch:protein_nonviral")
        root_control_results = {}

    split_results = status.get("split_results")
    if not isinstance(split_results, dict) or set(split_results) != set(CANDIDATES):
        failures.append("remote_split_result_set_mismatch:protein_nonviral")
        split_results = {}
    manifest = load_json(
        remote_root / "SPLIT_REQUESTS.json", failures,
        "SPLIT_REQUESTS:protein_nonviral",
    )
    if (
        manifest.get("strategy") != PROTEIN_NONVIRAL_SPLIT_STRATEGY
        or manifest.get("query_prefix") != "remote-protein_nonviral"
        or manifest.get("candidate_ids") != list(CANDIDATES)
        or manifest.get("validation_control_ids") != [control]
        or manifest.get("splits") != split_results
        or manifest.get("validation_failures") != []
    ):
        failures.append("remote_split_manifest_contract_mismatch:protein_nonviral")

    signatures: set[str] = set()
    child_control_details: list[dict] = []
    child_control_results: list[dict] = []
    aggregate_hit_lines: list[str] = []
    seen_hit_lines: set[str] = set()
    archive_manifest_lines: list[str] = []
    for candidate in CANDIDATES:
        split_root = remote_root / "SPLITS" / candidate
        split_query_argument = (
            f"remote-protein_nonviral/SPLITS/{candidate}/SEARCH_QUERIES.faa"
        )
        split_query_bytes = candidate_records[candidate] + control_records[control]
        split_sha = hashlib.sha256(split_query_bytes).hexdigest()
        try:
            if (split_root / "SEARCH_QUERIES.faa").read_bytes() != split_query_bytes:
                failures.append(f"remote_split_query_mismatch:{candidate}")
        except OSError:
            failures.append(f"remote_split_query_missing:{candidate}")
        split_queries = _read_fasta_contract(split_query_bytes)
        split_lengths = {row["id"]: row["length"] for row in split_queries}
        child_expected = {
            "query_file": split_query_argument,
            "query_file_sha256": split_sha,
            "candidate_ids": [candidate],
            "validation_control_ids": [control],
            "validation_controls": [{"id": control, **control_specs[control]}],
            "queries": split_queries,
            "split_candidate_id": candidate,
        }
        if load_json(
            split_root / "EXPECTED_QUERIES.json", failures,
            f"EXPECTED_QUERIES:protein_nonviral:{candidate}",
        ) != child_expected:
            failures.append(f"remote_split_expected_query_mismatch:{candidate}")
        try:
            observed_argv = shlex.split(
                (split_root / "COMMAND.txt").read_text(errors="strict")
            )
        except (OSError, UnicodeError, ValueError):
            observed_argv = []
        if observed_argv != _protein_nonviral_split_argv(split_query_argument):
            failures.append(f"remote_split_command_mismatch:{candidate}")
        try:
            if (split_root / "QUERY_SHA256.txt").read_text(errors="strict") != (
                f"{split_sha}  {split_query_argument}\n"
            ):
                failures.append(f"remote_split_query_sha_mismatch:{candidate}")
        except (OSError, UnicodeError):
            failures.append(f"remote_split_query_sha_missing:{candidate}")
        child = load_json(
            split_root / "SEARCH_STATUS.json", failures,
            f"SEARCH_STATUS:protein_nonviral:{candidate}",
        )
        failures.extend(
            _validate_split_attempt_archives(split_root, mode, candidate, child)
        )
        child_contract = {
            "mode": mode,
            "database": "nr",
            "query_file": split_query_argument,
            "query_sha256": split_sha,
            "query_count": 2,
            "query_ids": [candidate, control],
            "validation_control_ids": [control],
            "expected_query_lengths": split_lengths,
            "request_strategy": "protein_nonviral_candidate_control_split_v1",
            "split_candidate_id": candidate,
            "technical_complete": True,
            "command_completed_successfully": True,
            "result_archive_valid": True,
        }
        for key, expected_value in child_contract.items():
            if child.get(key) != expected_value:
                failures.append(f"remote_split_status_mismatch:{candidate}:{key}")
        child_per_query = child.get("per_query")
        child_controls = child.get("validation_control_results")
        if not isinstance(child_per_query, dict) or set(child_per_query) != {
            candidate, control
        }:
            failures.append(f"remote_split_per_query_mismatch:{candidate}")
            child_per_query = {}
        if not isinstance(child_controls, dict) or set(child_controls) != {control}:
            failures.append(f"remote_split_control_result_mismatch:{candidate}")
            child_controls = {}
        control_result = child_controls.get(control, {})
        if (
            not isinstance(control_result, dict)
            or control_result.get("validated") is not True
            or any(
                control_result.get(key) != value
                for key, value in control_specs[control].items()
            )
        ):
            failures.append(f"remote_split_control_failed:{candidate}")
        if root_per_query.get(candidate) != child_per_query.get(candidate):
            failures.append(f"remote_split_candidate_aggregation_mismatch:{candidate}")
        archive_path = split_root / "RESULTS.asn"
        archive_hash = (
            sha(archive_path)
            if archive_path.is_file() and archive_path.stat().st_size > 0 else ""
        )
        child_status_hash = (
            sha(split_root / "SEARCH_STATUS.json")
            if (split_root / "SEARCH_STATUS.json").is_file() else ""
        )
        summary = split_results.get(candidate, {})
        if (
            not isinstance(summary, dict)
            or summary.get("relative_root") != f"SPLITS/{candidate}"
            or summary.get("query_sha256") != split_sha
            or summary.get("result_archive_sha256") != archive_hash
            or summary.get("search_status_sha256") != child_status_hash
            or summary.get("technical_complete") is not True
            or summary.get("attempt_count") != child.get("attempt_count")
        ):
            failures.append(f"remote_split_summary_mismatch:{candidate}")
        observed_signature = _validate_and_hash_split_result(
            split_root,
            split_root / "EXPECTED_QUERIES.json",
            mode,
            candidate,
            failures,
        )
        signature = summary.get("database_signature_sha256") if isinstance(summary, dict) else ""
        if signature != observed_signature or len(observed_signature) != 64:
            failures.append(f"remote_split_signature_invalid:{candidate}")
        else:
            signatures.add(observed_signature)
        archive_manifest_lines.append(
            f"{archive_hash}  SPLITS/{candidate}/RESULTS.asn\n"
        )

        try:
            raw_lines = (split_root / "HITS.tsv").read_text(errors="strict").splitlines()
        except (OSError, UnicodeError):
            raw_lines = []
            failures.append(f"remote_split_hits_unreadable:{candidate}")
        parsed_rows: list[dict[str, str]] = []
        for raw in raw_lines:
            values = raw.split("\t")
            if len(values) != len(REMOTE_HIT_FIELDS):
                failures.append(f"remote_split_hit_shape_mismatch:{candidate}")
                continue
            row = dict(zip(REMOTE_HIT_FIELDS, values))
            parsed_rows.append(row)
            include = row["qseqid"] == candidate
            if row["qseqid"] == control:
                tokens = {row["saccver"]}
                tokens.update(token for token in row["sallacc"].split(";") if token)
                for identifier in row["sallseqid"].split(";"):
                    tokens.update(token for token in identifier.split("|") if token)
                try:
                    include = bool(
                        control_specs[control]["expected_accession"] in tokens
                        and float(row["pident"]) >= control_specs[control]["min_identity"]
                        and float(row["qcovs"]) >= control_specs[control]["min_query_coverage"]
                    )
                except ValueError:
                    include = False
            elif row["qseqid"] != candidate:
                failures.append(f"remote_split_hit_query_mismatch:{candidate}")
            if include and raw not in seen_hit_lines:
                seen_hit_lines.add(raw)
                aggregate_hit_lines.append(raw)
        recomputed_per_query = _summarize_split_hits(
            parsed_rows, split_lengths, failures,
            f"{mode}:{candidate}",
        )
        if child_per_query != recomputed_per_query:
            failures.append(f"remote_split_per_query_summary_mismatch:{candidate}")
        if root_per_query.get(candidate) != recomputed_per_query.get(candidate):
            failures.append(f"remote_split_candidate_raw_aggregation_mismatch:{candidate}")
        recomputed_controls = _split_control_results_from_hits(
            parsed_rows, control_specs
        )
        if child_controls != recomputed_controls:
            failures.append(f"remote_split_control_summary_mismatch:{candidate}")
        child_control_details.append(recomputed_per_query.get(control, {}))
        child_control_results.append(recomputed_controls.get(control, {}))

    if len(signatures) != 1:
        failures.append("remote_split_cross_request_signature_mismatch:protein_nonviral")
    if (
        len(child_control_details) != len(CANDIDATES)
        or any(value != child_control_details[0] for value in child_control_details[1:])
        or root_per_query.get(control) != child_control_details[0]
    ):
        failures.append("remote_split_control_summary_aggregation_mismatch")
    if (
        len(child_control_results) != len(CANDIDATES)
        or any(value != child_control_results[0] for value in child_control_results[1:])
        or root_control_results.get(control) != child_control_results[0]
    ):
        failures.append("remote_split_control_result_aggregation_mismatch")
    try:
        expected_hits = "".join(line + "\n" for line in aggregate_hit_lines)
        if (remote_root / "HITS.tsv").read_text(errors="strict") != expected_hits:
            failures.append("remote_split_hit_aggregation_mismatch")
    except (OSError, UnicodeError):
        failures.append("remote_split_aggregate_hits_unreadable")
    try:
        if (remote_root / "RESULT_ARCHIVES.sha256").read_text(errors="strict") != "".join(
            archive_manifest_lines
        ):
            failures.append("remote_split_archive_manifest_mismatch")
    except (OSError, UnicodeError):
        failures.append("remote_split_archive_manifest_missing")
    return failures


def validate_nt_nonviral_split_contract(
    collected: Path, status: dict,
) -> list[str]:
    """Bind all three candidate+both-control nt complement requests."""
    mode = "nt_nonviral"
    failures: list[str] = []
    remote_root = collected / "panax-remote-nt_nonviral"
    preflight = collected / "panax-query-preflight" / "panax_three_contigs.fna"
    source_root = Path(__file__).resolve().parent
    control_paths = (
        source_root / "remote_partition_controls.fna",
        source_root / "remote_nonpanax_control.fna",
    )
    try:
        candidate_bytes = preflight.read_bytes()
        candidate_records = _fasta_record_payloads(candidate_bytes)
        control_payloads = [path.read_bytes() for path in control_paths]
        control_records: dict[str, bytes] = {}
        for payload in control_payloads:
            for control_id, record in _fasta_record_payloads(payload).items():
                if control_id in control_records:
                    raise ValueError("duplicate nt split control")
                control_records[control_id] = record
    except (OSError, UnicodeError, ValueError) as exc:
        return [f"remote_nt_split_query_source_invalid:{type(exc).__name__}"]
    control_specs = MODE_CONTROL_SPECS[mode]
    control_ids = sorted(control_specs)
    if set(candidate_records) != set(CANDIDATES) or set(control_records) != set(
        control_specs
    ):
        return ["remote_nt_split_query_source_set_mismatch"]
    prefix = "remote-nt_nonviral"
    root_query_argument = f"{prefix}/SEARCH_QUERIES.fna"
    root_bytes = candidate_bytes
    if not root_bytes.endswith(b"\n"):
        root_bytes += b"\n"
    for payload in control_payloads:
        root_bytes += payload
        if not root_bytes.endswith(b"\n"):
            root_bytes += b"\n"
    root_sha = hashlib.sha256(root_bytes).hexdigest()
    try:
        if (remote_root / "SEARCH_QUERIES.fna").read_bytes() != root_bytes:
            failures.append("remote_constructed_query_mismatch:nt_nonviral")
        if (remote_root / "CANDIDATE_QUERIES.fna").read_bytes() != candidate_bytes:
            failures.append("remote_split_candidate_copy_mismatch:nt_nonviral")
        for path, payload in zip(control_paths, control_payloads):
            if (remote_root / path.name).read_bytes() != payload:
                failures.append(f"remote_split_control_copy_mismatch:nt_nonviral:{path.name}")
    except OSError:
        failures.append("remote_nt_split_copied_query_material_missing")
    root_queries = _read_fasta_contract(root_bytes)
    root_lengths = {row["id"]: row["length"] for row in root_queries}
    root_expected = {
        "query_file": root_query_argument,
        "query_file_sha256": root_sha,
        "candidate_ids": sorted(CANDIDATES),
        "validation_control_ids": control_ids,
        "validation_controls": [
            {"id": control_id, **control_specs[control_id]}
            for control_id in control_ids
        ],
        "queries": root_queries,
        "request_strategy": "nt_nonviral_candidate_controls_splits_v1",
        "split_candidate_ids": list(CANDIDATES),
    }
    if load_json(
        remote_root / "EXPECTED_QUERIES.json", failures,
        "EXPECTED_QUERIES:nt_nonviral",
    ) != root_expected:
        failures.append("remote_expected_query_contract_mismatch:nt_nonviral")
    try:
        if (remote_root / "QUERY_SHA256.txt").read_text(errors="strict") != (
            f"{root_sha}  {root_query_argument}\n"
        ):
            failures.append("remote_query_sha_manifest_mismatch:nt_nonviral")
    except (OSError, UnicodeError):
        failures.append("remote_query_sha_manifest_missing:nt_nonviral")
    status_contract = {
        "mode": mode,
        "database": "nt",
        "request_strategy": "nt_nonviral_candidate_controls_splits_v1",
        "split_candidate_ids": list(CANDIDATES),
        "query_file": root_query_argument,
        "query_sha256": root_sha,
        "query_count": len(root_queries),
        "query_ids": [row["id"] for row in root_queries],
        "validation_control_ids": control_ids,
        "expected_query_lengths": root_lengths,
        "technical_complete": True,
        "command_completed_successfully": True,
        "result_archive_valid": True,
    }
    for key, expected_value in status_contract.items():
        if status.get(key) != expected_value:
            failures.append(f"remote_status_contract_mismatch:{mode}:{key}")
    if status.get("split_validation_failures") != []:
        failures.append("remote_split_validation_failures_present:nt_nonviral")
    root_per_query = status.get("per_query")
    root_controls = status.get("validation_control_results")
    if not isinstance(root_per_query, dict) or set(root_per_query) != {
        *CANDIDATES, *control_ids
    }:
        failures.append("remote_split_root_per_query_mismatch:nt_nonviral")
        root_per_query = {}
    if not isinstance(root_controls, dict) or set(root_controls) != set(control_ids):
        failures.append("remote_split_root_control_result_mismatch:nt_nonviral")
        root_controls = {}
    split_results = status.get("split_results")
    if not isinstance(split_results, dict) or set(split_results) != set(CANDIDATES):
        failures.append("remote_split_result_set_mismatch:nt_nonviral")
        split_results = {}
    manifest = load_json(
        remote_root / "SPLIT_REQUESTS.json", failures,
        "SPLIT_REQUESTS:nt_nonviral",
    )
    if (
        manifest.get("strategy") != "nt_nonviral_candidate_controls_splits_v1"
        or manifest.get("query_prefix") != prefix
        or manifest.get("candidate_ids") != list(CANDIDATES)
        or manifest.get("validation_control_ids") != control_ids
        or manifest.get("splits") != split_results
        or manifest.get("validation_failures") != []
    ):
        failures.append("remote_split_manifest_contract_mismatch:nt_nonviral")

    signatures: set[str] = set()
    child_control_details = {control_id: [] for control_id in control_ids}
    child_control_results = {control_id: [] for control_id in control_ids}
    archive_manifest_lines: list[str] = []
    aggregate_hit_lines: list[str] = []
    seen_hit_lines: set[str] = set()
    for candidate in CANDIDATES:
        split_root = remote_root / "SPLITS" / candidate
        split_query_argument = f"{prefix}/SPLITS/{candidate}/SEARCH_QUERIES.fna"
        split_query_bytes = candidate_records[candidate] + b"".join(
            control_records[control_id] for control_id in control_records
        )
        split_sha = hashlib.sha256(split_query_bytes).hexdigest()
        try:
            if (split_root / "SEARCH_QUERIES.fna").read_bytes() != split_query_bytes:
                failures.append(f"remote_split_query_mismatch:nt_nonviral:{candidate}")
        except OSError:
            failures.append(f"remote_split_query_missing:nt_nonviral:{candidate}")
        split_queries = _read_fasta_contract(split_query_bytes)
        split_lengths = {row["id"]: row["length"] for row in split_queries}
        child_expected = {
            "query_file": split_query_argument,
            "query_file_sha256": split_sha,
            "candidate_ids": [candidate],
            "validation_control_ids": control_ids,
            "validation_controls": [
                {"id": control_id, **control_specs[control_id]}
                for control_id in control_ids
            ],
            "queries": split_queries,
            "split_candidate_id": candidate,
        }
        if load_json(
            split_root / "EXPECTED_QUERIES.json", failures,
            f"EXPECTED_QUERIES:nt_nonviral:{candidate}",
        ) != child_expected:
            failures.append(f"remote_split_expected_query_mismatch:nt_nonviral:{candidate}")
        try:
            argv = shlex.split((split_root / "COMMAND.txt").read_text(errors="strict"))
        except (OSError, UnicodeError, ValueError):
            argv = []
        request = REMOTE_REQUEST_CONTRACT[mode]
        expected_argv = [
            "blastn", "-remote", "-query", split_query_argument,
            "-db", "nt", "-evalue", "1e-5", "-max_target_seqs", "100",
            "-max_hsps", "1", "-outfmt", "11", "-out",
            f"{prefix}/SPLITS/{candidate}/RESULTS.asn", *request["extra"],
        ]
        if argv != expected_argv:
            failures.append(f"remote_split_command_mismatch:nt_nonviral:{candidate}")
        try:
            if (split_root / "QUERY_SHA256.txt").read_text(errors="strict") != (
                f"{split_sha}  {split_query_argument}\n"
            ):
                failures.append(f"remote_split_query_sha_mismatch:nt_nonviral:{candidate}")
        except (OSError, UnicodeError):
            failures.append(f"remote_split_query_sha_missing:nt_nonviral:{candidate}")
        child = load_json(
            split_root / "SEARCH_STATUS.json", failures,
            f"SEARCH_STATUS:nt_nonviral:{candidate}",
        )
        failures.extend(
            _validate_split_attempt_archives(split_root, mode, candidate, child)
        )
        child_contract = {
            "mode": mode, "database": "nt",
            "query_file": split_query_argument, "query_sha256": split_sha,
            "query_count": 1 + len(control_ids),
            "query_ids": [row["id"] for row in split_queries],
            "validation_control_ids": control_ids,
            "expected_query_lengths": split_lengths,
            "request_strategy": "nt_nonviral_candidate_controls_split_v1",
            "split_candidate_id": candidate,
            "technical_complete": True,
            "command_completed_successfully": True,
            "result_archive_valid": True,
        }
        for key, expected_value in child_contract.items():
            if child.get(key) != expected_value:
                failures.append(f"remote_split_status_mismatch:nt_nonviral:{candidate}:{key}")
        child_per_query = child.get("per_query")
        controls = child.get("validation_control_results")
        if not isinstance(child_per_query, dict) or set(child_per_query) != {
            candidate, *control_ids
        }:
            failures.append(f"remote_split_per_query_mismatch:nt_nonviral:{candidate}")
            child_per_query = {}
        if not isinstance(controls, dict) or set(controls) != set(control_ids):
            failures.append(f"remote_split_control_set_mismatch:nt_nonviral:{candidate}")
            controls = {}
        for control_id in control_ids:
            result = controls.get(control_id, {})
            if (
                not isinstance(result, dict)
                or result.get("validated") is not True
                or any(
                    result.get(key) != value
                    for key, value in control_specs[control_id].items()
                )
            ):
                failures.append(
                    f"remote_split_control_failed:nt_nonviral:{candidate}:{control_id}"
                )
        if root_per_query.get(candidate) != child_per_query.get(candidate):
            failures.append(f"remote_split_candidate_aggregation_mismatch:nt_nonviral:{candidate}")
        archive_path = split_root / "RESULTS.asn"
        archive_hash = (
            sha(archive_path)
            if archive_path.is_file() and archive_path.stat().st_size > 0 else ""
        )
        child_status_hash = (
            sha(split_root / "SEARCH_STATUS.json")
            if (split_root / "SEARCH_STATUS.json").is_file() else ""
        )
        summary = split_results.get(candidate, {})
        if (
            not isinstance(summary, dict)
            or summary.get("relative_root") != f"SPLITS/{candidate}"
            or summary.get("query_sha256") != split_sha
            or summary.get("result_archive_sha256") != archive_hash
            or summary.get("search_status_sha256") != child_status_hash
            or summary.get("technical_complete") is not True
            or summary.get("attempt_count") != child.get("attempt_count")
        ):
            failures.append(f"remote_split_summary_mismatch:nt_nonviral:{candidate}")
        observed_signature = _validate_and_hash_split_result(
            split_root,
            split_root / "EXPECTED_QUERIES.json",
            mode,
            candidate,
            failures,
        )
        signature = summary.get("database_signature_sha256") if isinstance(summary, dict) else ""
        if signature != observed_signature or len(observed_signature) != 64:
            failures.append(f"remote_split_signature_invalid:nt_nonviral:{candidate}")
        else:
            signatures.add(observed_signature)
        archive_manifest_lines.append(
            f"{archive_hash}  SPLITS/{candidate}/RESULTS.asn\n"
        )
        try:
            raw_lines = (split_root / "HITS.tsv").read_text(
                errors="strict"
            ).splitlines()
        except (OSError, UnicodeError):
            raw_lines = []
            failures.append(f"remote_split_hits_unreadable:nt_nonviral:{candidate}")
        validated_controls: set[str] = set()
        parsed_rows: list[dict[str, str]] = []
        for raw in raw_lines:
            values = raw.split("\t")
            if len(values) != len(REMOTE_HIT_FIELDS):
                failures.append(
                    f"remote_split_hit_shape_mismatch:nt_nonviral:{candidate}"
                )
                continue
            row = dict(zip(REMOTE_HIT_FIELDS, values))
            parsed_rows.append(row)
            query_id = row["qseqid"]
            include = query_id == candidate
            if query_id in control_specs:
                spec = control_specs[query_id]
                tokens = {row["saccver"]}
                tokens.update(
                    token for token in row["sallacc"].split(";") if token
                )
                for identifier in row["sallseqid"].split(";"):
                    tokens.update(
                        token for token in identifier.split("|") if token
                    )
                try:
                    include = bool(
                        spec["expected_accession"] in tokens
                        and float(row["pident"]) >= spec["min_identity"]
                        and float(row["qcovs"]) >= spec["min_query_coverage"]
                    )
                except ValueError:
                    include = False
                if include:
                    validated_controls.add(query_id)
            elif query_id != candidate:
                failures.append(
                    f"remote_split_hit_query_mismatch:nt_nonviral:"
                    f"{candidate}:{query_id}"
                )
            if include and raw not in seen_hit_lines:
                seen_hit_lines.add(raw)
                aggregate_hit_lines.append(raw)
        for control_id in control_ids:
            if control_id not in validated_controls:
                failures.append(
                    f"remote_split_control_hit_failed:nt_nonviral:"
                    f"{candidate}:{control_id}"
                )
        recomputed_per_query = _summarize_split_hits(
            parsed_rows, split_lengths, failures,
            f"{mode}:{candidate}",
        )
        if child_per_query != recomputed_per_query:
            failures.append(
                f"remote_split_per_query_summary_mismatch:nt_nonviral:{candidate}"
            )
        if root_per_query.get(candidate) != recomputed_per_query.get(candidate):
            failures.append(
                f"remote_split_candidate_raw_aggregation_mismatch:nt_nonviral:"
                f"{candidate}"
            )
        recomputed_controls = _split_control_results_from_hits(
            parsed_rows, control_specs
        )
        if controls != recomputed_controls:
            failures.append(
                f"remote_split_control_summary_mismatch:nt_nonviral:{candidate}"
            )
        for control_id in control_ids:
            child_control_details[control_id].append(
                recomputed_per_query.get(control_id, {})
            )
            child_control_results[control_id].append(
                recomputed_controls.get(control_id, {})
            )
    if len(signatures) != 1:
        failures.append("remote_split_cross_request_signature_mismatch:nt_nonviral")
    for control_id in control_ids:
        details = child_control_details[control_id]
        results = child_control_results[control_id]
        if (
            len(details) != len(CANDIDATES)
            or any(value != details[0] for value in details[1:])
            or root_per_query.get(control_id) != details[0]
        ):
            failures.append(f"remote_split_control_summary_aggregation_mismatch:{control_id}")
        if (
            len(results) != len(CANDIDATES)
            or any(value != results[0] for value in results[1:])
            or root_controls.get(control_id) != results[0]
        ):
            failures.append(f"remote_split_control_result_aggregation_mismatch:{control_id}")
    try:
        if (remote_root / "RESULT_ARCHIVES.sha256").read_text(errors="strict") != "".join(
            archive_manifest_lines
        ):
            failures.append("remote_split_archive_manifest_mismatch:nt_nonviral")
    except (OSError, UnicodeError):
        failures.append("remote_split_archive_manifest_missing:nt_nonviral")
    try:
        expected_hits = "".join(line + "\n" for line in aggregate_hit_lines)
        if (remote_root / "HITS.tsv").read_text(errors="strict") != expected_hits:
            failures.append("remote_split_hit_aggregation_mismatch:nt_nonviral")
    except (OSError, UnicodeError):
        failures.append("remote_split_aggregate_hits_unreadable:nt_nonviral")
    return failures


def validate_remote_contract(
    collected: Path, mode: str, status: dict,
) -> list[str]:
    """Bind one remote result to its exact trusted query and client argv."""

    if mode == "protein_nonviral":
        return validate_protein_nonviral_split_contract(collected, status)
    if mode == "nt_nonviral":
        return validate_nt_nonviral_split_contract(collected, status)

    failures: list[str] = []
    request = REMOTE_REQUEST_CONTRACT[mode]
    remote_root = collected / f"panax-remote-{mode}"
    preflight_root = collected / "panax-query-preflight"
    source_root = Path(__file__).resolve().parent

    chunks: list[bytes] = []
    for source_kind, filename in request["query_parts"]:
        source_path = (
            preflight_root / filename if source_kind == "preflight"
            else source_root / filename
        )
        try:
            chunks.append(source_path.read_bytes())
        except OSError:
            failures.append(
                f"remote_query_source_missing:{mode}:{source_kind}:{filename}"
            )
    expected_query_bytes = b"".join(chunks)
    query_argument = request["query_argument"]
    if len(chunks) == len(request["query_parts"]):
        query_sha256 = hashlib.sha256(expected_query_bytes).hexdigest()
        try:
            expected_queries = _read_fasta_contract(expected_query_bytes)
        except (UnicodeError, ValueError) as exc:
            failures.append(f"remote_query_fasta_invalid:{mode}:{type(exc).__name__}")
            expected_queries = []
    else:
        query_sha256 = ""
        expected_queries = []

    if query_argument.startswith(f"remote-{mode}/") and expected_query_bytes:
        constructed_path = remote_root / Path(query_argument).name
        try:
            if constructed_path.read_bytes() != expected_query_bytes:
                failures.append(f"remote_constructed_query_mismatch:{mode}")
        except OSError:
            failures.append(f"remote_constructed_query_missing:{mode}")

    command_path = remote_root / "COMMAND.txt"
    try:
        observed_argv = shlex.split(command_path.read_text(errors="strict"))
    except (OSError, UnicodeError, ValueError):
        observed_argv = []
        failures.append(f"remote_command_unreadable:{mode}")
    if observed_argv != expected_remote_argv(mode):
        failures.append(f"remote_command_contract_mismatch:{mode}")

    query_sha_path = remote_root / "QUERY_SHA256.txt"
    expected_query_sha_line = f"{query_sha256}  {query_argument}\n"
    try:
        if query_sha_path.read_text(errors="strict") != expected_query_sha_line:
            failures.append(f"remote_query_sha_manifest_mismatch:{mode}")
    except (OSError, UnicodeError):
        failures.append(f"remote_query_sha_manifest_missing:{mode}")

    expected_path = remote_root / "EXPECTED_QUERIES.json"
    expected_payload = load_json(
        expected_path, failures, f"EXPECTED_QUERIES:{mode}",
    )
    control_specs = MODE_CONTROL_SPECS[mode]
    expected_control_ids = sorted(control_specs)
    expected_payload_contract = {
        "query_file": query_argument,
        "query_file_sha256": query_sha256,
        "candidate_ids": sorted(CANDIDATES),
        "validation_control_ids": expected_control_ids,
        "validation_controls": [
            {"id": control, **control_specs[control]}
            for control in expected_control_ids
        ],
        "queries": expected_queries,
    }
    if expected_payload != expected_payload_contract:
        failures.append(f"remote_expected_query_contract_mismatch:{mode}")

    expected_lengths = {
        str(row["id"]): int(row["length"]) for row in expected_queries
    }
    status_contract = {
        "mode": mode,
        "database": request["database"],
        "query_file": query_argument,
        "query_sha256": query_sha256,
        "query_count": len(expected_queries),
        "query_ids": [row["id"] for row in expected_queries],
        "validation_control_ids": expected_control_ids,
        "expected_query_lengths": expected_lengths,
    }
    for key, expected_value in status_contract.items():
        if status.get(key) != expected_value:
            failures.append(f"remote_status_contract_mismatch:{mode}:{key}")

    observed_control_results = status.get("validation_control_results", {})
    if not isinstance(observed_control_results, dict):
        failures.append(f"remote_status_control_results_malformed:{mode}")
    else:
        for control, spec in control_specs.items():
            observed = observed_control_results.get(control, {})
            if not isinstance(observed, dict):
                failures.append(
                    f"remote_status_control_spec_mismatch:{mode}:{control}"
                )
                continue
            for key in ("expected_accession", "min_query_coverage", "min_identity"):
                if observed.get(key) != spec.get(key):
                    failures.append(
                        f"remote_status_control_spec_mismatch:{mode}:{control}:{key}"
                    )
    return failures


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


def load_tsv(path: Path, failures: list[str], label: str) -> list[dict[str, str]]:
    if not path.is_file():
        failures.append(f"missing_tsv:{label}")
        return []
    try:
        with path.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
    except Exception as exc:
        failures.append(f"invalid_tsv:{label}:{type(exc).__name__}")
        return []
    if not reader.fieldnames:
        failures.append(f"headerless_tsv:{label}")
        return []
    return rows


def load_headerless_hits(path: Path, failures: list[str], label: str) -> list[dict[str, str]]:
    if not path.is_file():
        failures.append(f"missing_tsv:{label}")
        return []
    rows: list[dict[str, str]] = []
    try:
        with path.open(newline="", errors="replace") as handle:
            for number, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
                if not values:
                    continue
                if len(values) != len(REMOTE_HIT_FIELDS):
                    failures.append(f"malformed_hit_row:{label}:{number}:{len(values)}")
                    continue
                rows.append(dict(zip(REMOTE_HIT_FIELDS, values)))
    except Exception as exc:
        failures.append(f"invalid_tsv:{label}:{type(exc).__name__}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    fields = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def completed_remote_near_identical(
    candidate: str, remote_status: dict[str, dict], remote_complete: dict[str, bool],
) -> tuple[bool, bool, bool]:
    """Return protein/nt near-hit flags from complete remote modes only."""
    protein_near = False
    nucleotide_near = False
    all_complete = True
    for mode, status in remote_status.items():
        if not remote_complete.get(mode, False):
            all_complete = False
            continue
        per_query = status.get("per_query", {})
        details = per_query.get(candidate, {}) if isinstance(per_query, dict) else {}
        if not isinstance(details, dict):
            details = {}
        if mode.startswith("protein_") and details.get(
            "near_identical_qcov80_pident90_count", 0
        ):
            protein_near = True
        if mode.startswith("nt_") and details.get(
            "near_identical_qcov80_pident95_count", 0
        ):
            nucleotide_near = True
    return protein_near, nucleotide_near, all_complete


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

    local_root = args.collected / "panax-local-evidence"
    panel_rows = load_tsv(
        local_root / "current_nr_top_hit_proteins.tsv", failures,
        "CURRENT_NR_TOP_HIT_PROTEINS",
    )
    panel_contract_rows = load_tsv(
        local_root / "CURRENT_NR_REFERENCE_CONTRACT.tsv", failures,
        "CURRENT_NR_REFERENCE_CONTRACT",
    )
    panel_schema_valid = bool(panel_rows and tuple(panel_rows[0]) == CURRENT_PANEL_FIELDS)
    copied_schema_valid = bool(
        panel_contract_rows and tuple(panel_contract_rows[0]) == CURRENT_PANEL_FIELDS
    )
    if panel_rows and not panel_schema_valid:
        failures.append(f"current_nr_panel_schema_mismatch:{list(panel_rows[0])}")
    if panel_contract_rows and not copied_schema_valid:
        failures.append(
            f"current_nr_copied_contract_schema_mismatch:{list(panel_contract_rows[0])}"
        )
    if panel_rows != panel_contract_rows:
        failures.append("current_nr_panel_copied_contract_row_mismatch")
    if any(
        not row.get(field, "").strip()
        for row in panel_rows for field in CURRENT_PANEL_FIELDS
    ):
        failures.append("current_nr_panel_has_blank_contract_field")
    panel_contract_sha = canonical_panel_sha(panel_rows) if panel_schema_valid else ""
    if local.get("current_nr_panel_contract_sha256") != panel_contract_sha:
        failures.append("current_nr_local_panel_contract_hash_mismatch")
    if tree.get("current_nr_panel_contract_sha256") != panel_contract_sha:
        failures.append("current_nr_tree_panel_contract_hash_mismatch")
    panel_by_accession = {row.get("accession", ""): row for row in panel_rows}
    if len(panel_rows) != 6 or len(panel_by_accession) != 6 or "" in panel_by_accession:
        failures.append(f"current_nr_panel_not_exactly_six:{sorted(panel_by_accession)}")
    expected_rank_contract = {
        "PNX_Picorna_A1": [1, 2],
        "PNX_Picorna_A2": [1, 2],
        "PNX_Picorna_B": [1, 2, 3, 4],
    }
    observed_rank_contract: dict[str, list[int]] = {query: [] for query in CANDIDATES}
    for accession, row in panel_by_accession.items():
        try:
            rank = int(row.get("distinct_rank", ""))
        except ValueError:
            failures.append(f"invalid_current_nr_rank:{accession}")
            continue
        for query in row.get("expected_queries", "").split(";"):
            if query not in observed_rank_contract:
                failures.append(f"invalid_current_nr_query:{accession}:{query}")
                continue
            observed_rank_contract[query].append(rank)
    if {query: sorted(ranks) for query, ranks in observed_rank_contract.items()} != expected_rank_contract:
        failures.append(f"current_nr_rank_contract_mismatch:{observed_rank_contract}")

    local_provenance = load_json(
        local_root / "CURATED_REFERENCE_PROVENANCE.json", failures,
        "CURATED_REFERENCE_PROVENANCE", object_only=False,
    )
    if not isinstance(local_provenance, list):
        failures.append("current_nr_local_provenance_not_list")
        local_provenance = []
    local_current = {
        row.get("accession", ""): row for row in local_provenance
        if isinstance(row, dict) and row.get("role") == "current_nr_top_hit_context"
    }
    domain_rows = load_tsv(local_root / "DOMAIN_GATE.tsv", failures, "DOMAIN_GATE")
    domain_by_id = {row.get("sequence_id", ""): row for row in domain_rows}
    if set(local_current) != set(panel_by_accession):
        failures.append(
            f"current_nr_local_provenance_set_mismatch:{sorted(local_current)}"
        )
    for accession, panel_row in panel_by_accession.items():
        provenance = local_current.get(accession, {})
        domain = domain_by_id.get(accession, {})
        if provenance.get("context_group") != panel_row.get("context_group"):
            failures.append(f"current_nr_local_context_mismatch:{accession}")
        if provenance.get("sequence_sha256") != panel_row.get("sequence_sha256"):
            failures.append(f"current_nr_local_hash_mismatch:{accession}")
        if provenance.get("expected_sequence_sha256") != panel_row.get("sequence_sha256"):
            failures.append(f"current_nr_expected_hash_mismatch:{accession}")
        if domain.get("role") != "current_nr_top_hit_context":
            failures.append(f"current_nr_domain_role_mismatch:{accession}")
        if domain.get("context_group") != panel_row.get("context_group"):
            failures.append(f"current_nr_domain_context_mismatch:{accession}")
        if domain.get("full_sequence_sha256") != panel_row.get("sequence_sha256"):
            failures.append(f"current_nr_domain_hash_mismatch:{accession}")
        if domain.get("domain_gate_pass") != "true":
            failures.append(f"current_nr_domain_gate_failed:{accession}")

    tree_contract = tree.get("current_nr_reference_contract", {})
    if (
        tree.get("reference_count") != 25
        or tree.get("expected_tip_count") != 28
        or tree.get("current_nr_reference_count") != 6
        or tree.get("untrimmed_bound_to_raw_cores") is not True
        or tree.get("trimmed_is_ordered_column_subset") is not True
        or set(tree_contract) != set(panel_by_accession)
    ):
        failures.append("expanded_phylogeny_contract_missing_or_incomplete")
    for accession, panel_row in panel_by_accession.items():
        tree_row = tree_contract.get(accession, {})
        if tree_row.get("role") != "current_nr_top_hit_context":
            failures.append(f"current_nr_tree_role_mismatch:{accession}")
        if tree_row.get("context_group") != panel_row.get("context_group"):
            failures.append(f"current_nr_tree_context_mismatch:{accession}")
        if tree_row.get("sequence_sha256") != panel_row.get("sequence_sha256"):
            failures.append(f"current_nr_tree_hash_mismatch:{accession}")
        if tree_row.get("PF00680_core_sha256") != domain_by_id.get(accession, {}).get("core_sha256"):
            failures.append(f"current_nr_tree_core_hash_mismatch:{accession}")
        if tree_row.get("present_in_both_alignments_and_trees") is not True:
            failures.append(f"current_nr_tree_tip_missing:{accession}")

    remote_rows: list[dict[str, object]] = []
    remote_status: dict[str, dict] = {}
    remote_complete: dict[str, bool] = {}
    for mode in REMOTE_MODES:
        status = load_json(
            args.collected / f"panax-remote-{mode}" / "SEARCH_STATUS.json",
            failures, f"SEARCH_STATUS:{mode}",
        )
        remote_status[mode] = status
        contract_failures = validate_remote_contract(args.collected, mode, status)
        failures.extend(contract_failures)
        expected_controls = set(MODE_CONTROL_SPECS[mode])
        expected_ids = set(CANDIDATES) | expected_controls
        raw_query_ids = status.get("query_ids", [])
        raw_control_ids = status.get("validation_control_ids", [])
        raw_control_results = status.get("validation_control_results", {})
        raw_per_query = status.get("per_query", {})
        query_ids_valid = bool(
            isinstance(raw_query_ids, list)
            and all(isinstance(value, str) for value in raw_query_ids)
        )
        control_ids_valid = bool(
            isinstance(raw_control_ids, list)
            and all(isinstance(value, str) for value in raw_control_ids)
        )
        observed_ids = set(raw_query_ids) if query_ids_valid else set()
        observed_controls = (
            set(raw_control_ids) if control_ids_valid else set()
        )
        control_results = (
            raw_control_results if isinstance(raw_control_results, dict) else {}
        )
        per_query = raw_per_query if isinstance(raw_per_query, dict) else {}
        status_shape_valid = bool(
            query_ids_valid
            and control_ids_valid
            and isinstance(raw_control_results, dict)
            and isinstance(raw_per_query, dict)
            and all(isinstance(value, dict) for value in control_results.values())
            and all(isinstance(value, dict) for value in per_query.values())
        )
        if not status_shape_valid:
            failures.append(f"remote_status_nested_shape_mismatch:{mode}")
        controls_valid = bool(
            set(control_results) == expected_controls
            and all(
                isinstance(control_results.get(control), dict)
                and control_results[control].get("validated")
                for control in expected_controls
            )
        ) if expected_controls else control_results in ({}, None)
        observed_per_query = set(per_query)
        complete = bool(
            not contract_failures
            and status_shape_valid
            and status.get("technical_complete") and status.get("mode") == mode
            and observed_ids == expected_ids and observed_per_query == expected_ids
            and observed_controls == expected_controls
            and controls_valid
        )
        remote_complete[mode] = complete
        if not complete:
            failures.append(f"remote_search_not_complete:{mode}")
        for query in sorted(expected_ids):
            details = per_query.get(query, {})
            # A failed remote request is missing evidence, not a biological
            # zero-hit result. Keep hit-derived presentation fields blank
            # unless the complete archive/control/command contract passed.
            presented = details if complete else {}
            remote_rows.append({
                "mode": mode, "database": status.get("database", ""), "query": query,
                "technical_complete": str(complete).lower(), "hit_count": presented.get("hit_count", ""),
                "near_identical_qcov80_pident90_count": presented.get("near_identical_qcov80_pident90_count", ""),
                "near_identical_qcov80_pident95_count": presented.get("near_identical_qcov80_pident95_count", ""),
                "top_accession": (presented.get("top_hit") or {}).get("saccver", ""),
                "top_identity": (presented.get("top_hit") or {}).get("pident", ""),
                "top_qcov": (presented.get("top_hit") or {}).get("qcovs", ""),
                "top_evalue": (presented.get("top_hit") or {}).get("evalue", ""),
                "top_title": (presented.get("top_hit") or {}).get("stitle", ""),
            })
        if expected_controls and complete:
            for control in sorted(expected_controls):
                if not per_query.get(control, {}).get("hit_count"):
                    failures.append(f"pipeline_control_no_hit:{mode}:{control}")
    write_tsv(args.out / "REMOTE_COMPLETENESS.tsv", remote_rows)

    viral_hits = load_headerless_hits(
        args.collected / "panax-remote-protein_viral" / "HITS.tsv",
        failures, "PROTEIN_VIRAL_HITS",
    )
    if not remote_complete.get("protein_viral", False):
        # Preserve the raw file in its checksummed artifact, but never use rows
        # from an incomplete command/archive contract as rank evidence.
        viral_hits = []
    ranked_distinct: dict[str, list[dict[str, str]]] = {}
    for query in CANDIDATES:
        try:
            ordered = sorted(
                (row for row in viral_hits if row.get("qseqid") == query),
                key=lambda row: float(row["bitscore"]), reverse=True,
            )
        except (KeyError, ValueError):
            failures.append(f"invalid_protein_viral_bitscore:{query}")
            ordered = []
        distinct: list[dict[str, str]] = []
        seen_accessions: set[str] = set()
        for row in ordered:
            accession = row.get("saccver", "")
            if accession and accession not in seen_accessions:
                seen_accessions.add(accession)
                distinct.append(row)
        ranked_distinct[query] = distinct
    panel_validation_rows: list[dict[str, object]] = []
    for accession, panel_row in panel_by_accession.items():
        try:
            rank = int(panel_row.get("distinct_rank", ""))
        except ValueError:
            continue
        for query in panel_row.get("expected_queries", "").split(";"):
            ranked = ranked_distinct.get(query, [])
            observed = ranked[rank - 1] if 0 < rank <= len(ranked) else {}
            matched = observed.get("saccver") == accession
            if not matched:
                failures.append(
                    f"current_nr_remote_rank_mismatch:{query}:rank{rank}:"
                    f"expected={accession}:observed={observed.get('saccver','missing')}"
                )
            panel_validation_rows.append({
                "query": query,
                "distinct_rank": rank,
                "expected_accession": accession,
                "observed_accession": observed.get("saccver", ""),
                "observed_pident": observed.get("pident", ""),
                "observed_qcovs": observed.get("qcovs", ""),
                "observed_evalue": observed.get("evalue", ""),
                "observed_bitscore": observed.get("bitscore", ""),
                "rank_contract_pass": str(matched).lower(),
            })
    write_tsv(args.out / "CURRENT_NR_REFERENCE_VALIDATION.tsv", panel_validation_rows)

    local_candidates = local.get("candidate_status", {})
    gate_rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        lc = local_candidates.get(candidate, {})
        remote_protein_near, remote_nucleotide_near, remote_all_complete = (
            completed_remote_near_identical(candidate, remote_status, remote_complete)
        )
        protein_near = bool(
            lc.get("current_refseq_near_identical_protein") or remote_protein_near
        )
        nucleotide_near = bool(
            lc.get("current_refseq_near_identical_nucleotide") or remote_nucleotide_near
        )
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
        "database_coverage_caveat": "Standard-task nr/nt coverage is split into an explicit viral Entrez partition and an indexed complement defined at run time as all[filter] NOT txid10239[ORGN], because the unfiltered remote service returned structurally invalid zero-statistic archives. Entrez indexing or taxonomy gaps can still affect that operational complement. NCBI nt also excludes bulk WGS and some project-based TSA/environmental sequence; therefore absence of a near-identical hit is not an exhaustive GenBank novelty proof.",
        "claim_boundary": "A passing result retains a hash-locked partial Picornavirales-like sequence candidate for read-support analysis. It does not establish a formal new species, true Panax host, active replication, root-rot association, causality, pathogenicity, transmission, or agricultural/medical effect.",
    }
    (args.out / "TECHNICAL_COMPLETENESS.json").write_text(json.dumps(completeness, indent=2) + "\n")
    (args.out / "CLAIM_BOUNDARY.md").write_text(
        "# Claim boundary\n\n"
        "A passing sequence gate supports only this statement: hash-locked partial RNA-sequence candidates with Picornavirales-like PF00680/RdRP evidence were recovered from Panax notoginseng-associated root RNA-seq data. It does not establish formal virus species, the true biological host, active replication, root-rot association or causation, pathogenicity, transmission, or agricultural/medical effects. No-hit and sequence divergence are not taxonomic novelty proofs.\n\n"
        "Standard-task nr/nt coverage is partitioned into an explicit viral Entrez search and an indexed complement defined at run time as `all[filter] NOT txid10239[ORGN]`, because the unfiltered remote service returned structurally invalid zero-statistic archives. Entrez indexing or taxonomy gaps can still affect that operational complement. The NCBI nucleotide collection also does not provide one universal remote alias covering all bulk WGS, TSA, and environmental project sequence; these coverage gaps are retained explicitly rather than hidden.\n"
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
