#!/usr/bin/env python3
"""Fail-closed validation and classification for Panax remote BLAST archives.

Exit status 20 means the archive/result structure is incomplete or invalid and
is eligible for a bounded remote retry.  Exit status 21 means the result is
structurally complete but a same-request positive control failed, which is a
deterministic partition/control failure.  Other nonzero statuses are reserved
for validator invocation or implementation errors.
"""

from __future__ import annotations

import json
import hashlib
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any


EXIT_STRUCTURAL = 20
EXIT_CONTROL = 21
CANDIDATES = {"PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B"}
MODE_CONTROL_CONTRACT: dict[str, dict[str, dict[str, Any]]] = {
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
MODE_SEARCH_CONTRACT: dict[str, dict[str, Any]] = {
    "protein_viral": {
        "program": "blastp", "database": "nr", "expect": 1e-5,
        "params": {"expect": 1e-5, "matrix": "BLOSUM62", "gap_open": 11,
                   "gap_extend": 1, "filter": "L;", "cbs": 2},
    },
    "protein_nonviral": {
        "program": "blastp", "database": "nr", "expect": 1e-5,
        "params": {"expect": 1e-5, "matrix": "BLOSUM62", "gap_open": 11,
                   "gap_extend": 1, "filter": "L;", "cbs": 2},
    },
    "protein_tsa": {
        "program": "blastp", "database": "tsa_nr", "expect": 1e-5,
        "params": {"expect": 1e-5, "matrix": "BLOSUM62", "gap_open": 11,
                   "gap_extend": 1, "filter": "L;", "cbs": 2},
    },
    "protein_environmental": {
        "program": "blastp", "database": "env_nr", "expect": 1e-5,
        "params": {"expect": 1e-5, "matrix": "BLOSUM62", "gap_open": 11,
                   "gap_extend": 1, "filter": "L;", "cbs": 2},
    },
    "nt_viral": {
        "program": "blastn", "database": "nt", "expect": 1e-5,
        "params": {"expect": 1e-5, "sc_match": 2, "sc_mismatch": -3,
                   "gap_open": 5, "gap_extend": 2, "filter": "L;m;"},
    },
    "nt_nonviral": {
        "program": "blastn", "database": "nt", "expect": 1e-5,
        "params": {"expect": 1e-5, "sc_match": 2, "sc_mismatch": -3,
                   "gap_open": 5, "gap_extend": 2, "filter": "L;m;"},
    },
    "nt_megablast": {
        "program": "blastn", "database": "nt", "expect": 1e-5,
        "params": {"expect": 1e-5, "sc_match": 1, "sc_mismatch": -2,
                   "gap_open": 0, "gap_extend": 0, "filter": "L;m;"},
    },
    "nt_panax": {
        "program": "blastn", "database": "nt", "expect": 1e-5,
        "params": {"expect": 1e-5, "sc_match": 2, "sc_mismatch": -3,
                   "gap_open": 5, "gap_extend": 2, "filter": "L;m;"},
    },
    "nt_tsa": {
        "program": "blastn", "database": "tsa_nt", "expect": 1e-5,
        "params": {"expect": 1e-5, "sc_match": 2, "sc_mismatch": -3,
                   "gap_open": 5, "gap_extend": 2, "filter": "L;m;"},
    },
}
HIT_FIELDS = [
    "qseqid",
    "saccver",
    "sallacc",
    "sallseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
    "staxids",
    "sscinames",
    "stitle",
    "qseq",
    "sseq",
]


def _accession_token(value: str) -> str:
    token = value.strip().strip("|")
    if "|" in token:
        parts = [part for part in token.split("|") if part]
        token = parts[-1] if parts else token
    return token


def _semicolon_tokens(value: Any) -> list[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def _tsv_sallacc_keys(row: dict[str, Any]) -> set[str]:
    return set(_semicolon_tokens(row["sallacc"]))


def _tsv_sallseqid_keys(row: dict[str, Any]) -> set[str]:
    return set(_semicolon_tokens(row["sallseqid"]))


def _tsv_accession_keys(row: dict[str, Any]) -> set[str]:
    values = {
        str(row["saccver"]),
        *_tsv_sallacc_keys(row),
        *_tsv_sallseqid_keys(row),
    }
    return {_accession_token(value) for value in values}


def _json_accession_keys(hit: dict[str, Any]) -> set[str]:
    descriptions = hit.get("description", [])
    if not isinstance(descriptions, list):
        return set()
    values: set[str] = set()
    for description in descriptions:
        if not isinstance(description, dict):
            continue
        for key in ("accession", "id"):
            value = description.get(key, "")
            if isinstance(value, str) and value.strip():
                values.add(_accession_token(value))
    return {value for value in values if value}


def _json_description_key_sets(hit: dict[str, Any]) -> tuple[set[str], set[str]]:
    descriptions = hit.get("description", [])
    if not isinstance(descriptions, list):
        return set(), set()
    accession_keys: set[str] = set()
    id_keys: set[str] = set()
    for description in descriptions:
        if not isinstance(description, dict):
            continue
        accession = description.get("accession")
        identifier = description.get("id")
        if isinstance(accession, str) and accession.strip():
            accession_keys.add(accession.strip())
        if isinstance(identifier, str) and identifier.strip():
            id_keys.add(identifier.strip())
    return accession_keys, id_keys


def _primary_saccver(raw_accession: str, raw_identifier: str) -> str:
    """Derive BLAST's primary ``saccver`` without discarding pipe context."""
    accession = raw_accession.strip()
    accession_base = re.sub(r"\.\d+$", "", accession)
    identifier_parts = [part for part in raw_identifier.split("|") if part]
    # BLAST renders PDB chain identifiers as ``pdb|ENTRY|CHAIN`` while
    # ``saccver``/description.accession uses ``ENTRY_CHAIN``.
    if len(identifier_parts) >= 3 and identifier_parts[0].lower() == "pdb":
        pdb_chain = f"{identifier_parts[1]}_{identifier_parts[2]}"
        if pdb_chain == accession:
            return accession
    for token in identifier_parts:
        if re.sub(r"\.\d+$", "", token) == accession_base:
            return token
    return ""


def _json_description_records(hit: dict[str, Any]) -> list[dict[str, Any]]:
    descriptions = hit.get("description", [])
    if not isinstance(descriptions, list):
        return []
    records: list[dict[str, Any]] = []
    for description in descriptions:
        if not isinstance(description, dict):
            continue
        raw_accession = description.get("accession")
        raw_identifier = description.get("id")
        if (
            not isinstance(raw_accession, str)
            or not raw_accession.strip()
            or not isinstance(raw_identifier, str)
            or not raw_identifier.strip()
        ):
            continue
        accession_keys: set[str] = set()
        for key in ("accession", "id"):
            value = description.get(key)
            if isinstance(value, str) and value.strip():
                accession_keys.add(_accession_token(value))
        title = description.get("title")
        if not accession_keys or not isinstance(title, str) or not title:
            continue
        raw_taxid = description.get("taxid")
        taxids: set[str] = set()
        if isinstance(raw_taxid, int) and not isinstance(raw_taxid, bool):
            if raw_taxid > 0:
                taxids.add(str(raw_taxid))
        if not taxids:
            continue
        primary_saccver = _primary_saccver(raw_accession, raw_identifier)
        if not primary_saccver:
            continue
        records.append({
            "accession_keys": accession_keys,
            "accession": _accession_token(raw_accession),
            "id": primary_saccver,
            "title": title,
            "taxids": taxids,
        })
    return records


def _expected_accession_present(expected: str, observed_keys: set[str]) -> bool:
    if re.search(r"\.\d+$", expected):
        return expected in observed_keys
    expected_base = re.sub(r"\.\d+$", "", expected)
    return any(re.sub(r"\.\d+$", "", key) == expected_base for key in observed_keys)


def _strict_integer(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("expected an integer JSON value")
    return value


def _strict_number(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("expected a numeric JSON value")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("expected a finite JSON value")
    return result


def validate(
    result_path: Path,
    expected_path: Path,
    mode: str,
    hits_path: Path,
) -> tuple[list[str], list[str], int]:
    """Return structural errors, control errors, and observed report count."""

    try:
        payload = json.loads(result_path.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return [f"unreadable or malformed result JSON: {exc}"], [], 0
    try:
        expected_payload = json.loads(expected_path.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"unreadable or malformed local expectation JSON: {exc}") from exc
    if not isinstance(payload, dict):
        return ["result JSON root is not an object"], [], 0
    if not isinstance(expected_payload, dict):
        raise ValueError("local expectation JSON root is not an object")
    queries = expected_payload.get("queries")
    if not isinstance(queries, list) or not queries:
        raise ValueError("local expectation JSON has no nonempty queries list")
    expected_ids = [item.get("id") for item in queries if isinstance(item, dict)]
    if len(expected_ids) != len(queries) or len(expected_ids) != len(set(expected_ids)):
        raise ValueError("local expectation JSON has malformed or duplicate query IDs")
    expected_lengths = {item["id"]: int(item["length"]) for item in queries}
    expected_hashes = {item["id"]: str(item.get("sequence_sha256", "")) for item in queries}
    query_path = Path(str(expected_payload.get("query_file", "")))
    if not query_path.is_file():
        raise ValueError(f"expected local query FASTA is absent: {query_path}")
    observed_query_file_hash = hashlib.sha256(query_path.read_bytes()).hexdigest()
    if observed_query_file_hash != expected_payload.get("query_file_sha256"):
        raise ValueError("local query FASTA file hash does not match expectation JSON")
    query_sequences: dict[str, str] = {}
    query_name: str | None = None
    for raw in query_path.read_text(errors="strict").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            query_name = line[1:].split()[0]
            if not query_name or query_name in query_sequences:
                raise ValueError("local query FASTA has a missing or duplicate ID")
            query_sequences[query_name] = ""
        elif query_name is None:
            raise ValueError("local query FASTA has sequence before its first header")
        else:
            query_sequences[query_name] += line.upper()
    if set(query_sequences) != set(expected_lengths):
        raise ValueError("local query FASTA ID set does not match expectation JSON")
    for query_id, sequence in query_sequences.items():
        if (
            len(sequence) != expected_lengths[query_id]
            or hashlib.sha256(sequence.encode()).hexdigest() != expected_hashes[query_id]
        ):
            raise ValueError(f"local query FASTA sequence contract failed: {query_id}")
    control_specs = expected_payload.get("validation_controls", [])
    control_ids = expected_payload.get("validation_control_ids", [])
    if mode not in MODE_CONTROL_CONTRACT:
        raise ValueError(f"unsupported remote-search mode in validator: {mode}")
    required_control_contract = MODE_CONTROL_CONTRACT[mode]
    search_contract = MODE_SEARCH_CONTRACT[mode]
    if set(expected_lengths) != CANDIDATES | set(required_control_contract):
        raise ValueError(
            "local expectation query set does not match the mode-specific "
            "A1/A2/B plus control contract"
        )
    if (
        not isinstance(control_specs, list)
        or not isinstance(control_ids, list)
        or {spec.get("id") for spec in control_specs if isinstance(spec, dict)} != set(control_ids)
        or not set(control_ids).issubset(expected_lengths)
        or set(control_ids) != set(required_control_contract)
    ):
        raise ValueError("local validation-control ID/spec contract is inconsistent")
    specs_by_id = {spec["id"]: spec for spec in control_specs}
    for control_id, required in required_control_contract.items():
        observed_spec = specs_by_id[control_id]
        for key, expected_value in required.items():
            observed_value = observed_spec.get(key)
            if key.startswith("min_"):
                try:
                    observed_value = float(observed_value)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"invalid local control threshold: {control_id}/{key}"
                    ) from exc
            if observed_value != expected_value:
                raise ValueError(
                    f"local validation-control contract mismatch: {control_id}/{key}"
                )
    reports = payload.get("BlastOutput2", [])
    structural_errors: list[str] = []
    control_errors: list[str] = []
    observed: dict[str, int] = {}
    json_hit_counts: dict[str, int] = {}
    json_hsp_counts: dict[str, int] = {}
    json_hsp_records: list[dict[str, Any]] = []
    request_signatures: list[str] = []
    database_signatures: list[tuple[int, int]] = []

    if not isinstance(reports, list) or not reports:
        structural_errors.append("missing BlastOutput2 reports")
    else:
        for index, item in enumerate(reports, 1):
            try:
                report = item["report"]
                search = report["results"]["search"]
                if not isinstance(report, dict):
                    raise TypeError("report payload is not an object")
                if not isinstance(search, dict):
                    raise TypeError("search payload is not an object")
            except (KeyError, TypeError) as exc:  # archive schema is external input
                structural_errors.append(
                    f"report {index} has no search payload: {exc}"
                )
                continue
            program = report.get("program")
            search_target = report.get("search_target")
            params = report.get("params")
            version = report.get("version")
            reference = report.get("reference")
            expected_version_program = search_contract["program"].upper()
            if (
                not isinstance(version, str)
                or re.fullmatch(
                    rf"{re.escape(expected_version_program)} \d+\.\d+\.\d+\+?",
                    version.strip(),
                ) is None
            ):
                structural_errors.append(
                    f"report {index} BLAST version/program mismatch: {version!r}"
                )
            if not isinstance(reference, str) or not reference.strip():
                structural_errors.append(f"report {index} has no BLAST reference")
            if not isinstance(params, dict) or set(params) != set(search_contract["params"]):
                observed_param_keys = sorted(params) if isinstance(params, dict) else []
                structural_errors.append(
                    f"report {index} parameter-key mismatch: "
                    f"observed={observed_param_keys}, "
                    f"expected={sorted(search_contract['params'])}"
                )
            request_signatures.append(json.dumps({
                "program": program,
                "version": version,
                "reference": reference,
                "search_target": search_target,
                "params": params,
            }, sort_keys=True, separators=(",", ":")))
            if program != search_contract["program"]:
                structural_errors.append(
                    f"report {index} program mismatch: observed={program!r}, "
                    f"expected={search_contract['program']!r}"
                )
            if (
                not isinstance(search_target, dict)
                or set(search_target) != {"db"}
                or search_target.get("db") != search_contract["database"]
            ):
                observed_db = (
                    search_target.get("db") if isinstance(search_target, dict) else None
                )
                structural_errors.append(
                    f"report {index} database mismatch: observed={observed_db!r}, "
                    f"expected={search_contract['database']!r}"
                )
            try:
                observed_expect = _strict_number(
                    params.get("expect") if isinstance(params, dict) else None
                )
            except ValueError:
                observed_expect = math.nan
            if not math.isclose(
                observed_expect, float(search_contract["expect"]),
                rel_tol=1e-12, abs_tol=0.0,
            ):
                structural_errors.append(
                    f"report {index} expect mismatch: observed={observed_expect!r}, "
                    f"expected={search_contract['expect']!r}"
                )
            if isinstance(params, dict):
                for key, expected_value in search_contract["params"].items():
                    observed_value = params.get(key)
                    if isinstance(expected_value, (int, float)):
                        try:
                            observed_value = _strict_number(observed_value)
                        except ValueError:
                            observed_value = math.nan
                        matches_parameter = math.isclose(
                            observed_value, float(expected_value),
                            rel_tol=1e-12, abs_tol=0.0,
                        )
                    else:
                        matches_parameter = observed_value == expected_value
                    if not matches_parameter:
                        structural_errors.append(
                            f"report {index} parameter mismatch for {key}: "
                            f"observed={observed_value!r}, expected={expected_value!r}"
                        )
            title = str(search.get("query_title", "")).strip()
            qid = title.split()[0] if title else ""
            try:
                qlen = _strict_integer(search.get("query_len"))
            except ValueError:
                qlen = 0
            if not qid:
                structural_errors.append(f"report {index} has no query title")
                continue
            if qid in observed:
                structural_errors.append(f"duplicate query report: {qid}")
                continue
            observed[qid] = qlen
            stat = search.get("stat") or {}
            if not isinstance(stat, dict):
                structural_errors.append(f"{qid} has a malformed statistics object")
                stat = {}
            database_values: dict[str, int] = {}
            for key in ("db_num", "db_len"):
                try:
                    value = _strict_integer(stat.get(key))
                except ValueError:
                    value = 0
                database_values[key] = value
                if value <= 0:
                    structural_errors.append(
                        f"{qid} has invalid database statistic: "
                        f"{key}={stat.get(key)!r}"
                    )
            if database_values["db_num"] > 0 and database_values["db_len"] > 0:
                database_signatures.append(
                    (database_values["db_num"], database_values["db_len"])
                )
            bad_stats = []
            for key in ("kappa", "lambda", "entropy"):
                try:
                    value = _strict_number(stat.get(key))
                except ValueError:
                    value = 0.0
                if not math.isfinite(value) or value <= 0:
                    bad_stats.append(f"{key}={stat.get(key)!r}")
            if bad_stats:
                structural_errors.append(
                    f'{qid} has invalid result statistics: {", ".join(bad_stats)}'
                )
            hits = search.get("hits", [])
            if not isinstance(hits, list):
                structural_errors.append(f"{qid} has a malformed hit list")
                hits = []
            if not hits and search.get("message") != "No hits found":
                structural_errors.append(
                    f"{qid} has no hits without the expected completion message"
                )
            if hits and search.get("message") not in (None, ""):
                structural_errors.append(
                    f"{qid} has hits together with a contradictory no-hit or error result message: "
                    f"{search.get('message')!r}"
                )
            if len(hits) > 100:
                structural_errors.append(
                    f"{qid} has {len(hits)} subjects despite max_target_seqs=100"
                )
            json_hit_counts[qid] = len(hits)
            hsp_count = 0
            for hit in hits:
                hsps = hit.get("hsps", []) if isinstance(hit, dict) else []
                accession_keys = _json_accession_keys(hit) if isinstance(hit, dict) else set()
                json_accessions, json_ids = (
                    _json_description_key_sets(hit) if isinstance(hit, dict)
                    else (set(), set())
                )
                description_records = (
                    _json_description_records(hit) if isinstance(hit, dict) else []
                )
                raw_descriptions = (
                    hit.get("description", []) if isinstance(hit, dict) else []
                )
                try:
                    subject_length = (
                        _strict_integer(hit.get("len")) if isinstance(hit, dict) else 0
                    )
                except ValueError:
                    subject_length = 0
                if (
                    not isinstance(raw_descriptions, list)
                    or not raw_descriptions
                    or len(description_records) != len(raw_descriptions)
                    or not accession_keys
                    or subject_length <= 0
                ):
                    structural_errors.append(
                        f"{qid} has a hit without valid subject identity/length"
                    )
                if not isinstance(hsps, list) or not hsps:
                    structural_errors.append(
                        f"{qid} has a hit without a nonempty HSP list"
                    )
                    continue
                if len(hsps) != 1:
                    structural_errors.append(
                        f"{qid} has {len(hsps)} HSPs for one subject despite max_hsps=1"
                    )
                for hsp in hsps:
                    hsp_count += 1
                    if not isinstance(hsp, dict):
                        structural_errors.append(f"{qid} has a malformed JSON HSP")
                        continue
                    try:
                        descriptions = hit.get("description", [])
                        raw_qseq, raw_sseq = hsp["qseq"], hsp["hseq"]
                        if (
                            not isinstance(raw_qseq, str)
                            or not isinstance(raw_sseq, str)
                            or not raw_qseq
                            or not raw_sseq
                        ):
                            raise ValueError("missing or non-string JSON HSP sequence")
                        record = {
                            "qseqid": qid,
                            "accession_keys": accession_keys,
                            "json_accessions": json_accessions,
                            "json_ids": json_ids,
                            "primary_id": (
                                description_records[0]["id"]
                                if description_records else ""
                            ),
                            "slen": subject_length,
                            "qstart": _strict_integer(hsp["query_from"]),
                            "qend": _strict_integer(hsp["query_to"]),
                            "sstart": _strict_integer(hsp["hit_from"]),
                            "send": _strict_integer(hsp["hit_to"]),
                            "length": _strict_integer(hsp["align_len"]),
                            "identity": _strict_integer(hsp["identity"]),
                            "evalue": _strict_number(hsp["evalue"]),
                            "bitscore": _strict_number(hsp["bit_score"]),
                            "qseq": raw_qseq.upper(),
                            "sseq": raw_sseq.upper(),
                            "descriptions": description_records,
                        }
                        if (
                            record["length"] <= 0
                            or not 0 <= record["identity"] <= record["length"]
                            or not math.isfinite(record["evalue"])
                            or record["evalue"] < 0
                            or record["evalue"] > search_contract["expect"]
                            or not math.isfinite(record["bitscore"])
                            or record["bitscore"] < 0
                        ):
                            raise ValueError("invalid JSON HSP numeric range")
                        json_hsp_records.append(record)
                    except (KeyError, TypeError, ValueError) as exc:
                        structural_errors.append(
                            f"{qid} has an incomplete JSON HSP: {exc}"
                        )
            json_hsp_counts[qid] = hsp_count

    if request_signatures and len(set(request_signatures)) != 1:
        structural_errors.append(
            "reports do not share one program/version/reference/database/parameter signature"
        )
    if database_signatures and len(set(database_signatures)) != 1:
        structural_errors.append(
            "reports do not share one database-size signature"
        )

    if observed != expected_lengths:
        structural_errors.append(
            f"query reports mismatch: observed={observed}, "
            f"expected={expected_lengths}"
        )

    hit_rows: list[dict[str, Any]] = []
    try:
        raw_hit_lines = hits_path.read_text(errors="replace").splitlines()
    except OSError as exc:
        structural_errors.append(f"unreadable formatted hit table: {exc}")
        raw_hit_lines = []
    for raw in raw_hit_lines:
        values = raw.split("\t")
        if len(values) != len(HIT_FIELDS):
            structural_errors.append(
                f"malformed BLAST result row with {len(values)} fields"
            )
            continue
        row: dict[str, Any] = dict(zip(HIT_FIELDS, values))
        if row["qseqid"] not in expected_lengths:
            structural_errors.append(
                f'unexpected query ID in BLAST result row: {row["qseqid"]}'
            )
            continue
        try:
            row["pident"] = float(row["pident"])
            row["qcovs"] = float(row["qcovs"])
            row["evalue"] = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            row_qlen = int(row["qlen"])
            row_slen = int(row["slen"])
            alignment_length = int(row["length"])
            qstart, qend = int(row["qstart"]), int(row["qend"])
            sstart, send = int(row["sstart"]), int(row["send"])
        except ValueError:
            structural_errors.append("malformed numeric value in BLAST result row")
            continue
        if (
            not math.isfinite(row["pident"])
            or not 0 <= row["pident"] <= 100
            or not math.isfinite(row["qcovs"])
            or not 0 <= row["qcovs"] <= 100
            or not math.isfinite(row["evalue"])
            or row["evalue"] < 0
            or row["evalue"] > search_contract["expect"]
            or not math.isfinite(row["bitscore"])
            or row["bitscore"] < 0
            or row_qlen <= 0
            or row_slen <= 0
            or alignment_length <= 0
        ):
            structural_errors.append(
                f'invalid numeric range in BLAST result row for {row["qseqid"]}'
            )
        if row_qlen != expected_lengths[row["qseqid"]]:
            structural_errors.append(
                f'query length mismatch in BLAST result row: {row["qseqid"]} '
                f"observed={row_qlen} expected={expected_lengths[row['qseqid']]}"
            )
        if not (1 <= min(qstart, qend) <= max(qstart, qend) <= row_qlen):
            structural_errors.append(
                f'query coordinates outside query bounds for {row["qseqid"]}'
            )
        if not (1 <= min(sstart, send) <= max(sstart, send) <= row_slen):
            structural_errors.append(
                f'subject coordinates outside subject bounds for {row["qseqid"]}'
            )
        qseq, sseq = str(row["qseq"]), str(row["sseq"])
        alignment_pattern = (
            r"[ACGTURYKMSWBDHVN-]+" if mode.startswith("nt_")
            else r"[A-Z*-]+"
        )
        if (
            not qseq
            or not sseq
            or len(qseq) != len(sseq)
            or len(qseq) != alignment_length
            or re.fullmatch(alignment_pattern, qseq.upper()) is None
            or re.fullmatch(alignment_pattern, sseq.upper()) is None
        ):
            structural_errors.append(
                f'malformed aligned query/subject sequence strings for {row["qseqid"]}'
            )
        else:
            if any(qchar == "-" and schar == "-" for qchar, schar in zip(qseq, sseq)):
                structural_errors.append(
                    f'double-gap alignment column for {row["qseqid"]}'
                )
            ungapped_query = len(qseq.replace("-", ""))
            ungapped_subject = len(sseq.replace("-", ""))
            if abs(qend - qstart) + 1 != ungapped_query:
                structural_errors.append(
                    f'query endpoint/alignment mismatch for {row["qseqid"]}'
                )
            else:
                source_query = query_sequences[row["qseqid"]]
                if qstart <= qend:
                    expected_query_fragment = source_query[qstart - 1:qend]
                elif mode.startswith("nt_"):
                    forward = source_query[qend - 1:qstart]
                    expected_query_fragment = forward.translate(
                        str.maketrans("ACGTN", "TGCAN")
                    )[::-1]
                else:
                    expected_query_fragment = ""
                if qseq.replace("-", "").upper() != expected_query_fragment:
                    structural_errors.append(
                        f'query sequence/archive mismatch for {row["qseqid"]}'
                    )
            if abs(send - sstart) + 1 != ungapped_subject:
                structural_errors.append(
                    f'subject endpoint/alignment mismatch for {row["qseqid"]}'
                )
            identities = sum(
                qchar.upper() == schar.upper() and qchar != "-"
                for qchar, schar in zip(qseq, sseq)
            )
            derived_pident = 100.0 * identities / alignment_length
            if abs(row["pident"] - derived_pident) > 0.01:
                structural_errors.append(
                    f'identity/alignment mismatch for {row["qseqid"]}: '
                    f'observed={row["pident"]}, derived={derived_pident:.6f}'
                )
            # The search command fixes max_hsps=1, so qcovs is the coverage of
            # this single HSP for its subject (reported by BLAST as an integer).
            derived_qcov = 100.0 * ungapped_query / row_qlen
            row["_derived_qcov"] = derived_qcov
            expected_qcov = math.floor(derived_qcov + 0.5)
            if not row["qcovs"].is_integer() or int(row["qcovs"]) != expected_qcov:
                structural_errors.append(
                    f'query-coverage/alignment mismatch for {row["qseqid"]}: '
                    f'observed={row["qcovs"]}, expected={expected_qcov}'
                )
        hit_rows.append(row)

    tsv_hit_counts = Counter(row["qseqid"] for row in hit_rows)
    for qid in expected_lengths:
        # BLAST tabular output has one row per HSP, whereas JSON ``hits`` has
        # one item per subject and can contain multiple HSPs.  Compare like
        # with like so a valid multi-HSP subject is not rejected.
        if json_hsp_counts.get(qid, 0) != tsv_hit_counts.get(qid, 0):
            structural_errors.append(
                f"JSON/TSV HSP-count mismatch for {qid}: "
                f"json={json_hsp_counts.get(qid, 0)}, "
                f"tsv={tsv_hit_counts.get(qid, 0)}"
            )

    for qid in expected_lengths:
        tsv_primary = [
            _accession_token(str(row["saccver"]))
            for row in hit_rows if row["qseqid"] == qid
        ]
        if len(tsv_primary) != len(set(tsv_primary)):
            structural_errors.append(
                f"duplicate primary subject accession in TSV for {qid}"
            )
        json_primary = [
            record["primary_id"]
            for record in json_hsp_records if record["qseqid"] == qid
        ]
        if len(json_primary) != len(set(json_primary)):
            structural_errors.append(
                f"duplicate primary subject accession in JSON for {qid}"
            )

    unused_json_hsps = set(range(len(json_hsp_records)))
    for row in hit_rows:
        row_accessions = _tsv_accession_keys(row)
        row_sallacc = _tsv_sallacc_keys(row)
        row_sallseqid = _tsv_sallseqid_keys(row)
        row_primary = _accession_token(str(row["saccver"]))
        row_identity = sum(
            qchar.upper() == schar.upper() and qchar != "-"
            for qchar, schar in zip(str(row["qseq"]), str(row["sseq"]))
        )
        row_taxids = {value for value in str(row["staxids"]).split(";") if value and value != "N/A"}
        matches = [
            index for index in unused_json_hsps
            if json_hsp_records[index]["qseqid"] == row["qseqid"]
            and row_primary == json_hsp_records[index]["primary_id"]
            and row_sallacc == json_hsp_records[index]["json_accessions"]
            and row_sallseqid == json_hsp_records[index]["json_ids"]
            and json_hsp_records[index]["slen"] == int(row["slen"])
            and json_hsp_records[index]["qstart"] == int(row["qstart"])
            and json_hsp_records[index]["qend"] == int(row["qend"])
            and json_hsp_records[index]["sstart"] == int(row["sstart"])
            and json_hsp_records[index]["send"] == int(row["send"])
            and json_hsp_records[index]["length"] == int(row["length"])
            and json_hsp_records[index]["identity"] == row_identity
            and math.isclose(
                json_hsp_records[index]["evalue"], float(row["evalue"]),
                rel_tol=0.02, abs_tol=0.0,
            )
            and abs(json_hsp_records[index]["bitscore"] - float(row["bitscore"])) <= 1.01
            and json_hsp_records[index]["qseq"] == str(row["qseq"]).upper()
            and json_hsp_records[index]["sseq"] == str(row["sseq"]).upper()
            and all((
                _expected_accession_present(
                    str(row["saccver"]), description["accession_keys"]
                )
                and description["title"] == str(row["stitle"])
                and (
                    bool(description["taxids"] & row_taxids)
                    if description["taxids"] or row_taxids
                    else True
                )
                for description in json_hsp_records[index]["descriptions"][:1]
            ))
            and row_taxids == {
                taxid
                for description in json_hsp_records[index]["descriptions"]
                for taxid in description["taxids"]
            }
        ]
        if not matches:
            structural_errors.append(
                f'JSON/TSV HSP identity mismatch for {row["qseqid"]}: '
                f'accessions={sorted(row_accessions)}'
            )
        else:
            matched_index = matches[0]
            row["_json_accession_keys"] = json_hsp_records[matched_index]["accession_keys"]
            unused_json_hsps.remove(matched_index)
    if unused_json_hsps:
        structural_errors.append(
            f"JSON contains {len(unused_json_hsps)} HSP(s) absent from the TSV binding"
        )

    # Structural integrity takes priority.  A missing control in an empty or
    # zero-statistic archive is not a deterministic control failure.
    if not structural_errors:
        for spec in control_specs:
            control = spec["id"]
            if json_hit_counts.get(control, 0) < 1:
                control_errors.append(f"positive control has no hit: {control}")
                continue
            accession = spec.get("expected_accession")
            if accession:
                matches = [
                    row
                    for row in hit_rows
                    if row["qseqid"] == control
                    and _expected_accession_present(
                        accession, _tsv_accession_keys(row)
                    )
                    and _expected_accession_present(
                        accession, set(row.get("_json_accession_keys", set()))
                    )
                    and row["pident"] >= float(spec["min_identity"])
                    and float(row.get("_derived_qcov", 0.0)) >= float(
                        spec["min_query_coverage"]
                    )
                ]
                if not matches:
                    control_errors.append(
                        "positive control did not recover a near-exact match: "
                        f"{control}"
                    )

    return structural_errors, control_errors, len(observed)


def main(argv: list[str]) -> int:
    if len(argv) != 5:
        print(
            "usage: validate_panax_remote_archive.py "
            "RESULTS.json EXPECTED_QUERIES.json MODE HITS.tsv",
            file=sys.stderr,
        )
        return 2
    structural, control, observed_count = validate(
        Path(argv[1]), Path(argv[2]), argv[3], Path(argv[4])
    )
    if structural:
        print("remote archive structural validation failed:", file=sys.stderr)
        for error in structural:
            print(f"- {error}", file=sys.stderr)
        return EXIT_STRUCTURAL
    if control:
        print("remote archive control validation failed:", file=sys.stderr)
        for error in control:
            print(f"- {error}", file=sys.stderr)
        return EXIT_CONTROL
    print(f"validated {observed_count} query reports with nonzero statistics")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
