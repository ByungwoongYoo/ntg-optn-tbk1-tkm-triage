#!/usr/bin/env python3
"""Aggregate six duplicate-aware proper-fragment Panax audits.

All thresholds are workflow-defined technical rules loaded from a hashable
JSON file. Exact fragment sharing is a batch-QC indicator; even a
``shadow_suspect`` result neither proves index hopping nor identifies its
mechanism. No output establishes a virus species, biological host, infection,
replication, phenotype association, pathogenicity, or transmission.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_RUNS = tuple(f"DRR8539{number:02d}" for number in range(7, 13))
ALIASES = dict(zip(EXPECTED_RUNS, ("ZhA", "ZhB", "ZhC", "BhA", "BhB", "BhC")))
REFERENCES = ("PNX_Picorna_A1_ref", "PNX_Picorna_A2_ref", "PNX_Picorna_B_ref")
HEX64 = re.compile(r"[0-9a-f]{64}")
CLAIM_BOUNDARY = (
    "Technical support for predefined partial Picornavirales-like sequence "
    "clusters only; no formal taxon, true host, infection, replication, "
    "phenotype association, pathogenicity, or transmission claim."
)
# All depth summaries are emitted to six decimal places.  In the worst case,
# independent rounding of the mean and the three weighted breadth terms can
# create a 5.5e-6 apparent violation of the depth lower bound.
DEPTH_ROUNDING_TOLERANCE = 6e-6


class AggregateError(ValueError):
    """Structural or semantic input failure."""


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: object, label: str, low: float = 0.0, high: float = math.inf) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise AggregateError(f"invalid_number:{label}:{value!r}") from exc
    if not math.isfinite(parsed) or not low <= parsed <= high:
        raise AggregateError(f"number_out_of_range:{label}:{parsed}")
    return parsed


def integer(value: object, label: str, low: int = 0) -> int:
    parsed = number(value, label, low)
    if not parsed.is_integer():
        raise AggregateError(f"not_an_integer:{label}:{parsed}")
    return int(parsed)


def load_thresholds(path: Path) -> dict[str, object]:
    try:
        raw = json.loads(path.read_text())
    except Exception as exc:
        raise AggregateError(f"invalid_threshold_json:{type(exc).__name__}") from exc
    if not isinstance(raw, dict):
        raise AggregateError("thresholds_must_be_object")
    selection = raw.get("fragment_selection")
    validated = raw.get("fragment_validated_positive")
    strong = raw.get("fragment_strong_positive")
    shadow = raw.get("cross_run_exact_fragment_shadow_qc")
    candidate = raw.get("candidate_fragment_support")
    structure = raw.get("candidate_structure_contract")
    if not all(
        isinstance(section, dict)
        for section in (selection, validated, strong, shadow, candidate, structure)
    ):
        raise AggregateError("missing_fragment_threshold_section")
    minimum_mapq = integer(selection.get("minimum_mapq_each_mate"), "minimum_mapq", 0)
    if minimum_mapq > 255:
        raise AggregateError("minimum_mapq_above_255")
    expected_selection = {
        "both_mates_required": True,
        "same_candidate_reference_required": True,
        "proper_pair_required": True,
        "primary_alignment_required": True,
        "qc_fail_excluded": True,
    }
    for field, expected in expected_selection.items():
        if selection.get(field) is not expected:
            raise AggregateError(f"invalid_fragment_selection:{field}")
    structure_filename = structure.get("filename")
    structure_sha = structure.get("sha256")
    structure_lengths = structure.get("reference_lengths_nt")
    if (
        not isinstance(structure_filename, str)
        or not structure_filename
        or Path(structure_filename).name != structure_filename
    ):
        raise AggregateError("invalid_candidate_structure_filename")
    if not isinstance(structure_sha, str) or not HEX64.fullmatch(structure_sha):
        raise AggregateError("invalid_candidate_structure_sha256")
    if not isinstance(structure_lengths, dict) or set(structure_lengths) != set(REFERENCES):
        raise AggregateError("invalid_candidate_structure_reference_set")
    pinned_lengths = {
        reference: integer(
            structure_lengths[reference], f"candidate_structure_length:{reference}", 1
        )
        for reference in REFERENCES
    }
    result: dict[str, object] = {
        "minimum_mapq": minimum_mapq,
        "validated_preduplicate_fragments": integer(
            validated.get("proper_fragments_preduplicate_minimum"), "validated_predup"
        ),
        "validated_preduplicate_breadth_1x": number(
            validated.get("proper_preduplicate_breadth_1x_minimum"), "validated_predup_b1", 0, 1
        ),
        "validated_nonduplicate_fragments": integer(
            validated.get("proper_fragments_nonduplicate_minimum"), "validated_nondup"
        ),
        "validated_nonduplicate_breadth_1x": number(
            validated.get("proper_nonduplicate_breadth_1x_minimum"), "validated_nondup_b1", 0, 1
        ),
        "strong_preduplicate_fragments": integer(
            strong.get("proper_fragments_preduplicate_minimum"), "strong_predup"
        ),
        "strong_preduplicate_breadth_1x": number(
            strong.get("proper_preduplicate_breadth_1x_minimum"), "strong_predup_b1", 0, 1
        ),
        "strong_preduplicate_breadth_5x": number(
            strong.get("proper_preduplicate_breadth_5x_minimum"), "strong_predup_b5", 0, 1
        ),
        "strong_nonduplicate_fragments": integer(
            strong.get("proper_fragments_nonduplicate_minimum"), "strong_nondup"
        ),
        "strong_nonduplicate_breadth_1x": number(
            strong.get("proper_nonduplicate_breadth_1x_minimum"), "strong_nondup_b1", 0, 1
        ),
        "shadow_min_unique": integer(
            shadow.get("minimum_unique_fingerprints_in_smaller_run"), "shadow_min_unique"
        ),
        "shadow_warning_shared": integer(
            shadow.get("warning_shared_combined_minimum"), "shadow_warning_shared"
        ),
        "shadow_warning_containment": number(
            shadow.get("warning_smaller_run_containment_minimum"), "shadow_warning_containment", 0, 1
        ),
        "shadow_suspect_shared": integer(
            shadow.get("suspect_shared_combined_minimum"), "shadow_suspect_shared"
        ),
        "shadow_suspect_containment": number(
            shadow.get("suspect_smaller_run_containment_minimum"), "shadow_suspect_containment", 0, 1
        ),
        "shadow_abundance_ratio": number(
            shadow.get("larger_to_smaller_abundance_ratio_minimum"), "shadow_ratio", 1
        ),
        "candidate_validated_runs": integer(
            candidate.get("minimum_independent_validated_runs"), "candidate_validated_runs", 1
        ),
        "candidate_strong_runs": integer(
            candidate.get("minimum_independent_strong_runs"), "candidate_strong_runs", 1
        ),
        "candidate_structure_filename": structure_filename,
        "candidate_structure_sha256": structure_sha,
        "candidate_reference_lengths": pinned_lengths,
        "scope": raw.get("scope", ""),
        "source_file_sha256": sha256_path(path),
    }
    for strong_key, validated_key in (
        ("strong_preduplicate_fragments", "validated_preduplicate_fragments"),
        ("strong_preduplicate_breadth_1x", "validated_preduplicate_breadth_1x"),
        ("strong_nonduplicate_fragments", "validated_nonduplicate_fragments"),
        ("strong_nonduplicate_breadth_1x", "validated_nonduplicate_breadth_1x"),
    ):
        if result[strong_key] < result[validated_key]:
            raise AggregateError(f"strong_threshold_below_validated:{strong_key}")
    if result["shadow_suspect_shared"] < result["shadow_warning_shared"]:
        raise AggregateError("shadow_suspect_shared_below_warning")
    if result["shadow_suspect_containment"] < result["shadow_warning_containment"]:
        raise AggregateError("shadow_suspect_containment_below_warning")
    if not 1 <= result["candidate_strong_runs"] <= result["candidate_validated_runs"] <= len(EXPECTED_RUNS):
        raise AggregateError("invalid_candidate_run_thresholds")
    return result


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise AggregateError(f"missing_tsv_header:{path}")
            fields = list(reader.fieldnames)
            if (
                any(not isinstance(field, str) or not field or field != field.strip() for field in fields)
                or len(fields) != len(set(fields))
            ):
                raise AggregateError(f"invalid_or_duplicate_tsv_header:{path}:{fields}")
            rows: list[dict[str, str]] = []
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise AggregateError(f"malformed_tsv_row:{path}:line_{reader.line_num}")
                rows.append(row)
            return fields, rows
    except csv.Error as exc:
        raise AggregateError(f"invalid_tsv:{path}:{exc}") from exc


def load_candidate_structure(
    path: Path, thresholds: dict[str, object]
) -> dict[str, int]:
    expected_sha = str(thresholds["candidate_structure_sha256"])
    observed_sha = sha256_path(path)
    if observed_sha != expected_sha:
        raise AggregateError(
            f"candidate_structure_sha256_mismatch:{observed_sha}!={expected_sha}"
        )
    _, rows = read_tsv(path)
    if len(rows) != len(REFERENCES):
        raise AggregateError(f"candidate_structure_row_count:{len(rows)}")
    required = {"reference", "reference_length_nt", "reference_sha256"}
    if any(not required.issubset(row) for row in rows):
        raise AggregateError("candidate_structure_missing_fields")
    by_reference = {row["reference"]: row for row in rows}
    if len(by_reference) != len(rows) or set(by_reference) != set(REFERENCES):
        raise AggregateError(
            f"candidate_structure_reference_set:{sorted(by_reference)}"
        )
    pinned = thresholds["candidate_reference_lengths"]
    assert isinstance(pinned, dict)
    lengths: dict[str, int] = {}
    for reference in REFERENCES:
        row = by_reference[reference]
        length = integer(
            row["reference_length_nt"], f"candidate_structure_length:{reference}", 1
        )
        if length != pinned[reference]:
            raise AggregateError(
                f"candidate_structure_length_mismatch:{reference}:"
                f"{length}!={pinned[reference]}"
            )
        if not HEX64.fullmatch(row["reference_sha256"]):
            raise AggregateError(f"candidate_structure_reference_sha256:{reference}")
        lengths[reference] = length
    return lengths


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_metric_row(
    row: dict[str, str], expected_run: str, thresholds: dict[str, object]
) -> dict[str, object]:
    run, alias, reference = row.get("run", ""), row.get("alias", ""), row.get("reference", "")
    if run != expected_run or alias != ALIASES.get(run) or reference not in REFERENCES:
        raise AggregateError(f"metric_identity_mismatch:{expected_run}:{run}:{alias}:{reference}")
    mapq = int(thresholds["minimum_mapq"])
    if integer(row.get("minimum_mate_mapq"), f"{run}:{reference}:mapq") != mapq:
        raise AggregateError(f"metric_mapq_mismatch:{run}:{reference}")
    pre_field = f"proper_fragments_mapq{mapq}_preduplicate"
    non_field = f"proper_fragments_mapq{mapq}_nonduplicate"
    k = 31
    required = {
        "reference_length_nt", pre_field, non_field,
        "duplicate_fragment_fraction",
        "panel_unique_kmer_fragments_nonduplicate",
        "panel_unique_kmer_fragment_fraction",
        f"panel_unique_{k}mers_available", f"panel_unique_{k}mers_observed",
        f"panel_unique_{k}mer_fraction_observed", "distinct_fragment_endpoints",
        "fr_orientation_fraction", "template_span_median", "template_span_p10",
        "template_span_p90", "softclipped_fragment_fraction",
        "max_softclipped_bases_per_fragment",
        *{
            f"proper_{duplicate}_{metric}"
            for duplicate in ("preduplicate", "nonduplicate")
            for metric in (
                "breadth_1x", "breadth_5x", "breadth_10x", "mean_depth",
                "median_depth", "max_zero_run", "max_internal_zero_run",
            )
        },
    }
    missing = sorted(required - set(row))
    if missing:
        raise AggregateError(f"missing_metric_fields:{run}:{reference}:{missing}")
    pre = integer(row[pre_field], f"{run}:{reference}:predup")
    non = integer(row[non_field], f"{run}:{reference}:nondup")
    reference_length = integer(
        row["reference_length_nt"], f"{run}:{reference}:reference_length", 1
    )
    if non > pre:
        raise AggregateError(f"nonduplicate_exceeds_preduplicate:{run}:{reference}")
    expected_fraction = "" if not pre else (pre - non) / pre
    observed_fraction = row["duplicate_fragment_fraction"]
    if expected_fraction == "":
        if observed_fraction != "":
            raise AggregateError(f"unexpected_duplicate_fraction:{run}:{reference}")
    elif not math.isclose(number(observed_fraction, f"{run}:{reference}:dupfrac", 0, 1), expected_fraction, abs_tol=1e-6):
        raise AggregateError(f"duplicate_fraction_mismatch:{run}:{reference}")
    diagnostic = integer(row["panel_unique_kmer_fragments_nonduplicate"], f"{run}:{reference}:panel_kmer")
    endpoints = integer(row["distinct_fragment_endpoints"], f"{run}:{reference}:endpoints")
    if diagnostic > non or endpoints > non:
        raise AggregateError(f"fragment_subcount_exceeds_nonduplicate:{run}:{reference}")
    for field, numerator, denominator in (
        ("panel_unique_kmer_fragment_fraction", diagnostic, non),
    ):
        observed = row[field]
        if denominator == 0:
            if observed != "":
                raise AggregateError(f"unexpected_empty_denominator_fraction:{run}:{reference}:{field}")
        elif not math.isclose(
            number(observed, f"{run}:{reference}:{field}", 0, 1),
            numerator / denominator, abs_tol=1e-6,
        ):
            raise AggregateError(f"fraction_mismatch:{run}:{reference}:{field}")
    available = integer(row[f"panel_unique_{k}mers_available"], f"{run}:{reference}:kmers_available")
    observed_kmers = integer(row[f"panel_unique_{k}mers_observed"], f"{run}:{reference}:kmers_observed")
    if observed_kmers > available:
        raise AggregateError(f"observed_kmers_exceed_available:{run}:{reference}")
    if (diagnostic == 0) != (observed_kmers == 0):
        raise AggregateError(
            f"diagnostic_fragment_kmer_mismatch:{run}:{reference}:"
            f"fragments={diagnostic}:observed_kmers={observed_kmers}"
        )
    if non == 0 and (diagnostic != 0 or observed_kmers != 0):
        raise AggregateError(f"zero_fragment_kmer_support:{run}:{reference}")
    kmer_fraction = row[f"panel_unique_{k}mer_fraction_observed"]
    if available == 0:
        if kmer_fraction != "":
            raise AggregateError(f"unexpected_kmer_fraction:{run}:{reference}")
    elif not math.isclose(
        number(kmer_fraction, f"{run}:{reference}:kmer_fraction", 0, 1),
        observed_kmers / available, abs_tol=1e-6,
    ):
        raise AggregateError(f"kmer_fraction_mismatch:{run}:{reference}")

    values: dict[str, float] = {}
    zero_runs: dict[str, int] = {}
    for duplicate in ("preduplicate", "nonduplicate"):
        breadths = []
        for depth in (1, 5, 10):
            field = f"proper_{duplicate}_breadth_{depth}x"
            values[field] = number(row[field], f"{run}:{reference}:{field}", 0, 1)
            breadths.append(values[field])
        if not breadths[0] + 1e-6 >= breadths[1] >= breadths[2] - 1e-6:
            raise AggregateError(f"nonmonotonic_breadth:{run}:{reference}:{duplicate}")
        for metric in ("mean_depth", "median_depth"):
            field = f"proper_{duplicate}_{metric}"
            values[field] = number(row[field], f"{run}:{reference}:{field}")
        for metric in ("max_zero_run", "max_internal_zero_run"):
            field = f"proper_{duplicate}_{metric}"
            zero_runs[field] = integer(row[field], f"{run}:{reference}:{field}", 0)
            if zero_runs[field] > reference_length:
                raise AggregateError(f"zero_run_exceeds_reference:{run}:{reference}:{field}")
        mean_field = f"proper_{duplicate}_mean_depth"
        minimum_mean_depth = breadths[0] + 4 * breadths[1] + 5 * breadths[2]
        if values[mean_field] + DEPTH_ROUNDING_TOLERANCE < minimum_mean_depth:
            raise AggregateError(
                f"mean_depth_below_breadth_lower_bound:{run}:{reference}:{duplicate}:"
                f"{values[mean_field]}<{minimum_mean_depth}"
            )
        max_zero_field = f"proper_{duplicate}_max_zero_run"
        max_internal_zero_field = f"proper_{duplicate}_max_internal_zero_run"
        max_zero = zero_runs[max_zero_field]
        if breadths[0] == 0:
            if max_zero != reference_length:
                raise AggregateError(
                    f"zero_breadth_max_zero_run_mismatch:{run}:{reference}:{duplicate}"
                )
        elif breadths[0] < 1:
            if max_zero < 1:
                raise AggregateError(
                    f"partial_breadth_requires_zero_run:{run}:{reference}:{duplicate}"
                )
        elif max_zero != 0:
            raise AggregateError(
                f"complete_breadth_forbids_zero_run:{run}:{reference}:{duplicate}"
            )
        if zero_runs[max_internal_zero_field] > max_zero:
            raise AggregateError(
                f"internal_zero_run_exceeds_maximum:{run}:{reference}:{duplicate}"
            )
        fragment_count = pre if duplicate == "preduplicate" else non
        if fragment_count == 0:
            zero_expected = [
                values[f"proper_{duplicate}_breadth_{depth}x"]
                for depth in (1, 5, 10)
            ] + [
                values[f"proper_{duplicate}_mean_depth"],
                values[f"proper_{duplicate}_median_depth"],
            ]
            if any(abs(value) > 1e-9 for value in zero_expected):
                raise AggregateError(
                    f"zero_fragment_nonzero_depth_metric:{run}:{reference}:{duplicate}"
                )
            if zero_runs[f"proper_{duplicate}_max_zero_run"] != reference_length:
                raise AggregateError(
                    f"zero_fragment_max_zero_run:{run}:{reference}:{duplicate}"
                )
            if zero_runs[f"proper_{duplicate}_max_internal_zero_run"] != 0:
                raise AggregateError(
                    f"zero_fragment_internal_zero_run:{run}:{reference}:{duplicate}"
                )
    for depth in (1, 5, 10):
        non_breadth = values[f"proper_nonduplicate_breadth_{depth}x"]
        pre_breadth = values[f"proper_preduplicate_breadth_{depth}x"]
        if non_breadth > pre_breadth + 1e-6:
            raise AggregateError(
                f"nonduplicate_breadth_exceeds_preduplicate:{run}:{reference}:{depth}x"
            )
    if values["proper_nonduplicate_mean_depth"] > values["proper_preduplicate_mean_depth"] + 1e-6:
        raise AggregateError(f"nonduplicate_depth_exceeds_preduplicate:{run}:{reference}")
    if values["proper_nonduplicate_median_depth"] > values["proper_preduplicate_median_depth"] + 1e-6:
        raise AggregateError(f"nonduplicate_median_exceeds_preduplicate:{run}:{reference}")
    if integer(
        row["proper_nonduplicate_max_zero_run"],
        f"{run}:{reference}:nonduplicate_max_zero_run",
    ) < integer(
        row["proper_preduplicate_max_zero_run"],
        f"{run}:{reference}:preduplicate_max_zero_run",
    ):
        raise AggregateError(f"nonduplicate_zero_run_below_preduplicate:{run}:{reference}")

    for field in ("fr_orientation_fraction", "softclipped_fragment_fraction"):
        observed = row[field]
        if non == 0:
            if observed != "":
                raise AggregateError(f"unexpected_fragment_fraction:{run}:{reference}:{field}")
        else:
            number(observed, f"{run}:{reference}:{field}", 0, 1)
    spans = [row[field] for field in ("template_span_p10", "template_span_median", "template_span_p90")]
    if non == 0:
        if any(value != "" for value in spans):
            raise AggregateError(f"unexpected_template_span:{run}:{reference}")
    else:
        parsed_spans = [number(value, f"{run}:{reference}:template_span", 1, reference_length) for value in spans]
        if not parsed_spans[0] <= parsed_spans[1] <= parsed_spans[2]:
            raise AggregateError(f"nonmonotonic_template_span:{run}:{reference}")
    max_softclip = integer(
        row["max_softclipped_bases_per_fragment"],
        f"{run}:{reference}:max_softclipped_bases", 0,
    )
    if non == 0 and max_softclip != 0:
        raise AggregateError(f"zero_fragment_softclip_maximum:{run}:{reference}")

    pre_validated = bool(
        pre >= thresholds["validated_preduplicate_fragments"]
        and values["proper_preduplicate_breadth_1x"] >= thresholds["validated_preduplicate_breadth_1x"]
    )
    non_validated = bool(
        non >= thresholds["validated_nonduplicate_fragments"]
        and values["proper_nonduplicate_breadth_1x"] >= thresholds["validated_nonduplicate_breadth_1x"]
    )
    validated = bool(
        pre_validated and non_validated
    )
    strong = bool(
        pre >= thresholds["strong_preduplicate_fragments"]
        and values["proper_preduplicate_breadth_1x"] >= thresholds["strong_preduplicate_breadth_1x"]
        and values["proper_preduplicate_breadth_5x"] >= thresholds["strong_preduplicate_breadth_5x"]
        and non >= thresholds["strong_nonduplicate_fragments"]
        and values["proper_nonduplicate_breadth_1x"] >= thresholds["strong_nonduplicate_breadth_1x"]
    )
    output: dict[str, object] = dict(row)
    output["fragment_preduplicate_validated_positive"] = str(pre_validated).lower()
    output["fragment_nonduplicate_validated_positive"] = str(non_validated).lower()
    output["fragment_validated_positive"] = str(validated).lower()
    output["fragment_strong_positive"] = str(strong).lower()
    return output


def validate_hash_rows(
    rows: list[dict[str, str]], run: str, metrics: dict[str, dict[str, object]], mapq: int
) -> list[dict[str, object]]:
    required = {
        "run", "alias", "reference", "minimum_mate_mapq", "pair_sequence_sha256",
        "endpoint_sha256", "fragment_fingerprint_sha256",
    }
    counts = {reference: 0 for reference in REFERENCES}
    endpoint_hashes = {reference: set() for reference in REFERENCES}
    output: list[dict[str, object]] = []
    for row in rows:
        if not required.issubset(row):
            raise AggregateError(f"missing_hash_fields:{run}")
        reference = row.get("reference", "")
        if row.get("run") != run or row.get("alias") != ALIASES[run] or reference not in REFERENCES:
            raise AggregateError(f"hash_identity_mismatch:{run}:{reference}")
        if integer(row.get("minimum_mate_mapq"), f"{run}:{reference}:hash_mapq") != mapq:
            raise AggregateError(f"hash_mapq_mismatch:{run}:{reference}")
        pair_hash, endpoint_hash, combined = (
            row.get("pair_sequence_sha256", ""), row.get("endpoint_sha256", ""),
            row.get("fragment_fingerprint_sha256", ""),
        )
        if not all(HEX64.fullmatch(value) for value in (pair_hash, endpoint_hash, combined)):
            raise AggregateError(f"invalid_fragment_hash:{run}:{reference}")
        expected = hashlib.sha256(
            (reference + "\n" + pair_hash + "\n" + endpoint_hash).encode()
        ).hexdigest()
        if combined != expected:
            raise AggregateError(f"combined_fragment_hash_mismatch:{run}:{reference}")
        counts[reference] += 1
        endpoint_hashes[reference].add(endpoint_hash)
        output.append(dict(row))
    non_field = f"proper_fragments_mapq{mapq}_nonduplicate"
    for reference in REFERENCES:
        expected_count = integer(
            metrics[reference][non_field], f"{run}:{reference}:nonduplicate_hash_count"
        )
        if counts[reference] != expected_count:
            raise AggregateError(
                f"hash_fragment_count_mismatch:{run}:{reference}:{counts[reference]}!={expected_count}"
            )
        expected_endpoints = integer(
            metrics[reference]["distinct_fragment_endpoints"],
            f"{run}:{reference}:distinct_fragment_endpoints",
        )
        if len(endpoint_hashes[reference]) != expected_endpoints:
            raise AggregateError(
                f"endpoint_hash_count_mismatch:{run}:{reference}:"
                f"{len(endpoint_hashes[reference])}!={expected_endpoints}"
            )
    return output


def load_inputs(
    root: Path, thresholds: dict[str, object], pinned_lengths: dict[str, int]
):
    tables = sorted(root.rglob("FRAGMENT_METRICS.tsv"))
    if len(tables) != len(EXPECTED_RUNS):
        raise AggregateError(f"fragment_metric_table_count:{len(tables)}!=6")
    mapq = int(thresholds["minimum_mapq"])
    by_run: dict[str, dict[str, dict[str, object]]] = {}
    all_metrics: list[dict[str, object]] = []
    all_hashes: list[dict[str, object]] = []
    all_continuity: list[dict[str, object]] = []
    metric_fields: list[str] | None = None
    hash_fields: list[str] | None = None
    continuity_fields: list[str] | None = None
    reference_lengths: dict[str, int] = {}
    panel_kmers_available: dict[str, set[int]] = {
        reference: set() for reference in REFERENCES
    }
    for table in tables:
        fields, raw_rows = read_tsv(table)
        if len(raw_rows) != len(REFERENCES):
            raise AggregateError(f"fragment_metric_row_count:{table}:{len(raw_rows)}")
        runs = {row.get("run", "") for row in raw_rows}
        if len(runs) != 1:
            raise AggregateError(f"mixed_metric_runs:{table}:{sorted(runs)}")
        run = next(iter(runs))
        if run not in EXPECTED_RUNS or run in by_run:
            raise AggregateError(f"unexpected_or_duplicate_metric_run:{run}")
        validated = [validate_metric_row(row, run, thresholds) for row in raw_rows]
        refs = [str(row["reference"]) for row in validated]
        if len(set(refs)) != len(refs) or set(refs) != set(REFERENCES):
            raise AggregateError(f"metric_reference_set:{run}:{refs}")
        for row in validated:
            reference = str(row["reference"])
            length = integer(
                row["reference_length_nt"], f"{run}:{reference}:metric_length", 1
            )
            if length != pinned_lengths[reference]:
                raise AggregateError(
                    f"metric_pinned_length_mismatch:{run}:{reference}:"
                    f"{length}!={pinned_lengths[reference]}"
                )
            panel_kmers_available[reference].add(
                integer(
                    row["panel_unique_31mers_available"],
                    f"{run}:{reference}:panel_unique_31mers_available",
                )
            )
        by_run[run] = {str(row["reference"]): row for row in validated}
        metric_fields = metric_fields or list(validated[0])
        if list(validated[0]) != metric_fields:
            raise AggregateError(f"metric_schema_mismatch:{run}")
        all_metrics.extend(validated)

        hash_path = table.parent / "FRAGMENT_HASHES.tsv"
        hfields, hrows = read_tsv(hash_path)
        if not set(("pair_sequence_sha256", "endpoint_sha256", "fragment_fingerprint_sha256")).issubset(hfields):
            raise AggregateError(f"hash_schema:{run}:{hfields}")
        hash_fields = hash_fields or hfields
        if hfields != hash_fields:
            raise AggregateError(f"hash_schema_mismatch:{run}")
        all_hashes.extend(validate_hash_rows(hrows, run, by_run[run], mapq))

        continuity_path = table.parent / "CONTINUITY_SCAN.tsv"
        cfields, crows = read_tsv(continuity_path)
        continuity_fields = continuity_fields or cfields
        if cfields != continuity_fields:
            raise AggregateError(f"continuity_schema_mismatch:{run}")
        seen: set[tuple[str, int]] = set()
        for row in crows:
            reference = row.get("reference", "")
            if row.get("run") != run or row.get("alias") != ALIASES[run] or reference not in REFERENCES:
                raise AggregateError(f"continuity_identity:{run}:{reference}")
            if integer(row.get("minimum_mate_mapq"), f"{run}:{reference}:continuity_mapq") != mapq:
                raise AggregateError(f"continuity_mapq:{run}:{reference}")
            length = integer(row.get("reference_length"), f"{run}:{reference}:length", 1)
            metric_length = integer(
                by_run[run][reference]["reference_length_nt"],
                f"{run}:{reference}:metric_reference_length", 1,
            )
            if length != metric_length:
                raise AggregateError(
                    f"metric_continuity_length_mismatch:{run}:{reference}:"
                    f"{metric_length}!={length}"
                )
            if reference in reference_lengths and reference_lengths[reference] != length:
                raise AggregateError(f"continuity_length_changed:{run}:{reference}:{length}")
            reference_lengths.setdefault(reference, length)
            boundary = integer(row.get("boundary_after_nt"), f"{run}:{reference}:boundary", 1)
            if boundary >= length or boundary % 250:
                raise AggregateError(f"continuity_boundary:{run}:{reference}:{boundary}:{length}")
            key = (reference, boundary)
            if key in seen:
                raise AggregateError(f"duplicate_continuity_boundary:{run}:{reference}:{boundary}")
            seen.add(key)
            if integer(
                row.get("direct_read_anchor_nt_each_side"),
                f"{run}:{reference}:{boundary}:direct_anchor", 1,
            ) != 25:
                raise AggregateError(
                    f"continuity_anchor_mismatch:{run}:{reference}:{boundary}"
                )
            spanning_field = f"proper_fragments_mapq{mapq}_nonduplicate_spanning"
            spanning = integer(row.get(spanning_field), f"{run}:{reference}:{boundary}:spanning")
            direct_fragments = integer(
                row.get("direct_read_spanning_fragments_anchor25"),
                f"{run}:{reference}:{boundary}:direct_fragments",
            )
            direct_reads = integer(
                row.get("direct_read_spanning_reads_anchor25"),
                f"{run}:{reference}:{boundary}:direct_reads",
            )
            non_field = f"proper_fragments_mapq{mapq}_nonduplicate"
            nonduplicate = integer(
                by_run[run][reference][non_field],
                f"{run}:{reference}:continuity_nonduplicate",
            )
            if (
                spanning > nonduplicate
                or direct_fragments > spanning
                or direct_reads < direct_fragments
                or direct_reads > 2 * direct_fragments
            ):
                raise AggregateError(f"continuity_subcount_mismatch:{run}:{reference}:{boundary}")
            all_continuity.append(dict(row))
        for reference in REFERENCES:
            lengths = {
                integer(
                    row["reference_length"],
                    f"{run}:{reference}:continuity_reference_length", 1,
                )
                for row in crows if row["reference"] == reference
            }
            if len(lengths) != 1:
                raise AggregateError(f"continuity_reference_length:{run}:{reference}")
            length = next(iter(lengths))
            expected = {(reference, boundary) for boundary in range(250, length, 250)}
            if {key for key in seen if key[0] == reference} != expected:
                raise AggregateError(f"continuity_boundary_set:{run}:{reference}")
    if set(by_run) != set(EXPECTED_RUNS):
        raise AggregateError(f"fragment_run_set:{sorted(by_run)}")
    for reference in REFERENCES:
        if reference_lengths.get(reference) != pinned_lengths[reference]:
            raise AggregateError(
                f"continuity_pinned_length_mismatch:{reference}:"
                f"{reference_lengths.get(reference)}!={pinned_lengths[reference]}"
            )
        if len(panel_kmers_available[reference]) != 1:
            raise AggregateError(
                f"panel_unique_kmer_available_changed:{reference}:"
                f"{sorted(panel_kmers_available[reference])}"
            )
    assert metric_fields and hash_fields and continuity_fields
    return by_run, all_metrics, all_hashes, all_continuity, metric_fields, hash_fields, continuity_fields


def build_sharing(
    hashes: list[dict[str, object]],
    by_run: dict[str, dict[str, dict[str, object]]],
    thresholds: dict[str, object],
) -> list[dict[str, object]]:
    by: dict[tuple[str, str], dict[str, set[str]]] = {}
    for run in EXPECTED_RUNS:
        for reference in REFERENCES:
            by[(run, reference)] = {"combined": set(), "sequence": set(), "endpoint": set()}
    for row in hashes:
        bucket = by[(str(row["run"]), str(row["reference"]))]
        bucket["combined"].add(str(row["fragment_fingerprint_sha256"]))
        bucket["sequence"].add(str(row["pair_sequence_sha256"]))
        bucket["endpoint"].add(str(row["endpoint_sha256"]))
    rows: list[dict[str, object]] = []
    for reference in REFERENCES:
        for run_a, run_b in itertools.combinations(EXPECTED_RUNS, 2):
            a, b = by[(run_a, reference)], by[(run_b, reference)]
            n_a, n_b = len(a["combined"]), len(b["combined"])
            non_field = f"proper_fragments_mapq{int(thresholds['minimum_mapq'])}_nonduplicate"
            fragments_a = integer(
                by_run[run_a][reference][non_field],
                f"{run_a}:{reference}:sharing_fragments",
            )
            fragments_b = integer(
                by_run[run_b][reference][non_field],
                f"{run_b}:{reference}:sharing_fragments",
            )
            if (fragments_a, run_a) <= (fragments_b, run_b):
                smaller_run, larger_run, smaller, larger = run_a, run_b, a, b
                smaller_fragments, larger_fragments = fragments_a, fragments_b
            else:
                smaller_run, larger_run, smaller, larger = run_b, run_a, b, a
                smaller_fragments, larger_fragments = fragments_b, fragments_a
            small_n, large_n = len(smaller["combined"]), len(larger["combined"])
            shared = len(smaller["combined"] & larger["combined"])
            containment = shared / small_n if small_n else 0.0
            union = len(a["combined"] | b["combined"])
            jaccard = shared / union if union else 0.0
            ratio = (
                larger_fragments / smaller_fragments
                if smaller_fragments else math.inf
            )
            fingerprint_eligible = small_n >= thresholds["shadow_min_unique"]
            if (
                fingerprint_eligible
                and ratio >= thresholds["shadow_abundance_ratio"]
                and shared >= thresholds["shadow_suspect_shared"]
                and containment >= thresholds["shadow_suspect_containment"]
            ):
                status = "shadow_suspect"
            elif (
                fingerprint_eligible
                and shared >= thresholds["shadow_warning_shared"]
                and containment >= thresholds["shadow_warning_containment"]
            ):
                status = "warning_shadow"
            elif small_n < thresholds["shadow_min_unique"]:
                status = "insufficient_smaller_run_fingerprints"
            else:
                status = "no_predeclared_shadow_signal"
            rows.append({
                "reference": reference, "run_a": run_a, "run_b": run_b,
                "unique_combined_a": n_a, "unique_combined_b": n_b,
                "smaller_run": smaller_run, "larger_run": larger_run,
                "smaller_unique_combined": small_n, "larger_unique_combined": large_n,
                "nonduplicate_fragments_a": fragments_a,
                "nonduplicate_fragments_b": fragments_b,
                "smaller_nonduplicate_fragments": smaller_fragments,
                "larger_nonduplicate_fragments": larger_fragments,
                "shared_combined_fingerprints": shared,
                "shared_sequence_hashes": len(a["sequence"] & b["sequence"]),
                "shared_endpoint_hashes": len(a["endpoint"] & b["endpoint"]),
                "smaller_run_containment": f"{containment:.6f}",
                "jaccard": f"{jaccard:.6f}",
                "larger_to_smaller_abundance_ratio": "" if not math.isfinite(ratio) else f"{ratio:.6f}",
                "shadow_qc_status": status,
            })
    return rows


def build_status(
    by_run: dict[str, dict[str, dict[str, object]]], sharing: list[dict[str, object]],
    thresholds: dict[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    suspect = {
        (str(row["smaller_run"]), str(row["reference"]))
        for row in sharing if row["shadow_qc_status"] == "shadow_suspect"
    }
    warning = {
        (str(row["smaller_run"]), str(row["reference"]))
        for row in sharing if row["shadow_qc_status"] == "warning_shadow"
    }
    insufficient = {
        (str(row["smaller_run"]), str(row["reference"]))
        for row in sharing
        if row["shadow_qc_status"] == "insufficient_smaller_run_fingerprints"
    }
    run_rows: list[dict[str, object]] = []
    for run in EXPECTED_RUNS:
        for reference in REFERENCES:
            row = by_run[run][reference]
            validated = row["fragment_validated_positive"] == "true"
            strong = row["fragment_strong_positive"] == "true"
            shadow = (
                "shadow_suspect" if (run, reference) in suspect else
                "warning_shadow" if (run, reference) in warning else
                "insufficient_for_shadow_assessment" if (run, reference) in insufficient else
                "no_predeclared_shadow_signal"
            )
            shadow_evaluable = (run, reference) not in insufficient
            independent_eligible = (
                (run, reference) not in suspect and shadow_evaluable
            )
            run_rows.append({
                "run": run, "alias": ALIASES[run], "reference": reference,
                "fragment_preduplicate_validated_positive": row["fragment_preduplicate_validated_positive"],
                "fragment_nonduplicate_validated_positive": row["fragment_nonduplicate_validated_positive"],
                "fragment_validated_positive": str(validated).lower(),
                "fragment_strong_positive": str(strong).lower(),
                "shadow_qc_status": shadow,
                "shadow_assessment_evaluable": str(shadow_evaluable).lower(),
                "eligible_as_independent_validated_support": str(
                    validated and independent_eligible
                ).lower(),
                "eligible_as_independent_strong_support": str(
                    strong and independent_eligible
                ).lower(),
            })
    candidate_rows: list[dict[str, object]] = []
    for reference in REFERENCES:
        relevant = [row for row in run_rows if row["reference"] == reference]
        validated_runs = [str(row["run"]) for row in relevant if row["fragment_validated_positive"] == "true"]
        strong_runs = [str(row["run"]) for row in relevant if row["fragment_strong_positive"] == "true"]
        suspect_runs = [str(row["run"]) for row in relevant if row["shadow_qc_status"] == "shadow_suspect"]
        warning_runs = [str(row["run"]) for row in relevant if row["shadow_qc_status"] == "warning_shadow"]
        insufficient_runs = [
            str(row["run"])
            for row in relevant
            if row["shadow_assessment_evaluable"] != "true"
        ]
        duplicate_sensitive = [
            str(row["run"]) for row in relevant
            if by_run[str(row["run"])][reference]["fragment_preduplicate_validated_positive"] == "true"
            and by_run[str(row["run"])][reference]["fragment_nonduplicate_validated_positive"] != "true"
        ]
        independent_validated = [str(row["run"]) for row in relevant if row["eligible_as_independent_validated_support"] == "true"]
        independent_strong = [str(row["run"]) for row in relevant if row["eligible_as_independent_strong_support"] == "true"]
        passed = bool(
            len(independent_validated) >= thresholds["candidate_validated_runs"]
            and len(independent_strong) >= thresholds["candidate_strong_runs"]
        )
        candidate_rows.append({
            "reference": reference,
            "validated_runs": ",".join(validated_runs), "strong_runs": ",".join(strong_runs),
            "shadow_suspect_runs": ",".join(suspect_runs),
            "warning_shadow_runs": ",".join(warning_runs),
            "shadow_assessment_insufficient_runs": ",".join(insufficient_runs),
            "duplicate_sensitive_runs": ",".join(duplicate_sensitive),
            "independent_validated_runs": ",".join(independent_validated),
            "independent_strong_runs": ",".join(independent_strong),
            "independent_validated_count": len(independent_validated),
            "independent_strong_count": len(independent_strong),
            "minimum_independent_validated_runs": thresholds["candidate_validated_runs"],
            "minimum_independent_strong_runs": thresholds["candidate_strong_runs"],
            "fragment_support_gate": "pass" if passed else "fail",
            "interpretation": "workflow_defined_technical_sequence_support_only",
        })
    return run_rows, candidate_rows


def write_checksums(out: Path) -> None:
    lines = [
        f"{sha256_path(path)}  {path.name}"
        for path in sorted(out.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (out / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")


def write_failure(
    out: Path, failure: str, threshold_sha: str = "", structure_sha: str = ""
) -> dict[str, object]:
    write_tsv(out / "FRAGMENT_METRICS.tsv", ["run", "reference"], [])
    write_tsv(out / "FRAGMENT_HASHES.tsv", ["run", "reference", "fragment_fingerprint_sha256"], [])
    write_tsv(out / "ALL_CONTINUITY_SCAN.tsv", ["run", "reference", "boundary_after_nt"], [])
    write_tsv(out / "RUN_REFERENCE_STATUS.tsv", ["run", "reference", "fragment_validated_positive"], [])
    write_tsv(out / "CROSS_RUN_FRAGMENT_SHARING.tsv", ["reference", "run_a", "run_b", "shadow_qc_status"], [])
    write_tsv(out / "CANDIDATE_FRAGMENT_GATE.tsv", ["reference", "fragment_support_gate"], [])
    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "technical_status": "technical_incomplete", "technical_complete": False,
        "audit_execution_status": "technical_incomplete",
        "overall_fragment_gate_status": "technical_incomplete",
        "failures": [failure], "fragment_support_gate": "technical_incomplete",
        "threshold_file_sha256": threshold_sha,
        "candidate_structure_sha256": structure_sha,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (out / "TECHNICAL_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    (out / "REPORT.md").write_text(
        "# Six-run proper-fragment audit\n\nTechnical status: **technical_incomplete**\n\n"
        f"- `{failure}`\n\n{CLAIM_BOUNDARY}\n"
    )
    write_checksums(out)
    return status


def run(
    root: Path, threshold_path: Path, out: Path,
    structure_path: Path | None = None,
) -> dict[str, object]:
    thresholds = load_thresholds(threshold_path)
    if structure_path is None:
        structure_path = threshold_path.with_name(
            str(thresholds["candidate_structure_filename"])
        )
    pinned_lengths = load_candidate_structure(structure_path, thresholds)
    by_run, metrics, hashes, continuity, metric_fields, hash_fields, continuity_fields = load_inputs(
        root, thresholds, pinned_lengths
    )
    sharing = build_sharing(hashes, by_run, thresholds)
    run_status, candidate = build_status(by_run, sharing, thresholds)
    write_tsv(out / "FRAGMENT_METRICS.tsv", metric_fields, metrics)
    write_tsv(out / "FRAGMENT_HASHES.tsv", hash_fields, hashes)
    write_tsv(out / "ALL_CONTINUITY_SCAN.tsv", continuity_fields, continuity)
    write_tsv(out / "RUN_REFERENCE_STATUS.tsv", list(run_status[0]), run_status)
    write_tsv(out / "CROSS_RUN_FRAGMENT_SHARING.tsv", list(sharing[0]), sharing)
    write_tsv(out / "CANDIDATE_FRAGMENT_GATE.tsv", list(candidate[0]), candidate)
    gate = "pass" if all(row["fragment_support_gate"] == "pass" for row in candidate) else "fail"
    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "technical_status": "pass", "technical_complete": True, "failures": [],
        "audit_execution_status": "complete",
        "overall_fragment_gate_status": gate,
        "fragment_support_gate": gate, "expected_runs": list(EXPECTED_RUNS),
        "expected_references": list(REFERENCES), "metric_row_count": len(metrics),
        "hash_row_count": len(hashes), "sharing_row_count": len(sharing),
        "candidate_structure_path": str(structure_path),
        "candidate_structure_sha256": sha256_path(structure_path),
        "workflow_defined_thresholds": thresholds, "claim_boundary": CLAIM_BOUNDARY,
        "shadow_interpretation": "QC suspicion only; not proof of index hopping or contamination.",
    }
    (out / "TECHNICAL_STATUS.json").write_text(json.dumps(status, indent=2) + "\n")
    report = [
        "# Six-run duplicate-aware proper-fragment audit", "",
        "Technical status: **pass**", "", f"Fragment support gate: **{gate}**", "",
        "| Reference | Independent validated runs | Independent strong runs | Gate |",
        "|---|---|---|---|",
    ]
    report.extend(
        f"| `{row['reference']}` | `{row['independent_validated_runs']}` | "
        f"`{row['independent_strong_runs']}` | `{row['fragment_support_gate']}` |"
        for row in candidate
    )
    report.extend([
        "", "Exact cross-run fragment sharing is a batch-QC indicator only; it cannot prove or exclude index hopping.",
        "", CLAIM_BOUNDARY,
    ])
    (out / "REPORT.md").write_text("\n".join(report) + "\n")
    write_checksums(out)
    return status


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--thresholds", required=True, type=Path)
    parser.add_argument(
        "--structure", type=Path,
        help=(
            "pinned candidate-structure TSV; defaults to the contract filename "
            "beside --thresholds"
        ),
    )
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")
    threshold_sha = sha256_path(args.thresholds) if args.thresholds.is_file() else ""
    default_structure = args.thresholds.with_name("panax_candidate_structure.tsv")
    failure_structure = args.structure or default_structure
    structure_sha = sha256_path(failure_structure) if failure_structure.is_file() else ""
    try:
        status = run(
            args.input_root, args.thresholds, args.out,
            structure_path=args.structure,
        )
    except (AggregateError, OSError, KeyError, AssertionError, ValueError, TypeError) as exc:
        status = write_failure(
            args.out, f"{type(exc).__name__}:{exc}", threshold_sha, structure_sha
        )
    print(json.dumps(status, indent=2))
    return 0 if status.get("technical_complete") and status.get("fragment_support_gate") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
