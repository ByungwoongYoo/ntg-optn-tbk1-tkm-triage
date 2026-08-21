#!/usr/bin/env python3
"""Fail-closed audit of mapping-seed-free Panax de novo assemblies.

The input BLAST table must use this exact outfmt 6 field order::

    qseqid sseqid pident length nident mismatch gapopen gaps qlen slen
    qstart qend sstart send qseq sseq evalue bitscore

Recovery is evaluated on one assembly contig at a time.  HSPs on the selected
contig may be combined only when their relative orientation and coordinate
order are coherent.  The source/non-source thresholds are read from the
workflow-defined threshold JSON; they are technical recovery rules rather
than biological detection limits or taxonomic criteria.

The provenance JSON is deliberately mandatory.  It must identify the run and
assembly hash and attest that the complete paired FASTQ was assembled without
candidate baiting, mapping seeds, reference guidance, or target-read
selection.  The attestation is auditable workflow provenance, not a way to
infer assembly independence from the FASTA alone.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


BLAST_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "nident", "mismatch",
    "gapopen", "gaps", "qlen", "slen", "qstart", "qend", "sstart",
    "send", "qseq", "sseq", "evalue", "bitscore",
)
RECOVERY_FIELDS = (
    "run", "query", "source_run", "run_role", "query_length_nt",
    "query_coverage_threshold", "identity_threshold",
    "maximum_internal_query_gap_nt", "maximum_internal_subject_gap_nt",
    "best_single_contig", "best_contig_length_nt",
    "best_single_contig_query_coverage", "coordinate_weighted_identity",
    "covered_query_nt", "max_internal_query_gap_nt",
    "max_internal_subject_gap_nt", "relative_orientation",
    "collinear_hsp_geometry", "hsp_count", "best_evalue", "total_bitscore",
    "gate_status", "recovery_status",
)
HIT_FIELDS = (
    "run", "query", "contig", "contig_length_nt", "contig_sha256",
    "single_contig_query_coverage", "coordinate_weighted_identity",
    "max_internal_query_gap_nt", "max_internal_subject_gap_nt",
    "relative_orientation",
    "collinear_hsp_geometry", "hsp_count", "selected_as_best",
)
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"
)
MIN_HIT_HSP_LENGTH = 100
MIN_HIT_PIDENT = 70.0
MAX_HIT_EVALUE = 1e-5
CLAIM_BOUNDARY = (
    "This audit tests mapping-seed-free, single-contig recovery of predefined "
    "partial RNA-sequence queries. It makes no formal taxonomic, biological "
    "host, active-replication, phenotype-association, causality, pathogenicity, "
    "transmission, or agricultural/medical-effect claim."
)


class AuditError(ValueError):
    """Input or provenance failure that makes the audit technically incomplete."""


def reverse_complement(sequence: str) -> str:
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class Hsp:
    query: str
    contig: str
    pident: float
    alignment_length: int
    nident: int
    mismatch: int
    gapopen: int
    gaps: int
    query_length: int
    contig_length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    aligned_query: str
    aligned_subject: str
    evalue: float
    bitscore: float

    @property
    def qlo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhi(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def subject_midpoint(self) -> float:
        return (self.sstart + self.send) / 2.0

    @property
    def relative_orientation(self) -> str:
        query_direction = 1 if self.qend >= self.qstart else -1
        subject_direction = 1 if self.send >= self.sstart else -1
        return "plus" if query_direction == subject_direction else "minus"

    @property
    def query_span(self) -> int:
        return self.qhi - self.qlo + 1

    @property
    def recovered_query_nt(self) -> int:
        return sum(
            query_base != "-" and subject_base != "-"
            for query_base, subject_base in zip(
                self.aligned_query, self.aligned_subject
            )
        )

    @staticmethod
    def longest_gap(sequence: str) -> int:
        longest = current = 0
        for base in sequence:
            if base == "-":
                current += 1
                longest = max(longest, current)
            else:
                current = 0
        return longest

    @property
    def internal_query_gap(self) -> int:
        return self.longest_gap(self.aligned_subject)

    @property
    def internal_subject_gap(self) -> int:
        return self.longest_gap(self.aligned_query)


@dataclass(frozen=True)
class ChainState:
    """A query-disjoint, subject-monotone chain of BLAST HSPs."""

    hsps: tuple[Hsp, ...]
    covered_query_nt: int
    identical_alignment_nt: int
    alignment_columns: int
    total_bitscore: float
    max_internal_query_gap: int
    max_internal_subject_gap: int

    @property
    def weighted_identity(self) -> float:
        if not self.alignment_columns:
            return 0.0
        return self.identical_alignment_nt / self.alignment_columns


@dataclass(frozen=True)
class ContigSummary:
    query: str
    contig: str
    query_length: int
    contig_length: int
    query_coverage: float
    weighted_identity: float
    covered_query_nt: int
    max_internal_gap: int
    max_internal_subject_gap: int
    orientation: str
    collinear: bool
    hsp_count: int
    best_evalue: float
    total_bitscore: float


@dataclass(frozen=True)
class AuditOutputs:
    recovery_rows: list[dict[str, object]]
    hit_rows: list[dict[str, object]]
    hit_contigs: list[FastaRecord]
    technical_status: str
    failures: list[str]
    provenance: dict[str, object]
    thresholds: dict[str, object]


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def read_fasta(path: Path, label: str) -> dict[str, FastaRecord]:
    if not path.is_file():
        raise AuditError(f"missing_{label}_fasta:{path}")
    records: dict[str, FastaRecord] = {}
    identifier: str | None = None
    description = ""
    chunks: list[str] = []

    def store() -> None:
        if identifier is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence:
            raise AuditError(f"empty_{label}_sequence:{identifier}")
        unexpected = sorted(set(sequence) - IUPAC_DNA)
        if unexpected:
            raise AuditError(
                f"invalid_{label}_sequence:{identifier}:{''.join(unexpected)}"
            )
        records[identifier] = FastaRecord(identifier, description, sequence)

    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            store()
            description = line[1:].strip()
            identifier = description.split()[0] if description else ""
            if not identifier:
                raise AuditError(f"empty_{label}_fasta_id:line_{number}")
            if identifier in records:
                raise AuditError(f"duplicate_{label}_fasta_id:{identifier}")
            chunks = []
        elif identifier is None:
            raise AuditError(f"sequence_before_{label}_fasta_header:line_{number}")
        else:
            chunks.append("".join(line.split()))
    store()
    if not records:
        raise AuditError(f"empty_{label}_fasta:{path}")
    return records


def require_number(
    value: object,
    label: str,
    minimum: float,
    maximum: float,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AuditError(f"invalid_threshold_type:{label}")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise AuditError(f"invalid_threshold_value:{label}:{value}")
    return number


def load_thresholds(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise AuditError(f"missing_thresholds:{path}")
    try:
        value = json.loads(path.read_text())
    except Exception as exc:  # pragma: no cover - exact decoder text is unstable
        raise AuditError(f"invalid_thresholds_json:{type(exc).__name__}") from exc
    if not isinstance(value, dict):
        raise AuditError("thresholds_must_be_object")
    section = value.get("de_novo_recovery")
    if not isinstance(section, dict):
        raise AuditError("missing_de_novo_recovery_thresholds")
    maximum_internal_gap = require_number(
        section.get("maximum_internal_query_gap_nt"),
        "maximum_internal_query_gap_nt", 0.0, 1_000_000.0,
    )
    if not maximum_internal_gap.is_integer():
        raise AuditError(
            "invalid_threshold_value:maximum_internal_query_gap_nt:"
            f"{maximum_internal_gap}"
        )
    maximum_internal_subject_gap = require_number(
        section.get("maximum_internal_subject_gap_nt"),
        "maximum_internal_subject_gap_nt", 0.0, 1_000_000.0,
    )
    if not maximum_internal_subject_gap.is_integer():
        raise AuditError(
            "invalid_threshold_value:maximum_internal_subject_gap_nt:"
            f"{maximum_internal_subject_gap}"
        )
    normalized = {
        "source_query_coverage": require_number(
            section.get("source_run_single_contig_query_coverage_minimum"),
            "source_query_coverage", 0.0, 1.0,
        ),
        "source_identity": require_number(
            section.get("source_run_identity_minimum"),
            "source_identity", 0.0, 1.0,
        ),
        "non_source_query_coverage": require_number(
            section.get("non_source_single_contig_query_coverage_minimum"),
            "non_source_query_coverage", 0.0, 1.0,
        ),
        "non_source_identity": require_number(
            section.get("non_source_identity_minimum"),
            "non_source_identity", 0.0, 1.0,
        ),
        "maximum_internal_query_gap_nt": int(maximum_internal_gap),
        "maximum_internal_subject_gap_nt": int(maximum_internal_subject_gap),
        "candidate_support_rule": section.get("candidate_support_rule", ""),
        "scope": value.get("scope", ""),
        "source_file_sha256": sha256_path(path),
    }
    if not isinstance(normalized["candidate_support_rule"], str):
        raise AuditError("invalid_candidate_support_rule")
    return normalized


def load_structure(
    path: Path,
    queries: dict[str, FastaRecord],
) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise AuditError(f"missing_structure_manifest:{path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference", "source_run", "reference_length_nt"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise AuditError(
                f"structure_manifest_fields:{reader.fieldnames}:required={sorted(required)}"
            )
        rows = list(reader)
    manifest: dict[str, dict[str, str]] = {}
    for number, row in enumerate(rows, 2):
        query = row.get("reference", "")
        if not query:
            raise AuditError(f"empty_structure_reference:line_{number}")
        if query in manifest:
            raise AuditError(f"duplicate_structure_reference:{query}")
        if not row.get("source_run"):
            raise AuditError(f"empty_structure_source_run:{query}")
        try:
            declared_length = int(row.get("reference_length_nt", ""))
        except ValueError as exc:
            raise AuditError(f"invalid_structure_reference_length:{query}") from exc
        if query not in queries:
            raise AuditError(f"unexpected_structure_reference:{query}")
        if declared_length != len(queries[query].sequence):
            raise AuditError(
                f"structure_query_length_mismatch:{query}:"
                f"{declared_length}!={len(queries[query].sequence)}"
            )
        manifest[query] = row
    if set(manifest) != set(queries):
        raise AuditError(
            "structure_query_set_mismatch:"
            f"missing={sorted(set(queries)-set(manifest))}:"
            f"unexpected={sorted(set(manifest)-set(queries))}"
        )
    return manifest


def load_provenance(
    path: Path, run: str, assembly: Path,
    query_path: Path, structure_path: Path,
) -> dict[str, object]:
    if not path.is_file():
        raise AuditError(f"missing_assembly_provenance:{path}")
    try:
        value = json.loads(path.read_text())
    except Exception as exc:  # pragma: no cover
        raise AuditError(f"invalid_assembly_provenance_json:{type(exc).__name__}") from exc
    if not isinstance(value, dict):
        raise AuditError("assembly_provenance_must_be_object")
    required = {
        "run", "assembly_sha256", "input_scope", "candidate_baiting",
        "mapping_seeded", "reference_guided", "target_read_selection",
        "assembler", "assembler_version", "assembly_file",
        "assembly_exit_code", "search_exit_code",
        "fastq_sha256_manifest", "fastq_sha256_manifest_sha256",
        "assembly_method_manifest", "assembly_method_manifest_sha256",
        "retained_assembly_file", "retained_assembly_manifest",
        "retained_assembly_compression", "candidate_query_sha256",
        "candidate_structure_sha256",
    }
    if set(value) < required:
        raise AuditError(
            f"assembly_provenance_missing_fields:{sorted(required-set(value))}"
        )
    if value["run"] != run:
        raise AuditError(f"assembly_provenance_run_mismatch:{value['run']}!={run}")
    observed_sha = sha256_path(assembly)
    if value["assembly_sha256"] != observed_sha:
        raise AuditError(
            "assembly_provenance_sha256_mismatch:"
            f"{value['assembly_sha256']}!={observed_sha}"
        )
    assembly_from_manifest = (path.parent / str(value["assembly_file"])).resolve()
    if assembly_from_manifest != assembly.resolve():
        raise AuditError(
            f"assembly_provenance_path_mismatch:{assembly_from_manifest}!={assembly.resolve()}"
        )
    if value["input_scope"] != "complete_paired_fastq":
        raise AuditError(f"assembly_input_scope_not_complete:{value['input_scope']}")
    for field in (
        "candidate_baiting", "mapping_seeded", "reference_guided",
        "target_read_selection",
    ):
        if value[field] is not False:
            raise AuditError(f"assembly_not_mapping_seed_free:{field}={value[field]!r}")
    for field in ("assembler", "assembler_version"):
        if not isinstance(value[field], str) or not value[field].strip():
            raise AuditError(f"assembly_provenance_empty_field:{field}")
    if value["assembler"] != "MEGAHIT":
        raise AuditError(f"unexpected_assembler:{value['assembler']}")
    for field in ("assembly_exit_code", "search_exit_code"):
        if str(value[field]) != "0":
            raise AuditError(f"nonzero_complete_provenance_exit:{field}={value[field]}")
    for field, source in (
        ("candidate_query_sha256", query_path),
        ("candidate_structure_sha256", structure_path),
    ):
        if value[field] != sha256_path(source):
            raise AuditError(f"provenance_input_sha256_mismatch:{field}")
    if value["retained_assembly_compression"] != "gzip_nondeterministic_metadata_disabled":
        raise AuditError("unexpected_retained_assembly_compression")
    for field in ("retained_assembly_file", "retained_assembly_manifest"):
        if not isinstance(value[field], str) or not value[field].strip():
            raise AuditError(f"assembly_provenance_empty_field:{field}")
    evidence_root = path.parent.parent.resolve()
    for path_field, sha_field in (
        ("fastq_sha256_manifest", "fastq_sha256_manifest_sha256"),
        ("assembly_method_manifest", "assembly_method_manifest_sha256"),
    ):
        evidence_path = (path.parent / str(value[path_field])).resolve()
        try:
            evidence_path.relative_to(evidence_root)
        except ValueError as exc:
            raise AuditError(f"provenance_manifest_escapes_evidence_root:{path_field}") from exc
        if not evidence_path.is_file():
            raise AuditError(f"missing_provenance_manifest:{path_field}:{evidence_path}")
        observed_manifest_sha = sha256_path(evidence_path)
        if value[sha_field] != observed_manifest_sha:
            raise AuditError(
                f"provenance_manifest_sha256_mismatch:{path_field}:"
                f"{value[sha_field]}!={observed_manifest_sha}"
            )
    return value


def parse_int(text: str, label: str, row_number: int, minimum: int = 1) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise AuditError(f"invalid_blast_{label}:row_{row_number}:{text!r}") from exc
    if value < minimum:
        raise AuditError(f"invalid_blast_{label}:row_{row_number}:{value}")
    return value


def parse_float(
    text: str,
    label: str,
    row_number: int,
    minimum: float,
    maximum: float | None = None,
) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise AuditError(f"invalid_blast_{label}:row_{row_number}:{text!r}") from exc
    if not math.isfinite(value) or value < minimum or (
        maximum is not None and value > maximum
    ):
        raise AuditError(f"invalid_blast_{label}:row_{row_number}:{value}")
    return value


def parse_blast(
    path: Path,
    queries: dict[str, FastaRecord],
    assembly: dict[str, FastaRecord],
) -> list[Hsp]:
    if not path.is_file():
        raise AuditError(f"missing_blast_table:{path}")
    hsps: list[Hsp] = []
    seen: set[tuple[str, ...]] = set()
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = tuple(line.split("\t"))
        if len(fields) != len(BLAST_FIELDS):
            raise AuditError(
                f"malformed_blast_row:{number}:"
                f"expected_{len(BLAST_FIELDS)}_fields_observed_{len(fields)}"
            )
        if fields in seen:
            raise AuditError(f"duplicate_blast_row:{number}")
        seen.add(fields)
        row = dict(zip(BLAST_FIELDS, fields))
        query = row["qseqid"]
        contig = row["sseqid"]
        if query not in queries:
            raise AuditError(f"unexpected_blast_query:row_{number}:{query}")
        if contig not in assembly:
            raise AuditError(f"unexpected_blast_contig:row_{number}:{contig}")
        hsp = Hsp(
            query=query,
            contig=contig,
            pident=parse_float(row["pident"], "pident", number, 0.0, 100.0),
            alignment_length=parse_int(row["length"], "length", number),
            nident=parse_int(row["nident"], "nident", number, 0),
            mismatch=parse_int(row["mismatch"], "mismatch", number, 0),
            gapopen=parse_int(row["gapopen"], "gapopen", number, 0),
            gaps=parse_int(row["gaps"], "gaps", number, 0),
            query_length=parse_int(row["qlen"], "qlen", number),
            contig_length=parse_int(row["slen"], "slen", number),
            qstart=parse_int(row["qstart"], "qstart", number),
            qend=parse_int(row["qend"], "qend", number),
            sstart=parse_int(row["sstart"], "sstart", number),
            send=parse_int(row["send"], "send", number),
            aligned_query=row["qseq"].upper(),
            aligned_subject=row["sseq"].upper(),
            evalue=parse_float(row["evalue"], "evalue", number, 0.0),
            bitscore=parse_float(row["bitscore"], "bitscore", number, 0.0),
        )
        actual_query_length = len(queries[query].sequence)
        actual_contig_length = len(assembly[contig].sequence)
        if hsp.query_length != actual_query_length:
            raise AuditError(
                f"blast_query_length_mismatch:row_{number}:{query}:"
                f"{hsp.query_length}!={actual_query_length}"
            )
        if hsp.contig_length != actual_contig_length:
            raise AuditError(
                f"blast_contig_length_mismatch:row_{number}:{contig}:"
                f"{hsp.contig_length}!={actual_contig_length}"
            )
        if hsp.qhi > hsp.query_length or hsp.shi > hsp.contig_length:
            raise AuditError(f"blast_coordinate_out_of_bounds:row_{number}")
        if (
            len(hsp.aligned_query) != hsp.alignment_length
            or len(hsp.aligned_subject) != hsp.alignment_length
        ):
            raise AuditError(f"blast_aligned_sequence_length_mismatch:row_{number}")
        if any(
            base not in IUPAC_DNA and base != "-"
            for base in hsp.aligned_query + hsp.aligned_subject
        ):
            raise AuditError(f"invalid_blast_aligned_sequence:row_{number}")
        if any(
            query_base == "-" and subject_base == "-"
            for query_base, subject_base in zip(
                hsp.aligned_query, hsp.aligned_subject
            )
        ):
            raise AuditError(f"double_gap_blast_column:row_{number}")
        if hsp.aligned_query.replace("-", "") == "" or hsp.aligned_subject.replace("-", "") == "":
            raise AuditError(f"empty_blast_aligned_residue_set:row_{number}")
        if len(hsp.aligned_query.replace("-", "")) != hsp.query_span:
            raise AuditError(f"blast_query_span_sequence_mismatch:row_{number}")
        if len(hsp.aligned_subject.replace("-", "")) != hsp.shi - hsp.slo + 1:
            raise AuditError(f"blast_subject_span_sequence_mismatch:row_{number}")
        expected_query = queries[query].sequence[hsp.qlo - 1:hsp.qhi]
        if hsp.qend < hsp.qstart:
            expected_query = reverse_complement(expected_query)
        expected_subject = assembly[contig].sequence[hsp.slo - 1:hsp.shi]
        if hsp.send < hsp.sstart:
            expected_subject = reverse_complement(expected_subject)
        if hsp.aligned_query.replace("-", "") != expected_query:
            raise AuditError(f"blast_query_sequence_mismatch:row_{number}")
        if hsp.aligned_subject.replace("-", "") != expected_subject:
            raise AuditError(f"blast_subject_sequence_mismatch:row_{number}")
        observed_nident = sum(
            query_base == subject_base and query_base != "-"
            for query_base, subject_base in zip(
                hsp.aligned_query, hsp.aligned_subject
            )
        )
        observed_mismatch = sum(
            query_base != subject_base
            and query_base != "-" and subject_base != "-"
            for query_base, subject_base in zip(
                hsp.aligned_query, hsp.aligned_subject
            )
        )
        observed_gaps = sum(
            query_base == "-" or subject_base == "-"
            for query_base, subject_base in zip(
                hsp.aligned_query, hsp.aligned_subject
            )
        )
        observed_gapopen = 0
        previous_gap_side = ""
        for query_base, subject_base in zip(
            hsp.aligned_query, hsp.aligned_subject
        ):
            gap_side = "query" if query_base == "-" else (
                "subject" if subject_base == "-" else ""
            )
            if gap_side and gap_side != previous_gap_side:
                observed_gapopen += 1
            previous_gap_side = gap_side
        if (
            (hsp.nident, hsp.mismatch, hsp.gaps, hsp.gapopen)
            != (observed_nident, observed_mismatch, observed_gaps, observed_gapopen)
        ):
            raise AuditError(f"blast_alignment_count_mismatch:row_{number}")
        if hsp.nident + hsp.mismatch + hsp.gaps != hsp.alignment_length:
            raise AuditError(f"blast_alignment_partition_mismatch:row_{number}")
        expected_pident = 100.0 * hsp.nident / hsp.alignment_length
        if not math.isclose(hsp.pident, expected_pident, abs_tol=0.001):
            raise AuditError(f"blast_pident_mismatch:row_{number}")
        hsps.append(hsp)
    return hsps


def chain_link_gaps(left: Hsp, right: Hsp, orientation: str) -> tuple[int, int] | None:
    """Return query/subject gaps for a strict non-overlapping monotone link."""
    if (
        orientation not in {"plus", "minus"}
        or left.relative_orientation != orientation
        or right.relative_orientation != orientation
        or right.qlo <= left.qhi
    ):
        return None
    query_gap = right.qlo - left.qhi - 1
    if orientation == "plus":
        if right.slo <= left.shi:
            return None
        subject_gap = right.slo - left.shi - 1
    else:
        if right.shi >= left.slo:
            return None
        subject_gap = left.slo - right.shi - 1
    return query_gap, subject_gap


def extend_chain(state: ChainState, hsp: Hsp, gaps: tuple[int, int]) -> ChainState:
    query_gap, subject_gap = gaps
    return ChainState(
        hsps=state.hsps + (hsp,),
        covered_query_nt=state.covered_query_nt + hsp.recovered_query_nt,
        identical_alignment_nt=state.identical_alignment_nt + hsp.nident,
        alignment_columns=state.alignment_columns + hsp.alignment_length,
        total_bitscore=state.total_bitscore + hsp.bitscore,
        max_internal_query_gap=max(
            state.max_internal_query_gap, query_gap, hsp.internal_query_gap
        ),
        max_internal_subject_gap=max(
            state.max_internal_subject_gap, subject_gap,
            hsp.internal_subject_gap,
        ),
    )


def prune_chain_states(
    states: Sequence[ChainState], identity_threshold: float
) -> list[ChainState]:
    """Keep states not dominated for coverage and identity-threshold surplus."""
    def surplus(state: ChainState) -> float:
        return (
            state.identical_alignment_nt
            - identity_threshold * state.alignment_columns
        )

    ordered = sorted(
        states,
        key=lambda state: (
            -state.covered_query_nt,
            -surplus(state),
            -state.weighted_identity,
            -state.total_bitscore,
            len(state.hsps),
        ),
    )
    kept: list[ChainState] = []
    maximum_surplus_seen = -math.inf
    for candidate in ordered:
        candidate_surplus = surplus(candidate)
        # ``ordered`` already guarantees every seen state has at least this
        # much query coverage.  A running maximum therefore performs the same
        # two-dimensional Pareto test in O(k log k), avoiding the quadratic
        # frontier scan that can stall on repetitive assemblies.
        if maximum_surplus_seen >= candidate_surplus - 1e-9:
            continue
        kept.append(candidate)
        maximum_surplus_seen = max(maximum_surplus_seen, candidate_surplus)
    return kept


def build_chain_states(
    hsps: Sequence[Hsp],
    orientation: str,
    identity_threshold: float,
    maximum_internal_query_gap: int | None,
    maximum_internal_subject_gap: int | None,
) -> list[ChainState]:
    ordered = sorted(
        (hsp for hsp in hsps if hsp.relative_orientation == orientation),
        key=lambda hsp: (
            hsp.qlo, hsp.qhi, hsp.slo, hsp.shi, -hsp.bitscore
        ),
    )
    ending: list[list[ChainState]] = []
    all_states: list[ChainState] = []
    for index, hsp in enumerate(ordered):
        candidates = [ChainState(
            hsps=(hsp,),
            covered_query_nt=hsp.recovered_query_nt,
            identical_alignment_nt=hsp.nident,
            alignment_columns=hsp.alignment_length,
            total_bitscore=hsp.bitscore,
            max_internal_query_gap=hsp.internal_query_gap,
            max_internal_subject_gap=hsp.internal_subject_gap,
        )]
        for previous_index, previous in enumerate(ordered[:index]):
            gaps = chain_link_gaps(previous, hsp, orientation)
            if gaps is None:
                continue
            query_gap, subject_gap = gaps
            if (
                maximum_internal_query_gap is not None
                and query_gap > maximum_internal_query_gap
            ) or (
                maximum_internal_subject_gap is not None
                and subject_gap > maximum_internal_subject_gap
            ):
                continue
            candidates.extend(
                extend_chain(state, hsp, gaps)
                for state in ending[previous_index]
            )
        ending.append(prune_chain_states(candidates, identity_threshold))
        all_states.extend(ending[-1])
    return all_states


def chain_summary(
    state: ChainState, query: str, contig: str,
    query_length: int, contig_length: int,
) -> ContigSummary:
    orientations = {hsp.relative_orientation for hsp in state.hsps}
    orientation = next(iter(orientations)) if len(orientations) == 1 else "mixed"
    return ContigSummary(
        query=query,
        contig=contig,
        query_length=query_length,
        contig_length=contig_length,
        query_coverage=state.covered_query_nt / query_length,
        weighted_identity=state.weighted_identity,
        covered_query_nt=state.covered_query_nt,
        max_internal_gap=state.max_internal_query_gap,
        max_internal_subject_gap=state.max_internal_subject_gap,
        orientation=orientation,
        collinear=orientation in {"plus", "minus"},
        hsp_count=len(state.hsps),
        best_evalue=min(hsp.evalue for hsp in state.hsps),
        total_bitscore=state.total_bitscore,
    )


def summarize_contig(
    hsps: Sequence[Hsp],
    coverage_threshold: float,
    identity_threshold: float,
    maximum_internal_query_gap: int,
    maximum_internal_subject_gap: int,
) -> ContigSummary:
    if not hsps:
        raise AuditError("cannot_summarize_empty_hsp_group")
    query = hsps[0].query
    contig = hsps[0].contig
    query_length = hsps[0].query_length
    contig_length = hsps[0].contig_length
    if any(
        hsp.query != query or hsp.contig != contig
        or hsp.query_length != query_length or hsp.contig_length != contig_length
        for hsp in hsps
    ):
        raise AuditError(f"mixed_hsp_group:{query}:{contig}")

    constrained: list[ChainState] = []
    unconstrained: list[ChainState] = []
    for orientation in ("plus", "minus"):
        constrained.extend(build_chain_states(
            hsps, orientation, identity_threshold,
            maximum_internal_query_gap, maximum_internal_subject_gap,
        ))
        unconstrained.extend(build_chain_states(
            hsps, orientation, identity_threshold, None, None,
        ))
    eligible = [
        state for state in constrained
        if state.covered_query_nt / query_length >= coverage_threshold
        and state.weighted_identity >= identity_threshold
    ]
    pool = eligible if eligible else unconstrained
    if not pool:  # All parsed HSPs necessarily have plus/minus orientation.
        raise AuditError(f"no_chain_state:{query}:{contig}")
    best = sorted(
        pool,
        key=lambda state: (
            -state.covered_query_nt,
            -state.weighted_identity,
            state.max_internal_query_gap,
            state.max_internal_subject_gap,
            -state.total_bitscore,
            len(state.hsps),
        ),
    )[0]
    summary = chain_summary(best, query, contig, query_length, contig_length)
    if not eligible and len({hsp.relative_orientation for hsp in hsps}) > 1:
        return ContigSummary(
            **{
                **summary.__dict__,
                "orientation": "mixed",
                "collinear": False,
            }
        )
    return summary


def choose_best(
    summaries: Sequence[ContigSummary],
    coverage_threshold: float,
    identity_threshold: float,
    maximum_internal_query_gap: int,
    maximum_internal_subject_gap: int,
) -> ContigSummary | None:
    if not summaries:
        return None
    return sorted(
        summaries,
        key=lambda row: (
            not (
                row.collinear
                and row.query_coverage >= coverage_threshold
                and row.weighted_identity >= identity_threshold
                and row.max_internal_gap <= maximum_internal_query_gap
                and row.max_internal_subject_gap <= maximum_internal_subject_gap
            ),
            not row.collinear,
            -row.query_coverage,
            -row.weighted_identity,
            -row.total_bitscore,
            row.contig,
        ),
    )[0]


def recovery_thresholds(
    run: str,
    source_run: str,
    thresholds: dict[str, object],
) -> tuple[str, float, float]:
    if run == source_run:
        return (
            "source",
            float(thresholds["source_query_coverage"]),
            float(thresholds["source_identity"]),
        )
    return (
        "non_source",
        float(thresholds["non_source_query_coverage"]),
        float(thresholds["non_source_identity"]),
    )


def run_audit(
    run: str,
    query_path: Path,
    assembly_path: Path,
    blast_path: Path,
    structure_path: Path,
    threshold_path: Path,
    provenance_path: Path,
) -> AuditOutputs:
    if not run or any(char.isspace() for char in run):
        raise AuditError(f"invalid_run:{run!r}")
    queries = read_fasta(query_path, "query")
    assembly = read_fasta(assembly_path, "assembly")
    thresholds = load_thresholds(threshold_path)
    structure = load_structure(structure_path, queries)
    provenance = load_provenance(
        provenance_path, run, assembly_path, query_path, structure_path
    )
    hsps = parse_blast(blast_path, queries, assembly)

    grouped: dict[tuple[str, str], list[Hsp]] = {}
    for hsp in hsps:
        grouped.setdefault((hsp.query, hsp.contig), []).append(hsp)
    maximum_internal_query_gap = int(
        thresholds["maximum_internal_query_gap_nt"]
    )
    maximum_internal_subject_gap = int(
        thresholds["maximum_internal_subject_gap_nt"]
    )
    summaries: dict[tuple[str, str], ContigSummary] = {}
    for (query, contig), group in sorted(grouped.items()):
        _, coverage_threshold, identity_threshold = recovery_thresholds(
            run, structure[query]["source_run"], thresholds
        )
        summaries[(query, contig)] = summarize_contig(
            group,
            coverage_threshold,
            identity_threshold,
            maximum_internal_query_gap,
            maximum_internal_subject_gap,
        )

    recovery_rows: list[dict[str, object]] = []
    selected: dict[str, str] = {}
    for query, query_record in queries.items():
        query_summaries = [
            summary for (group_query, _), summary in summaries.items()
            if group_query == query
        ]
        source_run = structure[query]["source_run"]
        role, coverage_threshold, identity_threshold = recovery_thresholds(
            run, source_run, thresholds
        )
        best = choose_best(
            query_summaries,
            coverage_threshold,
            identity_threshold,
            maximum_internal_query_gap,
            maximum_internal_subject_gap,
        )
        passed = bool(
            best is not None
            and best.collinear
            and best.query_coverage >= coverage_threshold
            and best.weighted_identity >= identity_threshold
            and best.max_internal_gap
            <= maximum_internal_query_gap
            and best.max_internal_subject_gap
            <= maximum_internal_subject_gap
        )
        if best is not None:
            selected[query] = best.contig
        recovery_rows.append({
            "run": run,
            "query": query,
            "source_run": source_run,
            "run_role": role,
            "query_length_nt": len(query_record.sequence),
            "query_coverage_threshold": f"{coverage_threshold:.6f}",
            "identity_threshold": f"{identity_threshold:.6f}",
            "maximum_internal_query_gap_nt": thresholds[
                "maximum_internal_query_gap_nt"
            ],
            "maximum_internal_subject_gap_nt": thresholds[
                "maximum_internal_subject_gap_nt"
            ],
            "best_single_contig": best.contig if best else "",
            "best_contig_length_nt": best.contig_length if best else "",
            "best_single_contig_query_coverage": (
                f"{best.query_coverage:.6f}" if best else "0.000000"
            ),
            "coordinate_weighted_identity": (
                f"{best.weighted_identity:.6f}" if best else ""
            ),
            "covered_query_nt": best.covered_query_nt if best else 0,
            "max_internal_query_gap_nt": best.max_internal_gap if best else "",
            "max_internal_subject_gap_nt": (
                best.max_internal_subject_gap if best else ""
            ),
            "relative_orientation": best.orientation if best else "",
            "collinear_hsp_geometry": str(best.collinear).lower() if best else "false",
            "hsp_count": best.hsp_count if best else 0,
            "best_evalue": f"{best.best_evalue:.6g}" if best else "",
            "total_bitscore": f"{best.total_bitscore:.6f}" if best else "",
            "gate_status": "pass" if passed else "fail",
            "recovery_status": (
                f"{role}_recovered" if passed else f"{role}_not_recovered"
            ),
        })

    hit_rows: list[dict[str, object]] = []
    hit_contig_ids: set[str] = set()
    for (query, contig), summary in summaries.items():
        group = grouped[(query, contig)]
        eligible = any(
            hsp.alignment_length >= MIN_HIT_HSP_LENGTH
            and hsp.pident >= MIN_HIT_PIDENT
            and hsp.evalue <= MAX_HIT_EVALUE
            for hsp in group
        )
        if not eligible:
            continue
        hit_contig_ids.add(contig)
        hit_rows.append({
            "run": run,
            "query": query,
            "contig": contig,
            "contig_length_nt": summary.contig_length,
            "contig_sha256": sha256_sequence(assembly[contig].sequence),
            "single_contig_query_coverage": f"{summary.query_coverage:.6f}",
            "coordinate_weighted_identity": f"{summary.weighted_identity:.6f}",
            "max_internal_query_gap_nt": summary.max_internal_gap,
            "max_internal_subject_gap_nt": summary.max_internal_subject_gap,
            "relative_orientation": summary.orientation,
            "collinear_hsp_geometry": str(summary.collinear).lower(),
            "hsp_count": summary.hsp_count,
            "selected_as_best": str(selected.get(query) == contig).lower(),
        })
    hit_rows.sort(key=lambda row: (str(row["query"]), str(row["contig"])))
    hit_contigs = [
        record for identifier, record in assembly.items() if identifier in hit_contig_ids
    ]
    return AuditOutputs(
        recovery_rows=recovery_rows,
        hit_rows=hit_rows,
        hit_contigs=hit_contigs,
        technical_status="pass",
        failures=[],
        provenance=provenance,
        thresholds=thresholds,
    )


def declared_technical_failure_outputs(
    run: str,
    query_path: Path,
    structure_path: Path,
    threshold_path: Path,
    reasons: Sequence[str],
) -> AuditOutputs:
    """Create complete per-query rows after an assembler failure or timeout.

    Query, structure, and threshold inputs remain mandatory and are validated.
    Assembly, BLAST, and provenance inputs are intentionally not inspected in
    this path because they may not exist after a declared assembler failure.
    """
    if not reasons or any(not reason.strip() for reason in reasons):
        raise AuditError("empty_declared_technical_failure")
    queries = read_fasta(query_path, "query")
    thresholds = load_thresholds(threshold_path)
    structure = load_structure(structure_path, queries)
    rows: list[dict[str, object]] = []
    for query, record in queries.items():
        source_run = structure[query]["source_run"]
        role, coverage_threshold, identity_threshold = recovery_thresholds(
            run, source_run, thresholds
        )
        rows.append({
            "run": run,
            "query": query,
            "source_run": source_run,
            "run_role": role,
            "query_length_nt": len(record.sequence),
            "query_coverage_threshold": f"{coverage_threshold:.6f}",
            "identity_threshold": f"{identity_threshold:.6f}",
            "maximum_internal_query_gap_nt": thresholds[
                "maximum_internal_query_gap_nt"
            ],
            "maximum_internal_subject_gap_nt": thresholds[
                "maximum_internal_subject_gap_nt"
            ],
            "best_single_contig": "",
            "best_contig_length_nt": "",
            "best_single_contig_query_coverage": "0.000000",
            "coordinate_weighted_identity": "",
            "covered_query_nt": 0,
            "max_internal_query_gap_nt": "",
            "max_internal_subject_gap_nt": "",
            "relative_orientation": "",
            "collinear_hsp_geometry": "false",
            "hsp_count": 0,
            "best_evalue": "",
            "total_bitscore": "",
            "gate_status": "technical_incomplete",
            "recovery_status": f"{role}_technical_incomplete",
        })
    failures = [f"declared_assembly_technical_failure:{reason}" for reason in reasons]
    return AuditOutputs(
        recovery_rows=rows,
        hit_rows=[],
        hit_contigs=[],
        technical_status="technical_incomplete",
        failures=failures,
        provenance={"declared_technical_failures": list(reasons)},
        thresholds=thresholds,
    )


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, records: Sequence[FastaRecord]) -> None:
    with path.open("w") as handle:
        for record in records:
            handle.write(f">{record.description}\n")
            for index in range(0, len(record.sequence), 80):
                handle.write(record.sequence[index:index + 80] + "\n")


def write_checksums(out: Path) -> None:
    rows = []
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            rows.append(f"{sha256_path(path)}  {path.name}")
    (out / "SHA256SUMS.txt").write_text("\n".join(rows) + "\n")


def write_outputs(out: Path, run: str, result: AuditOutputs) -> dict[str, object]:
    write_tsv(out / "DE_NOVO_RECOVERY.tsv", RECOVERY_FIELDS, result.recovery_rows)
    write_tsv(out / "DE_NOVO_HIT_CONTIGS.tsv", HIT_FIELDS, result.hit_rows)
    write_fasta(out / "DE_NOVO_HIT_CONTIGS.fna", result.hit_contigs)
    recovered = sum(row.get("gate_status") == "pass" for row in result.recovery_rows)
    audit = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "run": run,
        "technical_status": result.technical_status,
        "technical_complete": result.technical_status == "pass",
        "failures": result.failures,
        "recovery_gate_status": (
            "technical_incomplete" if result.technical_status != "pass"
            else ("pass" if recovered else "fail")
        ),
        "recovered_query_count": recovered,
        "query_count": len(result.recovery_rows),
        "mapping_seed_free_provenance_validated": result.technical_status == "pass",
        "provenance": result.provenance,
        "workflow_defined_thresholds": result.thresholds,
        "candidate_hit_filter": {
            "minimum_hsp_length": MIN_HIT_HSP_LENGTH,
            "minimum_pident_percent": MIN_HIT_PIDENT,
            "maximum_evalue": MAX_HIT_EVALUE,
        },
        "claim_boundary": CLAIM_BOUNDARY,
    }
    (out / "DE_NOVO_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    report = [
        "# Mapping-seed-free de novo recovery audit", "",
        f"Technical status: **{audit['technical_status']}**", "",
        f"Run recovery status: **{audit['recovery_gate_status']}**", "",
        "| Query | Role | Best single contig | Query coverage | Identity | Gate |",
        "|---|---|---|---:|---:|---|",
    ]
    for row in result.recovery_rows:
        report.append(
            f"| `{row['query']}` | `{row['run_role']}` | "
            f"`{row['best_single_contig']}` | "
            f"{row['best_single_contig_query_coverage']} | "
            f"{row['coordinate_weighted_identity']} | `{row['gate_status']}` |"
        )
    report.extend(["", CLAIM_BOUNDARY])
    if result.failures:
        report.extend(["", "## Technical failures", ""])
        report.extend(f"- `{failure}`" for failure in result.failures)
    (out / "DE_NOVO_REPORT.md").write_text("\n".join(report) + "\n")
    write_checksums(out)
    return audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument(
        "--declared-technical-failure",
        action="append",
        default=[],
        metavar="REASON",
        help=(
            "Emit one technical_incomplete row per query without reading the "
            "assembly/BLAST/provenance inputs; intended for an assembler "
            "failure or timeout. May be repeated."
        ),
    )
    parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")
    try:
        if args.declared_technical_failure:
            result = declared_technical_failure_outputs(
                run=args.run,
                query_path=args.queries,
                structure_path=args.structure,
                threshold_path=args.thresholds,
                reasons=args.declared_technical_failure,
            )
        else:
            result = run_audit(
                run=args.run,
                query_path=args.queries,
                assembly_path=args.assembly,
                blast_path=args.blast,
                structure_path=args.structure,
                threshold_path=args.thresholds,
                provenance_path=args.provenance,
            )
    except AuditError as exc:
        result = AuditOutputs(
            recovery_rows=[],
            hit_rows=[],
            hit_contigs=[],
            technical_status="technical_incomplete",
            failures=[str(exc)],
            provenance={},
            thresholds={},
        )
    audit = write_outputs(args.out, args.run, result)
    print(json.dumps(audit, indent=2))
    return 0 if audit["technical_status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
