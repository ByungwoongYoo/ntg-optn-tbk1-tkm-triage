#!/usr/bin/env python3
"""Audit proper paired fragments, duplicates, and panel-unique exact k-mers.

Input SAM records must come from a candidate-only BAM after name repair and
samtools markdup. The script reports fragments rather than alignment records.
Sequence hashes are orientation-insensitive and omit read names.

``CONTINUITY_SCAN.tsv`` uses fixed boundaries after reference bases 250, 500,
750, and so on.  A fragment spans a boundary when its outer aligned endpoints
fall on opposite sides.  A direct-read spanner additionally requires one read
to align at least 25 reference bases on both sides without a deletion or
reference skip in either anchor.  These are technical continuity checks, not
biological detection thresholds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import statistics
import sys
from collections import defaultdict
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path


UNMAPPED = 0x4
REVERSE = 0x10
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0x100
QCFAIL = 0x200
DUPLICATE = 0x400
SUPPLEMENTARY = 0x800
PROPER_PAIR = 0x2
PAIRED = 0x1
EXCLUDED = UNMAPPED | SECONDARY | QCFAIL | SUPPLEMENTARY
DNA = re.compile(r"^[ACGTN]+$")
CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")
CIGAR_FULL = re.compile(r"(?:[1-9][0-9]*[MIDNSHP=X])+")
DEFAULT_MINIMUM_MAPQ = 30
CONTINUITY_STEP = 250
DIRECT_ANCHOR = 25

HASH_FIELDS = [
    "run", "alias", "reference", "minimum_mate_mapq",
    "pair_sequence_sha256", "endpoint_sha256", "fragment_fingerprint_sha256",
]


def continuity_fields(minimum_mapq: int) -> list[str]:
    return [
        "run", "alias", "reference", "reference_length",
        "minimum_mate_mapq", "boundary_after_nt",
        "direct_read_anchor_nt_each_side",
        f"proper_fragments_mapq{minimum_mapq}_nonduplicate_spanning",
        "direct_read_spanning_fragments_anchor25",
        "direct_read_spanning_reads_anchor25",
    ]


@dataclass(frozen=True)
class Alignment:
    qname: str
    flag: int
    reference: str
    position: int
    mapq: int
    cigar: str
    sequence: str

    @property
    def end(self) -> int:
        span = sum(
            int(length)
            for length, op in CIGAR_TOKEN.findall(self.cigar)
            if op in "MDN=X"
        )
        return self.position + span - 1

    @property
    def soft_clip(self) -> int:
        tokens = CIGAR_TOKEN.findall(self.cigar)
        left_index = 0
        while left_index < len(tokens) and tokens[left_index][1] in "HP":
            left_index += 1
        right_index = len(tokens) - 1
        while right_index >= 0 and tokens[right_index][1] in "HP":
            right_index -= 1
        left = (
            int(tokens[left_index][0])
            if left_index < len(tokens) and tokens[left_index][1] == "S"
            else 0
        )
        right = (
            int(tokens[right_index][0])
            if right_index >= 0 and tokens[right_index][1] == "S"
            else 0
        )
        return left + right

    @property
    def aligned_reference_intervals(self) -> list[tuple[int, int]]:
        """Return 1-based closed intervals covered by query bases."""
        reference_position = self.position
        intervals: list[tuple[int, int]] = []
        for length_text, op in CIGAR_TOKEN.findall(self.cigar):
            length = int(length_text)
            if op in "M=X":
                start = reference_position
                end = reference_position + length - 1
                # Insertions do not consume reference and therefore leave
                # adjacent reference-aligned blocks continuous.
                if intervals and intervals[-1][1] + 1 == start:
                    intervals[-1] = (intervals[-1][0], end)
                else:
                    intervals.append((start, end))
                reference_position += length
            elif op in "DN":
                reference_position += length
        return intervals


@dataclass(frozen=True)
class AuditResult:
    metrics: list[dict[str, object]]
    hashes: list[dict[str, object]]
    continuity: list[dict[str, object]]
    preduplicate_qnames: set[str]
    nonduplicate_qnames: set[str]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f"duplicate FASTA record: {name}")
            chunks = []
        elif name is None:
            raise ValueError("sequence before FASTA header")
        else:
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    if not records:
        raise ValueError("empty reference FASTA")
    for ref, sequence in records.items():
        if not DNA.fullmatch(sequence):
            raise ValueError(f"invalid reference sequence: {ref}")
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def canonical(sequence: str) -> str:
    rc = reverse_complement(sequence)
    return sequence if sequence <= rc else rc


def informative_panel_kmer(sequence: str) -> bool:
    """Apply a minimal, explicit ambiguity/low-complexity exclusion."""
    return "N" not in sequence and len(set(sequence)) >= 3


def panel_unique_kmers(references: dict[str, str], k: int) -> dict[str, set[str]]:
    """Return canonical A/C/G/T k-mers unique within the supplied panel.

    This deliberately makes no claim of uniqueness against a host, vector,
    or broader sequence database. K-mers containing ambiguous bases or fewer
    than three distinct nucleotide symbols are not eligible.
    """
    owners: dict[str, set[str]] = defaultdict(set)
    for reference, sequence in references.items():
        if len(sequence) < k:
            raise ValueError(f"reference is shorter than k={k}: {reference}")
        for index in range(len(sequence) - k + 1):
            kmer = sequence[index:index + k]
            if informative_panel_kmer(kmer):
                owners[canonical(kmer)].add(reference)
    output = {reference: set() for reference in references}
    for kmer, refs in owners.items():
        if len(refs) == 1:
            output[next(iter(refs))].add(kmer)
    return output


def open_sam(path: str):
    return nullcontext(sys.stdin) if path == "-" else open(path)


def parse_sam(path: str, references: dict[str, str]) -> dict[str, list[Alignment]]:
    groups: dict[str, list[Alignment]] = defaultdict(list)
    unpaired_qnames: set[str] = set()
    with open_sam(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("@"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM row {number}: {len(fields)} fields")
            try:
                flag = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid SAM flag at row {number}: {fields[1]!r}") from exc
            if not 0 <= flag <= 0xFFFF:
                raise ValueError(f"SAM flag outside 0..65535 at row {number}: {flag}")
            qname = fields[0]
            if not qname or any(character.isspace() for character in qname):
                raise ValueError(f"invalid SAM QNAME at row {number}: {qname!r}")
            reference = fields[2]
            if flag & EXCLUDED:
                continue
            if reference not in references:
                raise ValueError(f"unexpected mapped reference at row {number}: {reference}")
            paired = bool(flag & PAIRED)
            is_read1 = bool(flag & READ1)
            is_read2 = bool(flag & READ2)
            if paired and is_read1 == is_read2:
                raise ValueError(
                    f"primary record must have exactly one READ1/READ2 bit at row "
                    f"{number}: {fields[0]}"
                )
            sequence = fields[9].upper()
            if sequence == "*" or not DNA.fullmatch(sequence):
                raise ValueError(f"invalid SAM sequence at row {number}: {fields[0]}")
            cigar = fields[5]
            if not CIGAR_FULL.fullmatch(cigar):
                raise ValueError(f"invalid SAM CIGAR at row {number}: {cigar!r}")
            tokens = CIGAR_TOKEN.findall(cigar)
            clip_ops = [op for _, op in tokens]
            core_start = 0
            while core_start < len(clip_ops) and clip_ops[core_start] in "HS":
                core_start += 1
            core_end = len(clip_ops)
            while core_end > core_start and clip_ops[core_end - 1] in "HS":
                core_end -= 1
            if any(op in "HS" for op in clip_ops[core_start:core_end]):
                raise ValueError(f"internal clipping operation at row {number}: {cigar!r}")
            if (
                any(op == "H" for op in clip_ops[1:core_start])
                or any(op == "H" for op in clip_ops[core_end:-1])
            ):
                raise ValueError(f"hard clip must be outermost at row {number}: {cigar!r}")
            query_span = sum(int(length) for length, op in tokens if op in "MIS=X")
            reference_span = sum(int(length) for length, op in tokens if op in "MDN=X")
            if query_span != len(sequence):
                raise ValueError(
                    f"CIGAR/SEQ length mismatch at row {number}: "
                    f"query_span={query_span}, sequence={len(sequence)}"
                )
            if reference_span < 1:
                raise ValueError(f"CIGAR consumes no reference at row {number}: {cigar}")
            try:
                position = int(fields[3])
                mapq = int(fields[4])
            except ValueError as exc:
                raise ValueError(
                    f"invalid SAM coordinate or MAPQ at row {number}: "
                    f"POS={fields[3]!r}, MAPQ={fields[4]!r}"
                ) from exc
            if position < 1 or position + reference_span - 1 > len(references[reference]):
                raise ValueError(
                    f"alignment outside reference bounds at row {number}: "
                    f"{reference}:{position}-{position + reference_span - 1}"
                )
            if not 0 <= mapq <= 255:
                raise ValueError(f"MAPQ outside 0..255 at row {number}: {mapq}")
            if not paired:
                # Filtering unmapped records before samtools fixmate leaves a
                # legitimate mapped singleton when its unmapped mate was
                # removed. SAM does not assign meaning to pair-only flag bits
                # when PAIRED is unset, and this record cannot contribute to a
                # both-mate proper-fragment audit.
                unpaired_qnames.add(qname)
                continue
            groups[fields[0]].append(Alignment(
                qname=qname, flag=flag, reference=reference,
                position=position, mapq=mapq, cigar=cigar,
                sequence=sequence,
            ))
    collisions = set(groups).intersection(unpaired_qnames)
    if collisions:
        example = sorted(collisions)[0]
        raise ValueError(
            f"QNAME occurs in both paired and unpaired primary records: {example}"
        )
    return groups


def percentile(values: list[int], fraction: float) -> float | str:
    if not values:
        return ""
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, round((len(ordered) - 1) * fraction)))
    return ordered[index]


def orientation(left: Alignment, right: Alignment) -> str:
    ordered = sorted((left, right), key=lambda record: (record.position, record.end, record.flag))
    return "".join("R" if record.flag & REVERSE else "F" for record in ordered)


def covers_reference_window(alignment: Alignment, start: int, end: int) -> bool:
    """Return whether query-aligned bases continuously cover a closed window."""
    cursor = start
    for interval_start, interval_end in alignment.aligned_reference_intervals:
        if interval_end < cursor:
            continue
        if interval_start > cursor:
            return False
        cursor = max(cursor, interval_end + 1)
        if cursor > end:
            return True
    return False


def direct_read_spans_boundary(
    alignment: Alignment, boundary_after_nt: int, anchor: int = DIRECT_ANCHOR,
) -> bool:
    """Require query-aligned bases on both sides of a reference boundary."""
    start = boundary_after_nt - anchor + 1
    end = boundary_after_nt + anchor
    return start >= 1 and covers_reference_window(alignment, start, end)


def audit(
    references: dict[str, str], groups: dict[str, list[Alignment]], k: int,
    run: str, alias: str, minimum_mapq: int = DEFAULT_MINIMUM_MAPQ,
) -> AuditResult:
    if k < 1:
        raise ValueError("k must be positive")
    if not 0 <= minimum_mapq <= 255:
        raise ValueError("minimum MAPQ must be in 0..255")
    kmers = panel_unique_kmers(references, k)
    stats: dict[str, dict[str, object]] = {}
    continuity: dict[str, dict[int, dict[str, int]]] = {}
    for reference in references:
        stats[reference] = {
            "predup": 0, "nondup": 0, "diagnostic": 0,
            "observed_kmers": set(), "endpoints": set(), "spans": [],
            "fr": 0, "softclipped": 0, "max_softclip": 0,
        }
        continuity[reference] = {
            boundary: {"fragment": 0, "direct_fragment": 0, "direct_read": 0}
            for boundary in range(CONTINUITY_STEP, len(references[reference]), CONTINUITY_STEP)
        }
    hashes: list[dict[str, object]] = []
    preduplicate_qnames: set[str] = set()
    nonduplicate_qnames: set[str] = set()

    for qname, records in groups.items():
        first = [record for record in records if record.flag & READ1]
        second = [record for record in records if record.flag & READ2]
        if len(first) > 1 or len(second) > 1:
            raise ValueError(
                f"multiple primary alignments for one mate: {qname} "
                f"(READ1={len(first)}, READ2={len(second)})"
            )
        if len(first) != 1 or len(second) != 1:
            continue
        r1, r2 = first[0], second[0]
        if (
            r1.reference != r2.reference
            or r1.mapq < minimum_mapq or r2.mapq < minimum_mapq
            or not (r1.flag & PROPER_PAIR and r2.flag & PROPER_PAIR)
        ):
            continue
        reference = r1.reference
        bucket = stats[reference]
        bucket["predup"] = int(bucket["predup"]) + 1
        preduplicate_qnames.add(qname)
        if r1.flag & DUPLICATE or r2.flag & DUPLICATE:
            continue
        bucket["nondup"] = int(bucket["nondup"]) + 1
        nonduplicate_qnames.add(qname)

        left = min(r1.position, r2.position)
        right = max(r1.end, r2.end)
        orient = orientation(r1, r2)
        detailed_endpoints = sorted(
            f"{read.position}:{read.end}:"
            f"{'R' if read.flag & REVERSE else 'F'}:{read.cigar}"
            for read in (r1, r2)
        )
        endpoint_signature = reference + "\n" + "\n".join(detailed_endpoints)
        cast_endpoints = bucket["endpoints"]
        assert isinstance(cast_endpoints, set)
        cast_endpoints.add(endpoint_signature)
        cast_spans = bucket["spans"]
        assert isinstance(cast_spans, list)
        cast_spans.append(right - left + 1)
        if orient == "FR":
            bucket["fr"] = int(bucket["fr"]) + 1

        soft_clip = r1.soft_clip + r2.soft_clip
        if soft_clip:
            bucket["softclipped"] = int(bucket["softclipped"]) + 1
        bucket["max_softclip"] = max(int(bucket["max_softclip"]), soft_clip)

        observed: set[str] = set()
        allowed = kmers[reference]
        for read in (r1.sequence, r2.sequence):
            for index in range(max(0, len(read) - k + 1)):
                word = canonical(read[index:index + k])
                if word in allowed:
                    observed.add(word)
        if observed:
            bucket["diagnostic"] = int(bucket["diagnostic"]) + 1
            cast_observed = bucket["observed_kmers"]
            assert isinstance(cast_observed, set)
            cast_observed.update(observed)

        pair = sorted((canonical(r1.sequence), canonical(r2.sequence)))
        pair_hash = hashlib.sha256((pair[0] + "\n" + pair[1]).encode()).hexdigest()
        endpoint_hash = hashlib.sha256(endpoint_signature.encode()).hexdigest()
        fragment_hash = hashlib.sha256(
            (reference + "\n" + pair_hash + "\n" + endpoint_hash).encode()
        ).hexdigest()
        hashes.append({
            "run": run, "alias": alias, "reference": reference,
            "minimum_mate_mapq": minimum_mapq,
            "pair_sequence_sha256": pair_hash,
            "endpoint_sha256": endpoint_hash,
            "fragment_fingerprint_sha256": fragment_hash,
        })

        for boundary, continuity_bucket in continuity[reference].items():
            if left <= boundary < right:
                continuity_bucket["fragment"] += 1
            direct_reads = sum(
                direct_read_spans_boundary(read, boundary) for read in (r1, r2)
            )
            if direct_reads:
                continuity_bucket["direct_fragment"] += 1
                continuity_bucket["direct_read"] += direct_reads

    rows: list[dict[str, object]] = []
    threshold = f"mapq{minimum_mapq}"
    for reference in references:
        bucket = stats[reference]
        predup = int(bucket["predup"])
        nondup = int(bucket["nondup"])
        diagnostic = int(bucket["diagnostic"])
        observed = bucket["observed_kmers"]
        endpoints = bucket["endpoints"]
        spans = bucket["spans"]
        assert isinstance(observed, set) and isinstance(endpoints, set) and isinstance(spans, list)
        rows.append({
            "run": run,
            "alias": alias,
            "reference": reference,
            "reference_length_nt": len(references[reference]),
            "minimum_mate_mapq": minimum_mapq,
            f"proper_fragments_{threshold}_preduplicate": predup,
            f"proper_fragments_{threshold}_nonduplicate": nondup,
            "duplicate_fragment_fraction": f"{(predup - nondup) / predup:.6f}" if predup else "",
            "panel_unique_kmer_fragments_nonduplicate": diagnostic,
            "panel_unique_kmer_fragment_fraction": f"{diagnostic / nondup:.6f}" if nondup else "",
            f"panel_unique_{k}mers_available": len(kmers[reference]),
            f"panel_unique_{k}mers_observed": len(observed),
            f"panel_unique_{k}mer_fraction_observed": (
                f"{len(observed) / len(kmers[reference]):.6f}" if kmers[reference] else ""
            ),
            "distinct_fragment_endpoints": len(endpoints),
            "fr_orientation_fraction": f"{int(bucket['fr']) / nondup:.6f}" if nondup else "",
            "template_span_median": f"{statistics.median(spans):.3f}" if spans else "",
            "template_span_p10": percentile(spans, 0.10),
            "template_span_p90": percentile(spans, 0.90),
            "softclipped_fragment_fraction": (
                f"{int(bucket['softclipped']) / nondup:.6f}" if nondup else ""
            ),
            "max_softclipped_bases_per_fragment": bucket["max_softclip"],
        })
    hashes.sort(key=lambda row: (
        str(row["reference"]), str(row["pair_sequence_sha256"]), str(row["endpoint_sha256"])
    ))
    continuity_rows: list[dict[str, object]] = []
    for reference in references:
        for boundary, bucket in continuity[reference].items():
            continuity_rows.append({
                "run": run,
                "alias": alias,
                "reference": reference,
                "reference_length": len(references[reference]),
                "minimum_mate_mapq": minimum_mapq,
                "boundary_after_nt": boundary,
                "direct_read_anchor_nt_each_side": DIRECT_ANCHOR,
                f"proper_fragments_mapq{minimum_mapq}_nonduplicate_spanning": (
                    bucket["fragment"]
                ),
                "direct_read_spanning_fragments_anchor25": bucket["direct_fragment"],
                "direct_read_spanning_reads_anchor25": bucket["direct_read"],
            })
    return AuditResult(
        metrics=rows,
        hashes=hashes,
        continuity=continuity_rows,
        preduplicate_qnames=preduplicate_qnames,
        nonduplicate_qnames=nonduplicate_qnames,
    )


def write_tsv(
    path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None,
) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError(f"field names required for an empty table: {path}")
        fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_qnames(path: Path, qnames: set[str]) -> None:
    path.write_text("".join(f"{qname}\n" for qname in sorted(qnames)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True, help="SAM file or - for stdin")
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--run", required=True)
    parser.add_argument("--alias", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument(
        "--minimum-mapq", type=int, default=DEFAULT_MINIMUM_MAPQ,
        help="minimum MAPQ required independently for both mates (default: 30)",
    )
    parser.add_argument(
        "--write-qualifying-qnames", action="store_true",
        help="write preduplicate and nonduplicate QNAME lists for samtools view -N",
    )
    args = parser.parse_args()
    if args.k < 15:
        raise SystemExit("panel-unique k must be at least 15")
    if not 0 <= args.minimum_mapq <= 255:
        raise SystemExit("minimum MAPQ must be in 0..255")
    args.out.mkdir(parents=True, exist_ok=True)
    if any(args.out.iterdir()):
        raise SystemExit(f"output directory must be empty: {args.out}")
    references = read_fasta(args.reference)
    groups = parse_sam(args.sam, references)
    result = audit(
        references, groups, args.k, args.run, args.alias, args.minimum_mapq,
    )
    write_tsv(args.out / "FRAGMENT_METRICS.tsv", result.metrics)
    write_tsv(args.out / "FRAGMENT_HASHES.tsv", result.hashes, HASH_FIELDS)
    write_tsv(
        args.out / "CONTINUITY_SCAN.tsv", result.continuity,
        continuity_fields(args.minimum_mapq),
    )
    if args.write_qualifying_qnames:
        write_qnames(
            args.out / "QUALIFYING_QNAMES_PRE_DUPLICATE.txt",
            result.preduplicate_qnames,
        )
        write_qnames(
            args.out / "QUALIFYING_QNAMES_NON_DUPLICATE.txt",
            result.nonduplicate_qnames,
        )
    print("\n".join("\t".join(map(str, row.values())) for row in result.metrics))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
