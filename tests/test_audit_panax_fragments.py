import contextlib
import csv
import io
import random
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from virus_hunt.finalization import audit_panax_fragments as audit_mod


def random_reference(seed: int, length: int = 620) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def alignment(
    qname: str,
    flag: int,
    reference: str,
    position: int,
    sequence: str,
    mapq: int = 60,
    cigar: str | None = None,
) -> audit_mod.Alignment:
    return audit_mod.Alignment(
        qname=qname,
        flag=flag,
        reference=reference,
        position=position,
        mapq=mapq,
        cigar=cigar or f"{len(sequence)}M",
        sequence=sequence,
    )


def pair(
    qname: str,
    reference: str,
    sequence: str,
    mapq1: int = 60,
    mapq2: int = 60,
    duplicate: bool = False,
) -> list[audit_mod.Alignment]:
    duplicate_flag = audit_mod.DUPLICATE if duplicate else 0
    return [
        alignment(
            qname,
            99 | duplicate_flag,
            reference,
            210,
            sequence[209:309],
            mapq1,
        ),
        alignment(
            qname,
            147 | duplicate_flag,
            reference,
            400,
            audit_mod.reverse_complement(sequence[399:499]),
            mapq2,
        ),
    ]


def sam_row(
    qname: str,
    flag: int,
    sequence: str,
    cigar: str | None = None,
    position: int = 1,
    mapq: int = 60,
    reference: str = "A",
) -> str:
    cigar = cigar or f"{len(sequence)}M"
    return (
        f"{qname}\t{flag}\t{reference}\t{position}\t{mapq}\t{cigar}"
        f"\t=\t0\t0\t{sequence}\t*\n"
    )


class AuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.references = {
            "A": random_reference(11),
            "B": random_reference(29),
        }

    def test_fragment_counts_duplicates_panel_kmers_and_continuity(self) -> None:
        groups = {
            "good": pair("good", "A", self.references["A"]),
            "duplicate": pair(
                "duplicate", "A", self.references["A"], duplicate=True
            ),
            "low_mapq": pair(
                "low_mapq", "A", self.references["A"], mapq2=29
            ),
        }
        result = audit_mod.audit(
            self.references, groups, 31, "RUN1", "Alias1"
        )
        by_reference = {row["reference"]: row for row in result.metrics}
        a_row = by_reference["A"]
        self.assertEqual(a_row["minimum_mate_mapq"], 30)
        self.assertEqual(a_row["proper_fragments_mapq30_preduplicate"], 2)
        self.assertEqual(a_row["proper_fragments_mapq30_nonduplicate"], 1)
        self.assertEqual(a_row["duplicate_fragment_fraction"], "0.500000")
        self.assertEqual(a_row["panel_unique_kmer_fragments_nonduplicate"], 1)
        self.assertGreater(a_row["panel_unique_31mers_observed"], 0)
        self.assertEqual(a_row["distinct_fragment_endpoints"], 1)
        self.assertEqual(a_row["fr_orientation_fraction"], "1.000000")
        self.assertEqual(a_row["template_span_median"], "290.000")
        self.assertEqual(result.preduplicate_qnames, {"good", "duplicate"})
        self.assertEqual(result.nonduplicate_qnames, {"good"})
        self.assertEqual(len(result.hashes), 1)
        self.assertEqual(len(result.hashes[0]["fragment_fingerprint_sha256"]), 64)

        continuity = {
            (row["reference"], row["boundary_after_nt"]): row
            for row in result.continuity
        }
        boundary = continuity[("A", 250)]
        self.assertEqual(
            boundary["proper_fragments_mapq30_nonduplicate_spanning"], 1
        )
        self.assertEqual(boundary["direct_read_spanning_fragments_anchor25"], 1)
        self.assertEqual(boundary["direct_read_spanning_reads_anchor25"], 1)
        self.assertEqual(
            continuity[("A", 500)][
                "proper_fragments_mapq30_nonduplicate_spanning"
            ],
            0,
        )

    def test_minimum_mapq_applies_to_both_mates_and_changes_field_name(self) -> None:
        groups = {"q": pair("q", "A", self.references["A"], mapq1=25, mapq2=60)}
        default_result = audit_mod.audit(
            self.references, groups, 31, "RUN", "Alias"
        )
        self.assertEqual(
            default_result.metrics[0]["proper_fragments_mapq30_preduplicate"], 0
        )
        result = audit_mod.audit(
            self.references, groups, 31, "RUN", "Alias", minimum_mapq=20
        )
        self.assertEqual(result.metrics[0]["minimum_mate_mapq"], 20)
        self.assertEqual(result.metrics[0]["proper_fragments_mapq20_preduplicate"], 1)
        self.assertIn(
            "proper_fragments_mapq20_nonduplicate_spanning",
            result.continuity[0],
        )

    def test_hashes_are_mate_order_and_orientation_insensitive(self) -> None:
        original_pair = pair("q", "A", self.references["A"])
        original = audit_mod.audit(
            self.references, {"q": original_pair}, 31, "RUN", "Alias"
        ).hashes[0]
        reversed_pair = [
            alignment(
                record.qname,
                record.flag,
                record.reference,
                record.position,
                audit_mod.reverse_complement(record.sequence),
                record.mapq,
                record.cigar,
            )
            for record in reversed(original_pair)
        ]
        reversed_result = audit_mod.audit(
            self.references, {"q": reversed_pair}, 31, "RUN", "Alias"
        ).hashes[0]
        self.assertEqual(
            original["pair_sequence_sha256"],
            reversed_result["pair_sequence_sha256"],
        )
        self.assertEqual(
            original["fragment_fingerprint_sha256"],
            reversed_result["fragment_fingerprint_sha256"],
        )

        moved = list(original_pair)
        moved[1] = alignment(
            moved[1].qname,
            moved[1].flag,
            moved[1].reference,
            moved[1].position + 1,
            moved[1].sequence,
            moved[1].mapq,
            moved[1].cigar,
        )
        moved_hash = audit_mod.audit(
            self.references, {"q": moved}, 31, "RUN", "Alias"
        ).hashes[0]
        self.assertEqual(original["pair_sequence_sha256"], moved_hash["pair_sequence_sha256"])
        self.assertNotEqual(original["endpoint_sha256"], moved_hash["endpoint_sha256"])
        self.assertNotEqual(
            original["fragment_fingerprint_sha256"],
            moved_hash["fragment_fingerprint_sha256"],
        )

    def test_multiple_primary_alignment_for_one_mate_is_rejected(self) -> None:
        records = pair("bad", "A", self.references["A"])
        records.append(records[0])
        with self.assertRaisesRegex(ValueError, "multiple primary alignments"):
            audit_mod.audit(
                self.references, {"bad": records}, 31, "RUN", "Alias"
            )


class CigarAndContinuityTests(unittest.TestCase):
    def test_soft_clips_are_found_inside_outer_hard_clips(self) -> None:
        record = alignment(
            "q", 99, "A", 201, "A" * 110, cigar="5H10S90M10S5H"
        )
        self.assertEqual(record.soft_clip, 20)
        self.assertEqual(record.end, 290)
        self.assertTrue(audit_mod.direct_read_spans_boundary(record, 250))

    def test_deletion_breaks_direct_anchor_but_insertion_does_not(self) -> None:
        deletion = alignment(
            "q", 99, "A", 226, "A" * 50, cigar="25M1D25M"
        )
        insertion = alignment(
            "q", 99, "A", 226, "A" * 51, cigar="25M1I25M"
        )
        self.assertFalse(audit_mod.direct_read_spans_boundary(deletion, 250))
        self.assertTrue(audit_mod.direct_read_spans_boundary(insertion, 250))


class PanelKmerTests(unittest.TestCase):
    def test_shared_reverse_complement_ambiguous_and_low_complexity_are_excluded(self) -> None:
        shared = "ACGTTGCAAGTCGAT"
        references = {
            "A": shared + "NNNNN" + "AACCGGTTACGATGCTAGCA",
            "B": audit_mod.reverse_complement(shared) + "GTCAGTACCGATGCATCGTA",
            "C": "A" * 20 + "ACGTGCACTGATCGTAGCTA",
        }
        kmers = audit_mod.panel_unique_kmers(references, 15)
        self.assertNotIn(audit_mod.canonical(shared), kmers["A"])
        self.assertNotIn(audit_mod.canonical(shared), kmers["B"])
        self.assertTrue(all("N" not in kmer for values in kmers.values() for kmer in values))
        self.assertNotIn(audit_mod.canonical("A" * 15), kmers["C"])
        self.assertGreater(len(kmers["A"]), 0)
        self.assertGreater(len(kmers["B"]), 0)


class SamValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.references = {"A": random_reference(37)}
        self.tempdir = tempfile.TemporaryDirectory()
        self.sam_path = Path(self.tempdir.name) / "input.sam"

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def parse_row(self, row: str) -> dict[str, list[audit_mod.Alignment]]:
        self.sam_path.write_text(row)
        return audit_mod.parse_sam(str(self.sam_path), self.references)

    def test_single_record_with_both_mate_bits_is_rejected(self) -> None:
        row = sam_row("bad", 1 | 2 | 64 | 128, "A" * 50)
        with self.assertRaisesRegex(ValueError, "exactly one READ1/READ2"):
            self.parse_row(row)

    def test_fixmate_normalized_mapped_singletons_are_ignored(self) -> None:
        rows = "".join([
            sam_row("forward", 0, "A" * 50),
            sam_row("read_designated", audit_mod.READ1, "A" * 50),
            sam_row(
                "reverse_duplicate",
                audit_mod.REVERSE | audit_mod.DUPLICATE,
                "A" * 50,
            ),
        ])
        self.assertEqual(self.parse_row(rows), {})

    def test_unpaired_record_does_not_change_a_valid_pair(self) -> None:
        records = pair("good", "A", self.references["A"])
        rows = "".join(
            sam_row(
                record.qname,
                record.flag,
                record.sequence,
                cigar=record.cigar,
                position=record.position,
                mapq=record.mapq,
                reference=record.reference,
            )
            for record in records
        ) + sam_row("singleton", 0, "A" * 50)
        groups = self.parse_row(rows)
        result = audit_mod.audit(self.references, groups, 31, "RUN", "Alias")
        self.assertEqual(result.preduplicate_qnames, {"good"})
        self.assertEqual(result.nonduplicate_qnames, {"good"})
        self.assertEqual(len(result.hashes), 1)

    def test_same_qname_paired_and_unpaired_records_are_rejected(self) -> None:
        records = pair("collision", "A", self.references["A"])
        rows = "".join(
            sam_row(
                record.qname,
                record.flag,
                record.sequence,
                cigar=record.cigar,
                position=record.position,
                mapq=record.mapq,
                reference=record.reference,
            )
            for record in records
        ) + sam_row("collision", 0, "A" * 50)
        with self.assertRaisesRegex(ValueError, "both paired and unpaired"):
            self.parse_row(rows)

    def test_malformed_cigar_and_sequence_length_are_rejected(self) -> None:
        cases = [
            sam_row("bad", 99, "A" * 10, cigar="10Mgarbage"),
            sam_row("bad", 99, "A" * 10, cigar="0M"),
            sam_row("bad", 99, "A" * 10, cigar="9M"),
            sam_row("bad", 99, "A" * 100, cigar="10S5H90M"),
        ]
        for row in cases:
            with self.subTest(row=row):
                with self.assertRaises(ValueError):
                    self.parse_row(row)

    def test_out_of_bounds_and_invalid_mapq_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside reference bounds"):
            self.parse_row(sam_row("bad", 99, "A" * 100, position=550))
        with self.assertRaisesRegex(ValueError, "MAPQ outside"):
            self.parse_row(sam_row("bad", 99, "A" * 100, mapq=256))

    def test_qc_fail_primary_is_excluded(self) -> None:
        result = self.parse_row(
            sam_row("qc", 99 | audit_mod.QCFAIL, "A" * 100)
        )
        self.assertEqual(result, {})


class OutputTests(unittest.TestCase):
    def test_empty_hash_table_has_header_and_no_sentinel_row(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "hashes.tsv"
            audit_mod.write_tsv(path, [], audit_mod.HASH_FIELDS)
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, audit_mod.HASH_FIELDS)
                self.assertEqual(list(reader), [])

    def test_qname_output_is_sorted_and_can_be_empty(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            nonempty = Path(tempdir) / "qnames.txt"
            empty = Path(tempdir) / "empty.txt"
            audit_mod.write_qnames(nonempty, {"z", "a"})
            audit_mod.write_qnames(empty, set())
            self.assertEqual(nonempty.read_text(), "a\nz\n")
            self.assertEqual(empty.read_text(), "")

    def test_cli_writes_metrics_hashes_continuity_and_both_qname_sets(self) -> None:
        references = {"A": random_reference(41), "B": random_reference(43)}
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "references.fna"
            sam = root / "reads.sam"
            out = root / "out"
            fasta.write_text(
                "".join(f">{name}\n{sequence}\n" for name, sequence in references.items())
            )
            records = pair("q1", "A", references["A"])
            sam.write_text(
                "".join(
                    sam_row(
                        record.qname,
                        record.flag,
                        record.sequence,
                        cigar=record.cigar,
                        position=record.position,
                        mapq=record.mapq,
                        reference=record.reference,
                    )
                    for record in records
                )
            )
            argv = [
                "audit_panax_fragments.py",
                "--sam", str(sam),
                "--reference", str(fasta),
                "--run", "RUN1",
                "--alias", "Alias1",
                "--out", str(out),
                "--minimum-mapq", "30",
                "--write-qualifying-qnames",
            ]
            with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(audit_mod.main(), 0)

            self.assertEqual(
                (out / "QUALIFYING_QNAMES_PRE_DUPLICATE.txt").read_text(), "q1\n"
            )
            self.assertEqual(
                (out / "QUALIFYING_QNAMES_NON_DUPLICATE.txt").read_text(), "q1\n"
            )
            with (out / "FRAGMENT_METRICS.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["minimum_mate_mapq"], "30")
            self.assertEqual(rows[0]["proper_fragments_mapq30_nonduplicate"], "1")
            with (out / "FRAGMENT_HASHES.tsv").open(newline="") as handle:
                hashes = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(hashes), 1)
            self.assertEqual(len(hashes[0]["fragment_fingerprint_sha256"]), 64)
            with (out / "CONTINUITY_SCAN.tsv").open(newline="") as handle:
                continuity = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(continuity)
            self.assertEqual(continuity[0]["minimum_mate_mapq"], "30")


if __name__ == "__main__":
    unittest.main()
