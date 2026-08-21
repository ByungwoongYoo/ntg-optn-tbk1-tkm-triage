from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FINALIZATION = ROOT / "virus_hunt" / "finalization"
sys.path.insert(0, str(FINALIZATION))

import audit_panax_denovo as denovo  # noqa: E402


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier} fixture\n{sequence}\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class DenovoFixture:
    run = "DRR853910"

    def __init__(self, root: Path):
        self.root = root
        self.queries = root / "queries.fna"
        self.assembly = root / "assembly.fna"
        self.blast = root / "blast.tsv"
        self.structure = root / "structure.tsv"
        self.thresholds = root / "thresholds.json"
        self.provenance = root / "provenance.json"
        self.fastq_manifest = root / "FASTQ_SHA256.txt"
        self.method_manifest = root / "ASSEMBLY_METHOD.tsv"
        self.fastq_manifest.write_text("fixture fastq manifest\n")
        self.method_manifest.write_text("fixture assembly method\n")
        self.query_sequences = {"Query_A": "A" * 1000, "Query_B": "C" * 800}
        self.assembly_sequences = {
            "contig_A": "A" * 1200,
            "contig_B": "C" * 900,
            "contig_reverse": "G" * 900,
            "contig_gap": "A" * 2000,
            "contig_bridge": "A" * 6000,
            "unused": "G" * 500,
        }
        write_fasta(self.queries, self.query_sequences)
        write_fasta(self.assembly, self.assembly_sequences)
        self.structure.write_text(
            "reference\tsource_run\tsource_contig\treference_length_nt\n"
            "Query_A\tDRR853910\tsource_A\t1000\n"
            "Query_B\tDRR853912\tsource_B\t800\n"
        )
        self.thresholds.write_text(json.dumps({
            "scope": "workflow-defined technical rules",
            "de_novo_recovery": {
                "source_run_single_contig_query_coverage_minimum": 0.90,
                "source_run_identity_minimum": 0.98,
                "non_source_single_contig_query_coverage_minimum": 0.70,
                "non_source_identity_minimum": 0.90,
                "maximum_internal_query_gap_nt": 150,
                "maximum_internal_subject_gap_nt": 500,
                "candidate_support_rule": (
                    "at least one source-run recovery and at least one "
                    "non-source-run recovery"
                ),
            },
        }))
        self.write_provenance()

    def write_provenance(self, **overrides: object) -> None:
        value: dict[str, object] = {
            "run": self.run,
            "assembly_sha256": denovo.sha256_path(self.assembly),
            "assembly_file": self.assembly.name,
            "retained_assembly_file": "FULL_ASSEMBLY.fna.gz",
            "retained_assembly_manifest": "FULL_ASSEMBLY_MANIFEST.tsv",
            "retained_assembly_compression": "gzip_nondeterministic_metadata_disabled",
            "input_scope": "complete_paired_fastq",
            "candidate_baiting": False,
            "mapping_seeded": False,
            "reference_guided": False,
            "target_read_selection": False,
            "assembler": "MEGAHIT",
            "assembler_version": "1.2.9",
            "assembly_exit_code": "0",
            "search_exit_code": "0",
            "fastq_sha256_manifest": self.fastq_manifest.name,
            "fastq_sha256_manifest_sha256": denovo.sha256_path(self.fastq_manifest),
            "assembly_method_manifest": self.method_manifest.name,
            "assembly_method_manifest_sha256": denovo.sha256_path(self.method_manifest),
            "candidate_query_sha256": denovo.sha256_path(self.queries),
            "candidate_structure_sha256": denovo.sha256_path(self.structure),
        }
        value.update(overrides)
        self.provenance.write_text(json.dumps(value))

    def blast_row(
        self,
        query: str,
        contig: str,
        pident: float,
        length: int,
        qlen: int,
        slen: int,
        qstart: int,
        qend: int,
        sstart: int,
        send: int,
        evalue: str = "1e-50",
        bitscore: float = 500.0,
        aligned_query: str | None = None,
        aligned_subject: str | None = None,
    ) -> str:
        if aligned_query is None or aligned_subject is None:
            query_start, query_end = sorted((qstart, qend))
            subject_start, subject_end = sorted((sstart, send))
            aligned_query = self.query_sequences[query][query_start - 1:query_end]
            aligned_subject = self.assembly_sequences[contig][subject_start - 1:subject_end]
            if qend < qstart:
                aligned_query = denovo.reverse_complement(aligned_query)
            if send < sstart:
                aligned_subject = denovo.reverse_complement(aligned_subject)
            if len(aligned_query) != length or len(aligned_subject) != length:
                # Deliberately malformed fixtures can still exercise the
                # production parser's fail-closed span validation.
                aligned_query = aligned_query.ljust(length, "A")[:length]
                aligned_subject = aligned_subject.ljust(length, "A")[:length]
            nident = sum(
                query_base == subject_base and query_base != "-"
                for query_base, subject_base in zip(aligned_query, aligned_subject)
            )
            mismatch = sum(
                query_base != subject_base
                and query_base != "-" and subject_base != "-"
                for query_base, subject_base in zip(aligned_query, aligned_subject)
            )
        else:
            if len(aligned_query) != length or len(aligned_subject) != length:
                raise ValueError("aligned fixture strings must match alignment length")
            nident = sum(
                query_base == subject_base and query_base != "-"
                for query_base, subject_base in zip(aligned_query, aligned_subject)
            )
            mismatch = sum(
                query_base != subject_base
                and query_base != "-" and subject_base != "-"
                for query_base, subject_base in zip(aligned_query, aligned_subject)
            )
        gaps = sum(
            query_base == "-" or subject_base == "-"
            for query_base, subject_base in zip(aligned_query, aligned_subject)
        )
        gapopen = 0
        previous_gap_side = ""
        for query_base, subject_base in zip(aligned_query, aligned_subject):
            gap_side = "query" if query_base == "-" else (
                "subject" if subject_base == "-" else ""
            )
            if gap_side and gap_side != previous_gap_side:
                gapopen += 1
            previous_gap_side = gap_side
        pident = 100 * nident / length
        return "\t".join(map(str, [
            query, contig, pident, length, nident, mismatch, gapopen, gaps,
            qlen, slen, qstart, qend, sstart, send,
            aligned_query, aligned_subject, evalue, bitscore,
        ]))

    def args(self, out: Path, extra: list[str] | None = None) -> list[str]:
        return [
            "--run", self.run,
            "--queries", str(self.queries),
            "--assembly", str(self.assembly),
            "--blast", str(self.blast),
            "--structure", str(self.structure),
            "--thresholds", str(self.thresholds),
            "--provenance", str(self.provenance),
            "--out", str(out),
            *(extra or []),
        ]


class AuditPanaxDenovoTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = DenovoFixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def invoke(self, out: Path, extra: list[str] | None = None) -> int:
        with redirect_stdout(io.StringIO()):
            return denovo.main(self.fixture.args(out, extra))

    def test_source_and_non_source_single_contig_recovery(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_A", 99.0, 500, 1000, 1200,
                1, 500, 1, 500, bitscore=600,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_A", 98.2, 450, 1000, 1200,
                501, 950, 501, 950, bitscore=550,
            ),
            self.fixture.blast_row(
                "Query_B", "contig_B", 92.0, 600, 800, 900,
                1, 600, 1, 600, bitscore=500,
            ),
        ]) + "\n")
        out = self.root / "out"
        self.assertEqual(self.invoke(out), 0)

        rows = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}
        self.assertEqual(rows["Query_A"]["run_role"], "source")
        self.assertEqual(rows["Query_A"]["gate_status"], "pass")
        self.assertEqual(rows["Query_A"]["best_single_contig"], "contig_A")
        self.assertEqual(rows["Query_A"]["best_single_contig_query_coverage"], "0.950000")
        self.assertEqual(rows["Query_B"]["run_role"], "non_source")
        self.assertEqual(rows["Query_B"]["gate_status"], "pass")
        hit_fasta = (out / "DE_NOVO_HIT_CONTIGS.fna").read_text()
        self.assertIn(">contig_A fixture", hit_fasta)
        self.assertIn(">contig_B fixture", hit_fasta)
        self.assertNotIn(">unused", hit_fasta)
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "pass")
        self.assertTrue(audit["mapping_seed_free_provenance_validated"])

    def test_empty_blast_is_complete_but_not_recovered(self) -> None:
        self.fixture.blast.write_text("")
        out = self.root / "empty"
        self.assertEqual(self.invoke(out), 0)
        rows = read_tsv(out / "DE_NOVO_RECOVERY.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["gate_status"] for row in rows}, {"fail"})
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "pass")
        self.assertEqual(audit["recovery_gate_status"], "fail")
        self.assertEqual((out / "DE_NOVO_HIT_CONTIGS.fna").read_text(), "")

    def test_mixed_orientation_hsps_fail_collinearity(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_A", 99.0, 500, 1000, 1200,
                1, 500, 1, 500,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_A", 99.0, 450, 1000, 1200,
                501, 950, 1100, 651,
            ),
        ]) + "\n")
        out = self.root / "mixed"
        self.assertEqual(self.invoke(out), 0)
        row = read_tsv(out / "DE_NOVO_RECOVERY.tsv")[0]
        self.assertEqual(row["relative_orientation"], "mixed")
        self.assertEqual(row["collinear_hsp_geometry"], "false")
        self.assertEqual(row["gate_status"], "fail")

    def test_internal_query_gap_over_workflow_threshold_fails(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_B", "contig_B", 92.0, 300, 800, 900,
                1, 300, 1, 300,
            ),
            self.fixture.blast_row(
                "Query_B", "contig_B", 92.0, 340, 800, 900,
                461, 800, 461, 800,
            ),
        ]) + "\n")
        out = self.root / "large_gap"
        self.assertEqual(self.invoke(out), 0)
        rows = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}
        self.assertEqual(rows["Query_B"]["run_role"], "non_source")
        self.assertEqual(rows["Query_B"]["best_single_contig_query_coverage"], "0.800000")
        self.assertEqual(rows["Query_B"]["max_internal_query_gap_nt"], "160")
        self.assertEqual(rows["Query_B"]["maximum_internal_query_gap_nt"], "150")
        self.assertEqual(rows["Query_B"]["gate_status"], "fail")

    def test_internal_subject_gap_over_workflow_threshold_fails(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 475, 1000, 2000,
                1, 475, 1, 475,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 475, 1000, 2000,
                476, 950, 977, 1451,
            ),
        ]) + "\n")
        out = self.root / "large_subject_gap"
        self.assertEqual(self.invoke(out), 0)
        rows = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}
        self.assertEqual(rows["Query_A"]["best_single_contig_query_coverage"], "0.950000")
        self.assertEqual(rows["Query_A"]["max_internal_query_gap_nt"], "0")
        self.assertEqual(rows["Query_A"]["max_internal_subject_gap_nt"], "501")
        self.assertEqual(rows["Query_A"]["maximum_internal_subject_gap_nt"], "500")
        self.assertEqual(rows["Query_A"]["gate_status"], "fail")

    def test_overlapping_hsps_cannot_bridge_subject_gap(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.5, 450, 1000, 2000,
                1, 450, 1, 450,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.5, 501, 1000, 2000,
                200, 700, 451, 951,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.5, 450, 1000, 2000,
                451, 900, 952, 1401,
            ),
        ]) + "\n")
        out = self.root / "overlapping_subject_bridge"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_A"]
        self.assertEqual(row["best_single_contig_query_coverage"], "0.900000")
        self.assertEqual(row["max_internal_subject_gap_nt"], "501")
        self.assertEqual(row["hsp_count"], "2")
        self.assertEqual(row["gate_status"], "fail")

    def test_repeated_overlapping_hsps_cannot_create_false_full_chain(self) -> None:
        rows = [
            self.fixture.blast_row(
                "Query_A", "contig_bridge", 99.5, 400, 1000, 6000,
                1, 400, 1, 400,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_bridge", 99.5, 400, 1000, 6000,
                601, 1000, 5001, 5400,
            ),
        ]
        rows.extend(
            self.fixture.blast_row(
                "Query_A", "contig_bridge", 99.5, 401, 1000, 6000,
                300, 700, subject_start, subject_start + 400,
            )
            for subject_start in range(401, 4602, 400)
        )
        self.fixture.blast.write_text("\n".join(rows) + "\n")
        out = self.root / "repeated_overlap_bridge"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_A"]
        self.assertLess(float(row["best_single_contig_query_coverage"]), 0.9)
        self.assertEqual(row["gate_status"], "fail")

    def test_impossible_alignment_span_is_technical_incomplete(self) -> None:
        self.fixture.blast.write_text(self.fixture.blast_row(
            "Query_A", "contig_gap", 99.0, 950, 1000, 2000,
            1, 950, 1, 450,
        ) + "\n")
        out = self.root / "impossible_alignment_span"
        self.assertEqual(self.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("blast_subject_span_sequence_mismatch", audit["failures"][0])

    def test_aligned_strings_must_match_fasta_coordinates(self) -> None:
        self.fixture.blast.write_text(self.fixture.blast_row(
            "Query_A", "contig_B", 99.0, 900, 1000, 900,
            1, 900, 1, 900,
            aligned_query="A" * 900,
            aligned_subject="A" * 900,
        ) + "\n")
        out = self.root / "forged_aligned_subject"
        self.assertEqual(self.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("blast_subject_sequence_mismatch", audit["failures"][0])

    def test_gapped_alignment_uses_residue_coverage_and_gap_runs(self) -> None:
        aligned_query = "C" * 300 + "-" * 10 + "C" * 300
        aligned_subject = "C" * 150 + "-" * 10 + "C" * 450
        self.fixture.blast.write_text(self.fixture.blast_row(
            "Query_B", "contig_B", 0.0, 610, 800, 900,
            1, 600, 1, 600,
            aligned_query=aligned_query,
            aligned_subject=aligned_subject,
        ) + "\n")
        out = self.root / "gapped_alignment"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_B"]
        self.assertEqual(row["covered_query_nt"], "590")
        self.assertEqual(row["best_single_contig_query_coverage"], "0.737500")
        self.assertEqual(row["max_internal_query_gap_nt"], "10")
        self.assertEqual(row["max_internal_subject_gap_nt"], "10")
        self.assertEqual(row["gate_status"], "pass")

    def test_three_contiguous_hsps_have_zero_subject_gap(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 300, 1000, 2000,
                1, 300, 1, 300,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 300, 1000, 2000,
                301, 600, 301, 600,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 350, 1000, 2000,
                601, 950, 601, 950,
            ),
        ]) + "\n")
        out = self.root / "three_contiguous_hsps"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_A"]
        self.assertEqual(row["max_internal_subject_gap_nt"], "0")
        self.assertEqual(row["gate_status"], "pass")

    def test_contiguous_reverse_orientation_chain_passes(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_B", "contig_reverse", 92.0, 300, 800, 900,
                1, 300, 900, 601,
            ),
            self.fixture.blast_row(
                "Query_B", "contig_reverse", 92.0, 300, 800, 900,
                301, 600, 600, 301,
            ),
        ]) + "\n")
        out = self.root / "reverse_chain"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_B"]
        self.assertEqual(row["relative_orientation"], "minus")
        self.assertEqual(row["max_internal_subject_gap_nt"], "0")
        self.assertEqual(row["best_single_contig_query_coverage"], "0.750000")
        self.assertEqual(row["gate_status"], "pass")

    def test_qualifying_contig_is_selected_over_higher_coverage_gapped_contig(self) -> None:
        self.fixture.blast.write_text("\n".join([
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 475, 1000, 2000,
                1, 475, 1, 475,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_gap", 99.0, 475, 1000, 2000,
                476, 950, 977, 1451,
            ),
            self.fixture.blast_row(
                "Query_A", "contig_A", 99.0, 900, 1000, 1200,
                1, 900, 1, 900,
            ),
        ]) + "\n")
        out = self.root / "eligible_best"
        self.assertEqual(self.invoke(out), 0)
        row = {row["query"]: row for row in read_tsv(out / "DE_NOVO_RECOVERY.tsv")}["Query_A"]
        self.assertEqual(row["best_single_contig"], "contig_A")
        self.assertEqual(row["gate_status"], "pass")

    def test_provenance_hash_mismatch_is_fail_closed(self) -> None:
        self.fixture.blast.write_text("")
        self.fixture.write_provenance(assembly_sha256="0" * 64)
        out = self.root / "bad_provenance"
        self.assertEqual(self.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("assembly_provenance_sha256_mismatch", audit["failures"][0])
        self.assertEqual(len(read_tsv(out / "DE_NOVO_RECOVERY.tsv")), 0)

    def test_declared_assembler_timeout_preserves_query_rows(self) -> None:
        self.fixture.blast.write_text("")
        self.fixture.assembly = self.root / "assembly_was_not_created.fna"
        self.fixture.provenance = self.root / "provenance_was_not_created.json"
        out = self.root / "timeout"
        self.assertEqual(
            self.invoke(out, ["--declared-technical-failure", "megahit_timeout"]),
            1,
        )
        rows = read_tsv(out / "DE_NOVO_RECOVERY.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            {row["gate_status"] for row in rows}, {"technical_incomplete"}
        )
        audit = json.loads((out / "DE_NOVO_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("megahit_timeout", audit["failures"][0])


if __name__ == "__main__":
    unittest.main()
