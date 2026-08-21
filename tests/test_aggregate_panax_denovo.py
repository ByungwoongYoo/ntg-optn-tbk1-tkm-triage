from __future__ import annotations

import csv
import gzip
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

import aggregate_panax_denovo as aggregate  # noqa: E402
import audit_panax_denovo as denovo  # noqa: E402


QUERIES = {"Query_A": 1000, "Query_B": 800, "Query_C": 600}
SOURCES = {
    "Query_A": "DRR853910",
    "Query_B": "DRR853911",
    "Query_C": "DRR853912",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AggregateFixture:
    def __init__(self, root: Path):
        self.root = root
        self.collected = root / "collected"
        self.collected.mkdir()
        self.queries = root / "queries.fna"
        self.structure = root / "structure.tsv"
        self.thresholds = root / "thresholds.json"
        self.queries.write_text("".join(
            f">{query}\n{'A' * length}\n" for query, length in QUERIES.items()
        ))
        self.structure.write_text(
            "reference\tsource_run\treference_length_nt\n" + "".join(
                f"{query}\t{SOURCES[query]}\t{length}\n"
                for query, length in QUERIES.items()
            )
        )
        self.thresholds.write_text(json.dumps({
            "scope": "workflow-defined technical rules",
            "de_novo_recovery": {
                "source_run_single_contig_query_coverage_minimum": 0.9,
                "source_run_identity_minimum": 0.98,
                "non_source_single_contig_query_coverage_minimum": 0.7,
                "non_source_identity_minimum": 0.9,
                "maximum_internal_query_gap_nt": 150,
                "maximum_internal_subject_gap_nt": 500,
                "candidate_support_rule": (
                    "at least one source-run recovery and at least one "
                    "non-source-run recovery"
                ),
            },
        }))
        self.normalized = denovo.load_thresholds(self.thresholds)

    def write_run(
        self,
        run: str,
        *,
        technical_status: str = "pass",
        blocked_source: str | None = None,
        high_gap: tuple[str, str] | None = None,
    ) -> None:
        directory = self.collected / run
        directory.mkdir()
        de_novo_directory = directory / "de_novo"
        de_novo_directory.mkdir()
        fastq_manifest = directory / "FASTQ_SHA256.txt"
        method_manifest = de_novo_directory / "ASSEMBLY_METHOD.tsv"
        fastq_manifest.write_text("fixture fastq hashes\n")
        method_manifest.write_text("fixture assembly method\n")
        assembly_bytes = f">assembly_{run}\n{'A' * 100}\n".encode()
        retained = de_novo_directory / "FULL_ASSEMBLY.fna.gz"
        retained.write_bytes(gzip.compress(assembly_bytes, compresslevel=9, mtime=0))
        assembly_sha = hashlib.sha256(assembly_bytes).hexdigest()
        compressed_sha = hashlib.sha256(retained.read_bytes()).hexdigest()
        (de_novo_directory / "FULL_ASSEMBLY_MANIFEST.tsv").write_text(
            "uncompressed_sha256\tcompressed_sha256\t"
            "uncompressed_bytes\tcompressed_bytes\n"
            f"{assembly_sha}\t{compressed_sha}\t{len(assembly_bytes)}\t"
            f"{retained.stat().st_size}\n"
        )
        rows: list[dict[str, object]] = []
        for query, length in QUERIES.items():
            source = SOURCES[query]
            role = "source" if run == source else "non_source"
            coverage_threshold = 0.9 if role == "source" else 0.7
            identity_threshold = 0.98 if role == "source" else 0.9
            ordinary_pass = (
                (role == "source" or run == "DRR853907")
                and not (role == "source" and query == blocked_source)
            )
            incomplete = technical_status == "technical_incomplete"
            gap_failure = high_gap == (run, query)
            has_contig = incomplete or ordinary_pass or gap_failure
            gate = (
                "technical_incomplete" if incomplete
                else ("pass" if ordinary_pass and not gap_failure else "fail")
            )
            rows.append({
                "run": run,
                "query": query,
                "source_run": source,
                "run_role": role,
                "query_length_nt": length,
                "query_coverage_threshold": f"{coverage_threshold:.6f}",
                "identity_threshold": f"{identity_threshold:.6f}",
                "maximum_internal_query_gap_nt": 150,
                "maximum_internal_subject_gap_nt": 500,
                "best_single_contig": f"contig_{run}_{query}" if has_contig else "",
                "best_contig_length_nt": length + 100 if has_contig else "",
                "best_single_contig_query_coverage": "1.000000" if has_contig else "0.000000",
                "coordinate_weighted_identity": "1.000000" if has_contig else "",
                "covered_query_nt": length if has_contig else 0,
                "max_internal_query_gap_nt": 151 if gap_failure else (25 if has_contig else ""),
                "max_internal_subject_gap_nt": 25 if has_contig else "",
                "relative_orientation": "plus" if has_contig else "",
                "collinear_hsp_geometry": "true" if has_contig else "false",
                "hsp_count": 2 if has_contig else 0,
                "best_evalue": "1e-50" if has_contig else "",
                "total_bitscore": "500.000000" if has_contig else "",
                "gate_status": gate,
                "recovery_status": (
                    f"{role}_technical_incomplete" if incomplete
                    else f"{role}_{'recovered' if gate == 'pass' else 'not_recovered'}"
                ),
            })
        denovo.write_tsv(
            directory / "DE_NOVO_RECOVERY.tsv", denovo.RECOVERY_FIELDS, rows
        )
        provenance = ({
            "run": run,
            "assembly_sha256": assembly_sha,
            "assembly_file": "megahit/final.contigs.fa",
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
            "fastq_sha256_manifest": "../FASTQ_SHA256.txt",
            "fastq_sha256_manifest_sha256": denovo.sha256_path(fastq_manifest),
            "assembly_method_manifest": "ASSEMBLY_METHOD.tsv",
            "assembly_method_manifest_sha256": denovo.sha256_path(method_manifest),
            "candidate_query_sha256": denovo.sha256_path(self.queries),
            "candidate_structure_sha256": denovo.sha256_path(self.structure),
        } if technical_status == "pass" else {
            "declared_technical_failures": ["megahit_timeout"]
        })
        if technical_status == "pass":
            (de_novo_directory / "ASSEMBLY_PROVENANCE.json").write_text(
                json.dumps(provenance)
            )
        recovered = sum(row["gate_status"] == "pass" for row in rows)
        audit = {
            "run": run,
            "technical_status": technical_status,
            "technical_complete": technical_status == "pass",
            "mapping_seed_free_provenance_validated": technical_status == "pass",
            "failures": (
                [] if technical_status == "pass" else ["megahit_timeout"]
            ),
            "recovery_gate_status": (
                "technical_incomplete" if technical_status != "pass"
                else ("pass" if recovered else "fail")
            ),
            "recovered_query_count": recovered,
            "query_count": len(rows),
            "workflow_defined_thresholds": self.normalized,
            "provenance": provenance,
        }
        (directory / "DE_NOVO_AUDIT.json").write_text(json.dumps(audit))

    def write_six(self, **kwargs: object) -> None:
        for run in aggregate.EXPECTED_RUNS:
            self.write_run(run, **kwargs)

    def invoke(self, out: Path) -> int:
        with redirect_stdout(io.StringIO()):
            return aggregate.main([
                "--collected", str(self.collected),
                "--queries", str(self.queries),
                "--structure", str(self.structure),
                "--thresholds", str(self.thresholds),
                "--out", str(out),
            ])


class AggregatePanaxDenovoTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = AggregateFixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_exactly_six_complete_runs_pass_candidate_rule(self) -> None:
        self.fixture.write_six()
        out = self.root / "out"
        self.assertEqual(self.fixture.invoke(out), 0)
        self.assertEqual(len(read_tsv(out / "ALL_DE_NOVO_RECOVERY.tsv")), 18)
        self.assertEqual(
            {row["gate_status"] for row in read_tsv(out / "DE_NOVO_CANDIDATE_GATE.tsv")},
            {"pass"},
        )

    def test_missing_run_is_structural_technical_incomplete(self) -> None:
        for run in aggregate.EXPECTED_RUNS[:-1]:
            self.fixture.write_run(run)
        out = self.root / "missing"
        self.assertEqual(self.fixture.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AGGREGATE_AUDIT.json").read_text())
        self.assertEqual(audit["overall_gate_status"], "technical_incomplete")
        self.assertIn("per_run_recovery_table_count:5!=6", audit["failures"][0])

    def test_source_recovery_is_required(self) -> None:
        self.fixture.write_six(blocked_source="Query_A")
        out = self.root / "source_fail"
        self.assertEqual(self.fixture.invoke(out), 1)
        rows = {row["query"]: row for row in read_tsv(out / "DE_NOVO_CANDIDATE_GATE.tsv")}
        self.assertEqual(rows["Query_A"]["gate_status"], "fail")
        audit = json.loads((out / "DE_NOVO_AGGREGATE_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "pass")
        self.assertEqual(audit["overall_gate_status"], "fail")

    def test_declared_timeout_preserves_all_rows_and_partial_evidence(self) -> None:
        for run in aggregate.EXPECTED_RUNS:
            self.fixture.write_run(
                run,
                technical_status=(
                    "technical_incomplete" if run == "DRR853909" else "pass"
                ),
            )
        out = self.root / "incomplete"
        self.assertEqual(self.fixture.invoke(out), 1)
        combined = read_tsv(out / "ALL_DE_NOVO_RECOVERY.tsv")
        self.assertEqual(len(combined), 18)
        preserved = [
            row for row in combined
            if row["run"] == "DRR853909" and row["query"] == "Query_A"
        ][0]
        self.assertTrue(preserved["best_single_contig"])
        self.assertEqual(preserved["gate_status"], "technical_incomplete")
        self.assertEqual(
            {row["gate_status"] for row in read_tsv(out / "DE_NOVO_CANDIDATE_GATE.tsv")},
            {"technical_incomplete"},
        )

    def test_aggregate_recomputes_maximum_internal_gap(self) -> None:
        for run in aggregate.EXPECTED_RUNS:
            self.fixture.write_run(
                run,
                high_gap=("DRR853907", "Query_A"),
            )
        out = self.root / "gap"
        self.assertEqual(self.fixture.invoke(out), 1)
        rows = {row["query"]: row for row in read_tsv(out / "DE_NOVO_CANDIDATE_GATE.tsv")}
        self.assertEqual(rows["Query_A"]["gate_status"], "fail")

    def test_missing_retained_full_assembly_is_technical_incomplete(self) -> None:
        self.fixture.write_six()
        retained = (
            self.fixture.collected / "DRR853907" / "de_novo"
            / "FULL_ASSEMBLY.fna.gz"
        )
        retained.unlink()
        out = self.root / "missing_retained_assembly"
        self.assertEqual(self.fixture.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AGGREGATE_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("missing_retained_assembly_evidence", audit["failures"][0])

    def test_blank_gap_and_impossible_counts_are_fail_closed(self) -> None:
        self.fixture.write_six()
        table = self.fixture.collected / "DRR853907" / "DE_NOVO_RECOVERY.tsv"
        rows = read_tsv(table)
        rows[0]["max_internal_query_gap_nt"] = ""
        rows[0]["covered_query_nt"] = "0"
        rows[0]["hsp_count"] = "0"
        denovo.write_tsv(table, denovo.RECOVERY_FIELDS, rows)
        out = self.root / "impossible_recovery_metrics"
        self.assertEqual(self.fixture.invoke(out), 1)
        audit = json.loads((out / "DE_NOVO_AGGREGATE_AUDIT.json").read_text())
        self.assertEqual(audit["technical_status"], "technical_incomplete")
        self.assertIn("missing_best_contig_metric", audit["failures"][0])


if __name__ == "__main__":
    unittest.main()
