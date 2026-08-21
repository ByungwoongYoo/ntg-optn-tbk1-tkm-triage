import contextlib
import csv
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from virus_hunt.finalization import aggregate_panax_fragment_audit as aggregate


METRIC_FIELDS = [
    "run", "alias", "reference", "reference_length_nt", "minimum_mate_mapq",
    "proper_fragments_mapq30_preduplicate",
    "proper_fragments_mapq30_nonduplicate", "duplicate_fragment_fraction",
    "panel_unique_kmer_fragments_nonduplicate", "panel_unique_kmer_fragment_fraction",
    "panel_unique_31mers_available", "panel_unique_31mers_observed",
    "panel_unique_31mer_fraction_observed", "distinct_fragment_endpoints",
    "fr_orientation_fraction", "template_span_median", "template_span_p10",
    "template_span_p90", "softclipped_fragment_fraction",
    "max_softclipped_bases_per_fragment",
    "proper_preduplicate_breadth_1x", "proper_preduplicate_breadth_5x",
    "proper_preduplicate_breadth_10x", "proper_preduplicate_mean_depth",
    "proper_preduplicate_median_depth", "proper_preduplicate_max_zero_run",
    "proper_preduplicate_max_internal_zero_run",
    "proper_nonduplicate_breadth_1x", "proper_nonduplicate_breadth_5x",
    "proper_nonduplicate_breadth_10x", "proper_nonduplicate_mean_depth",
    "proper_nonduplicate_median_depth", "proper_nonduplicate_max_zero_run",
    "proper_nonduplicate_max_internal_zero_run",
]
HASH_FIELDS = [
    "run", "alias", "reference", "minimum_mate_mapq",
    "pair_sequence_sha256", "endpoint_sha256", "fragment_fingerprint_sha256",
]
CONTINUITY_FIELDS = [
    "run", "alias", "reference", "reference_length", "minimum_mate_mapq",
    "boundary_after_nt", "direct_read_anchor_nt_each_side",
    "proper_fragments_mapq30_nonduplicate_spanning",
    "direct_read_spanning_fragments_anchor25",
    "direct_read_spanning_reads_anchor25",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def hashes(run: str, reference: str, count: int, shared: list[tuple[str, str]] | None = None):
    rows = []
    shared = shared or []
    for index in range(count):
        if index < len(shared):
            pair_hash, endpoint_hash = shared[index]
        else:
            pair_hash = hashlib.sha256(f"pair:{run}:{reference}:{index}".encode()).hexdigest()
            endpoint_hash = hashlib.sha256(f"end:{run}:{reference}:{index}".encode()).hexdigest()
        combined = hashlib.sha256(
            (reference + "\n" + pair_hash + "\n" + endpoint_hash).encode()
        ).hexdigest()
        rows.append({
            "run": run, "alias": aggregate.ALIASES[run], "reference": reference,
            "minimum_mate_mapq": 30, "pair_sequence_sha256": pair_hash,
            "endpoint_sha256": endpoint_hash,
            "fragment_fingerprint_sha256": combined,
        })
    return rows


class FragmentAggregateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / "collected"
        self.root.mkdir()
        self.out = Path(self.temp.name) / "out"
        self.out.mkdir()
        self.thresholds = Path(__file__).parents[1] / "virus_hunt/finalization/panax_enhanced_audit_thresholds.json"
        self.structure = self.thresholds.with_name("panax_candidate_structure.tsv")
        self.threshold_values = aggregate.load_thresholds(self.thresholds)
        self.lengths = aggregate.load_candidate_structure(
            self.structure, self.threshold_values
        )

    def tearDown(self):
        self.temp.cleanup()

    def build(self, count_by: dict[tuple[str, str], int] | None = None, shared=None):
        count_by = count_by or {}
        shared = shared or {}
        for run in aggregate.EXPECTED_RUNS:
            directory = self.root / run / "fragment_audit"
            directory.mkdir(parents=True)
            metric_rows = []
            hash_rows = []
            continuity_rows = []
            for reference in aggregate.REFERENCES:
                count = count_by.get((run, reference), 30)
                reference_length = self.lengths[reference]
                if count >= 25:
                    breadth_targets = (0.95, 0.60, 0.30)
                elif count >= 10:
                    breadth_targets = (0.85, 0.10, 0.02)
                elif count:
                    breadth_targets = (0.10, 0.0, 0.0)
                else:
                    breadth_targets = (0.0, 0.0, 0.0)
                covered1, covered5, covered10 = (
                    round(reference_length * target)
                    for target in breadth_targets
                )
                breadth1, breadth5, breadth10 = (
                    covered / reference_length
                    for covered in (covered1, covered5, covered10)
                )
                # This is the exact summary of a realizable synthetic depth
                # profile with consecutive 10x, 5x, 1x, and terminal 0x blocks.
                mean_depth = (
                    covered1 + 4 * covered5 + 5 * covered10
                ) / reference_length
                median_depth = 5 if covered5 * 2 > reference_length else (
                    1 if covered1 * 2 > reference_length else 0
                )
                max_zero_run = reference_length - covered1
                kmer_available = 100
                kmer_observed = min(count, kmer_available)
                metric_rows.append({
                    "run": run, "alias": aggregate.ALIASES[run], "reference": reference,
                    "reference_length_nt": reference_length,
                    "minimum_mate_mapq": 30,
                    "proper_fragments_mapq30_preduplicate": count,
                    "proper_fragments_mapq30_nonduplicate": count,
                    "duplicate_fragment_fraction": "0.000000" if count else "",
                    "panel_unique_kmer_fragments_nonduplicate": count,
                    "panel_unique_kmer_fragment_fraction": "1.000000" if count else "",
                    "panel_unique_31mers_available": kmer_available,
                    "panel_unique_31mers_observed": kmer_observed,
                    "panel_unique_31mer_fraction_observed": f"{kmer_observed / kmer_available:.6f}",
                    "distinct_fragment_endpoints": count,
                    "fr_orientation_fraction": "1.000000" if count else "",
                    "template_span_median": "290.000" if count else "",
                    "template_span_p10": 280 if count else "",
                    "template_span_p90": 300 if count else "",
                    "softclipped_fragment_fraction": "0.000000" if count else "",
                    "max_softclipped_bases_per_fragment": 0,
                    "proper_preduplicate_breadth_1x": f"{breadth1:.6f}",
                    "proper_preduplicate_breadth_5x": f"{breadth5:.6f}",
                    "proper_preduplicate_breadth_10x": f"{breadth10:.6f}",
                    "proper_nonduplicate_breadth_1x": f"{breadth1:.6f}",
                    "proper_nonduplicate_breadth_5x": f"{breadth5:.6f}",
                    "proper_nonduplicate_breadth_10x": f"{breadth10:.6f}",
                    "proper_preduplicate_mean_depth": f"{mean_depth:.6f}",
                    "proper_preduplicate_median_depth": f"{median_depth:.3f}",
                    "proper_preduplicate_max_zero_run": max_zero_run,
                    "proper_preduplicate_max_internal_zero_run": 0,
                    "proper_nonduplicate_mean_depth": f"{mean_depth:.6f}",
                    "proper_nonduplicate_median_depth": f"{median_depth:.3f}",
                    "proper_nonduplicate_max_zero_run": max_zero_run,
                    "proper_nonduplicate_max_internal_zero_run": 0,
                })
                hash_rows.extend(hashes(run, reference, count, shared.get((run, reference))))
                for boundary in range(250, reference_length, 250):
                    continuity_rows.append({
                        "run": run, "alias": aggregate.ALIASES[run], "reference": reference,
                        "reference_length": reference_length, "minimum_mate_mapq": 30,
                        "boundary_after_nt": boundary,
                        "direct_read_anchor_nt_each_side": 25,
                        "proper_fragments_mapq30_nonduplicate_spanning": min(3, count),
                        "direct_read_spanning_fragments_anchor25": min(1, count),
                        "direct_read_spanning_reads_anchor25": min(1, count),
                    })
            write_tsv(directory / "FRAGMENT_METRICS.tsv", METRIC_FIELDS, metric_rows)
            write_tsv(directory / "FRAGMENT_HASHES.tsv", HASH_FIELDS, hash_rows)
            write_tsv(directory / "CONTINUITY_SCAN.tsv", CONTINUITY_FIELDS, continuity_rows)

    def test_complete_six_run_audit_passes(self):
        self.build()
        status = aggregate.run(self.root, self.thresholds, self.out)
        self.assertTrue(status["technical_complete"])
        self.assertEqual(status["fragment_support_gate"], "pass")
        with (self.out / "CANDIDATE_FRAGMENT_GATE.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row["independent_validated_count"] == "6" for row in rows))
        with (self.out / "CROSS_RUN_FRAGMENT_SHARING.tsv").open(newline="") as handle:
            sharing = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(sharing), 45)

    def test_predeclared_shadow_suspect_excludes_smaller_run(self):
        reference = aggregate.REFERENCES[0]
        small, large = aggregate.EXPECTED_RUNS[:2]
        shared_pairs = [
            (
                hashlib.sha256(f"shared-pair:{index}".encode()).hexdigest(),
                hashlib.sha256(f"shared-end:{index}".encode()).hexdigest(),
            )
            for index in range(8)
        ]
        counts = {(small, reference): 10, (large, reference): 200}
        shared = {(small, reference): shared_pairs, (large, reference): shared_pairs}
        self.build(counts, shared)
        aggregate.run(self.root, self.thresholds, self.out)
        with (self.out / "CROSS_RUN_FRAGMENT_SHARING.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        row = next(
            row for row in rows
            if row["reference"] == reference and {row["run_a"], row["run_b"]} == {small, large}
        )
        self.assertEqual(row["shadow_qc_status"], "shadow_suspect")
        with (self.out / "CANDIDATE_FRAGMENT_GATE.tsv").open(newline="") as handle:
            candidate = next(
                row for row in csv.DictReader(handle, delimiter="\t")
                if row["reference"] == reference
            )
        self.assertIn(small, candidate["shadow_suspect_runs"])
        self.assertNotIn(small, candidate["independent_validated_runs"])

    def test_missing_run_and_tampered_combined_hash_fail_closed(self):
        self.build()
        missing = self.root / aggregate.EXPECTED_RUNS[-1]
        for path in sorted(missing.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()
        missing.rmdir()
        with self.assertRaises(aggregate.AggregateError):
            aggregate.run(self.root, self.thresholds, self.out)

        self.temp.cleanup()
        self.setUp()
        self.build()
        table = next(self.root.rglob("FRAGMENT_HASHES.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["fragment_fingerprint_sha256"] = "0" * 64
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(aggregate.AggregateError, "combined_fragment_hash_mismatch"):
            aggregate.run(self.root, self.thresholds, self.out)

    def test_warning_does_not_require_suspect_abundance_ratio(self):
        reference = aggregate.REFERENCES[0]
        small, large = aggregate.EXPECTED_RUNS[:2]
        shared_pairs = [
            (
                hashlib.sha256(f"warning-pair:{index}".encode()).hexdigest(),
                hashlib.sha256(f"warning-end:{index}".encode()).hexdigest(),
            )
            for index in range(6)
        ]
        self.build(
            {(small, reference): 10, (large, reference): 15},
            {(small, reference): shared_pairs, (large, reference): shared_pairs},
        )
        aggregate.run(self.root, self.thresholds, self.out)
        with (self.out / "CROSS_RUN_FRAGMENT_SHARING.tsv").open(newline="") as handle:
            row = next(
                row for row in csv.DictReader(handle, delimiter="\t")
                if row["reference"] == reference
                and {row["run_a"], row["run_b"]} == {small, large}
            )
        self.assertEqual(row["shadow_qc_status"], "warning_shadow")
        self.assertEqual(row["smaller_nonduplicate_fragments"], "10")

    def test_endpoint_metric_must_match_unique_endpoint_hashes(self):
        self.build()
        table = next(self.root.rglob("FRAGMENT_METRICS.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["distinct_fragment_endpoints"] = "0"
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(aggregate.AggregateError, "endpoint_hash_count_mismatch"):
            aggregate.run(self.root, self.thresholds, self.out)

    def test_malformed_numeric_input_materializes_fail_closed_status(self):
        self.build()
        table = next(self.root.rglob("CONTINUITY_SCAN.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["reference_length"] = "620.5"
        write_tsv(table, fields, rows)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(aggregate.main([
                "--input-root", str(self.root),
                "--thresholds", str(self.thresholds),
                "--out", str(self.out),
            ]), 1)
        status = json.loads((self.out / "TECHNICAL_STATUS.json").read_text())
        self.assertEqual(status["technical_status"], "technical_incomplete")

    def test_insufficient_shadow_assessment_is_not_independent_support(self):
        target_reference = aggregate.REFERENCES[0]
        insufficient_run, strong_run = aggregate.EXPECTED_RUNS[:2]
        by_run = {}
        for run in aggregate.EXPECTED_RUNS:
            by_run[run] = {}
            for reference in aggregate.REFERENCES:
                validated = reference == target_reference and run in {
                    insufficient_run, strong_run
                }
                strong = reference == target_reference and run == strong_run
                by_run[run][reference] = {
                    "fragment_preduplicate_validated_positive": str(validated).lower(),
                    "fragment_nonduplicate_validated_positive": str(validated).lower(),
                    "fragment_validated_positive": str(validated).lower(),
                    "fragment_strong_positive": str(strong).lower(),
                }
        sharing = [{
            "smaller_run": insufficient_run,
            "reference": target_reference,
            "shadow_qc_status": "insufficient_smaller_run_fingerprints",
        }]
        run_rows, candidates = aggregate.build_status(
            by_run, sharing, self.threshold_values
        )
        run_row = next(
            row for row in run_rows
            if row["run"] == insufficient_run and row["reference"] == target_reference
        )
        self.assertEqual(run_row["shadow_assessment_evaluable"], "false")
        self.assertEqual(
            run_row["eligible_as_independent_validated_support"], "false"
        )
        candidate = next(
            row for row in candidates if row["reference"] == target_reference
        )
        self.assertEqual(candidate["fragment_support_gate"], "fail")
        self.assertEqual(candidate["independent_validated_count"], 1)
        self.assertIn(
            insufficient_run, candidate["shadow_assessment_insufficient_runs"]
        )

    def test_pinned_structure_and_metric_lengths_are_enforced(self):
        self.build()
        table = next(self.root.rglob("FRAGMENT_METRICS.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["reference_length_nt"] = str(int(rows[0]["reference_length_nt"]) + 1)
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(
            aggregate.AggregateError, "metric_pinned_length_mismatch"
        ):
            aggregate.run(self.root, self.thresholds, self.out)

        tampered = Path(self.temp.name) / "tampered_structure.tsv"
        tampered.write_text(self.structure.read_text() + "\n")
        with self.assertRaisesRegex(
            aggregate.AggregateError, "candidate_structure_sha256_mismatch"
        ):
            aggregate.run(
                self.root, self.thresholds, self.out, structure_path=tampered
            )

    def test_zero_fragments_require_zero_depth_and_observed_kmer_support(self):
        run, reference = aggregate.EXPECTED_RUNS[0], aggregate.REFERENCES[0]
        self.build({(run, reference): 0})
        table = next(
            path for path in self.root.rglob("FRAGMENT_METRICS.tsv")
            if run in str(path)
        )
        fields, rows = aggregate.read_tsv(table)
        row = next(row for row in rows if row["reference"] == reference)
        row["proper_nonduplicate_mean_depth"] = "1.000000"
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(
            aggregate.AggregateError, "zero_fragment_nonzero_depth_metric"
        ):
            aggregate.run(self.root, self.thresholds, self.out)
        row["proper_nonduplicate_mean_depth"] = "0.000000"
        row["panel_unique_31mers_observed"] = "1"
        row["panel_unique_31mer_fraction_observed"] = "0.010000"
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(
            aggregate.AggregateError, "diagnostic_fragment_kmer_mismatch"
        ):
            aggregate.run(self.root, self.thresholds, self.out)

    def test_depth_mean_lower_bound_allows_only_rounding_tolerance(self):
        self.build()
        table = next(self.root.rglob("FRAGMENT_METRICS.tsv"))
        _, rows = aggregate.read_tsv(table)
        base = rows[0]
        minimum_mean = (
            float(base["proper_preduplicate_breadth_1x"])
            + 4 * float(base["proper_preduplicate_breadth_5x"])
            + 5 * float(base["proper_preduplicate_breadth_10x"])
        )

        tolerated = dict(base)
        for duplicate in ("preduplicate", "nonduplicate"):
            tolerated[f"proper_{duplicate}_mean_depth"] = (
                f"{minimum_mean - 0.000005:.6f}"
            )
        aggregate.validate_metric_row(
            tolerated, tolerated["run"], self.threshold_values
        )

        invalid = dict(tolerated)
        invalid["proper_preduplicate_mean_depth"] = (
            f"{minimum_mean - 0.000008:.6f}"
        )
        with self.assertRaisesRegex(
            aggregate.AggregateError, "mean_depth_below_breadth_lower_bound"
        ):
            aggregate.validate_metric_row(
                invalid, invalid["run"], self.threshold_values
            )

    def test_breadth_and_zero_run_invariants(self):
        self.build()
        table = next(self.root.rglob("FRAGMENT_METRICS.tsv"))
        _, rows = aggregate.read_tsv(table)
        base = rows[0]
        reference_length = int(base["reference_length_nt"])

        partial_without_zero = dict(base)
        partial_without_zero["proper_preduplicate_max_zero_run"] = "0"
        with self.assertRaisesRegex(
            aggregate.AggregateError, "partial_breadth_requires_zero_run"
        ):
            aggregate.validate_metric_row(
                partial_without_zero,
                partial_without_zero["run"],
                self.threshold_values,
            )

        zero_breadth_short_run = dict(base)
        for depth in (1, 5, 10):
            zero_breadth_short_run[f"proper_preduplicate_breadth_{depth}x"] = "0.000000"
        zero_breadth_short_run["proper_preduplicate_mean_depth"] = "0.000000"
        zero_breadth_short_run["proper_preduplicate_median_depth"] = "0.000"
        zero_breadth_short_run["proper_preduplicate_max_zero_run"] = str(
            reference_length - 1
        )
        with self.assertRaisesRegex(
            aggregate.AggregateError, "zero_breadth_max_zero_run_mismatch"
        ):
            aggregate.validate_metric_row(
                zero_breadth_short_run,
                zero_breadth_short_run["run"],
                self.threshold_values,
            )

        complete_with_zero = dict(base)
        complete_with_zero["proper_preduplicate_breadth_1x"] = "1.000000"
        complete_with_zero["proper_preduplicate_mean_depth"] = "5.000000"
        complete_with_zero["proper_preduplicate_max_zero_run"] = "1"
        with self.assertRaisesRegex(
            aggregate.AggregateError, "complete_breadth_forbids_zero_run"
        ):
            aggregate.validate_metric_row(
                complete_with_zero,
                complete_with_zero["run"],
                self.threshold_values,
            )

        internal_above_maximum = dict(base)
        maximum = int(internal_above_maximum["proper_preduplicate_max_zero_run"])
        internal_above_maximum["proper_preduplicate_max_internal_zero_run"] = str(
            maximum + 1
        )
        with self.assertRaisesRegex(
            aggregate.AggregateError, "internal_zero_run_exceeds_maximum"
        ):
            aggregate.validate_metric_row(
                internal_above_maximum,
                internal_above_maximum["run"],
                self.threshold_values,
            )

    def test_diagnostic_fragments_and_observed_kmers_must_agree(self):
        self.build()
        table = next(self.root.rglob("FRAGMENT_METRICS.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["panel_unique_kmer_fragments_nonduplicate"] = "0"
        rows[0]["panel_unique_kmer_fragment_fraction"] = "0.000000"
        self.assertGreater(int(rows[0]["panel_unique_31mers_observed"]), 0)
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(
            aggregate.AggregateError, "diagnostic_fragment_kmer_mismatch"
        ):
            aggregate.run(self.root, self.thresholds, self.out)

    def test_direct_spanning_read_count_cannot_be_below_fragment_count(self):
        self.build()
        table = next(self.root.rglob("CONTINUITY_SCAN.tsv"))
        fields, rows = aggregate.read_tsv(table)
        rows[0]["direct_read_spanning_fragments_anchor25"] = "1"
        rows[0]["direct_read_spanning_reads_anchor25"] = "0"
        write_tsv(table, fields, rows)
        with self.assertRaisesRegex(
            aggregate.AggregateError, "continuity_subcount_mismatch"
        ):
            aggregate.run(self.root, self.thresholds, self.out)

    def test_duplicate_blank_and_width_mismatched_tsv_rows_fail(self):
        malformed = {
            "duplicate.tsv": "a\ta\n1\t2\n",
            "blank.tsv": "a\t\n1\t2\n",
            "extra.tsv": "a\tb\n1\t2\t3\n",
            "missing.tsv": "a\tb\n1\n",
        }
        for name, contents in malformed.items():
            with self.subTest(name=name):
                path = Path(self.temp.name) / name
                path.write_text(contents)
                with self.assertRaises(aggregate.AggregateError):
                    aggregate.read_tsv(path)


if __name__ == "__main__":
    unittest.main()
