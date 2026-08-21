import hashlib
import json
import shlex
import tempfile
import unittest
from pathlib import Path

from virus_hunt.finalization import aggregate_panax_nonviral_splits as AGG
from virus_hunt.finalization import finalize_panax_sequence_gate as FINALIZER


def write_checksums(root):
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append(f"{digest}  ./{path.relative_to(root).as_posix()}\n")
    (root / "SHA256SUMS.txt").write_text("".join(rows))


def rebind_root_archive_metadata(remote, status, candidate, payload):
    """Keep outer metadata coherent while replacing one final child ASN."""
    archive = remote / "SPLITS" / candidate / "RESULTS.asn"
    archive.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    status["split_results"][candidate]["result_archive_sha256"] = digest
    (remote / "SEARCH_STATUS.json").write_text(
        json.dumps(status, indent=2) + "\n"
    )
    manifest_path = remote / "SPLIT_REQUESTS.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["splits"] = status["split_results"]
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    (remote / "RESULT_ARCHIVES.sha256").write_text("".join(
        f"{status['split_results'][item]['result_archive_sha256']}  "
        f"SPLITS/{item}/RESULTS.asn\n"
        for item in AGG.CANDIDATES
    ))


def sync_child_status_metadata(remote, status, candidate, child):
    child_path = remote / "SPLITS" / candidate / "SEARCH_STATUS.json"
    child_path.write_text(json.dumps(child, indent=2) + "\n")
    status["split_results"][candidate]["search_status_sha256"] = hashlib.sha256(
        child_path.read_bytes()
    ).hexdigest()
    (remote / "SEARCH_STATUS.json").write_text(
        json.dumps(status, indent=2) + "\n"
    )
    manifest_path = remote / "SPLIT_REQUESTS.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["splits"] = status["split_results"]
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


class ProteinNonviralSplitAggregationTests(unittest.TestCase):
    query_prefix = "remote-protein_nonviral"

    def build_fixture(self, root, *, omit=None):
        candidate_fasta = root / "candidates.faa"
        control_fasta = root / "control.faa"
        candidate_fasta.write_text(
            ">PNX_Picorna_A1 source=mock\nAAAAAAAAAAAA\n"
            ">PNX_Picorna_A2 source=mock\nCCCCCCCCCCCC\n"
            ">PNX_Picorna_B source=mock\nDDDDDDDDDDDD\n"
        )
        control_fasta.write_bytes(
            (Path(AGG.__file__).parent / "remote_partition_controls.faa").read_bytes()
        )
        out = root / "out"
        out.mkdir()
        (out / "SPLITS").mkdir()
        (out / "CANDIDATE_QUERIES.faa").write_bytes(candidate_fasta.read_bytes())
        (out / "remote_partition_controls.faa").write_bytes(control_fasta.read_bytes())
        (out / "SEARCH_QUERIES.faa").write_bytes(
            candidate_fasta.read_bytes() + control_fasta.read_bytes()
        )
        candidates = AGG._fasta_records(candidate_fasta)
        control = AGG._fasta_records(control_fasta)[AGG.CONTROL]
        for candidate in AGG.CANDIDATES:
            if candidate == omit:
                continue
            split = out / "SPLITS" / candidate
            split.mkdir()
            query_bytes = candidates[candidate]["payload"] + control["payload"]
            query_arg = (
                f"{self.query_prefix}/SPLITS/{candidate}/SEARCH_QUERIES.faa"
            )
            query_sha = hashlib.sha256(query_bytes).hexdigest()
            candidate_spec = AGG._query_spec(candidate, candidates[candidate])
            control_spec = AGG._query_spec(AGG.CONTROL, control)
            lengths = {
                candidate: candidate_spec["length"],
                AGG.CONTROL: control_spec["length"],
            }
            (split / "SEARCH_QUERIES.faa").write_bytes(query_bytes)
            (split / "RESULTS.asn").write_bytes(
                f"archive:{candidate}\n".encode()
            )
            (split / "ATTEMPT_ARCHIVES").mkdir()
            (split / "ATTEMPT_ARCHIVES" / "attempt1.asn").write_bytes(
                (split / "RESULTS.asn").read_bytes()
            )
            expected = {
                "query_file": query_arg,
                "query_file_sha256": query_sha,
                "candidate_ids": [candidate],
                "validation_control_ids": [AGG.CONTROL],
                "validation_controls": [{"id": AGG.CONTROL, **AGG.CONTROL_SPEC}],
                "queries": [candidate_spec, control_spec],
                "split_candidate_id": candidate,
            }
            (split / "EXPECTED_QUERIES.json").write_text(
                json.dumps(expected, indent=2) + "\n"
            )
            (split / "COMMAND.txt").write_text(
                shlex.join(AGG._expected_argv(query_arg)) + "\n"
            )
            (split / "QUERY_SHA256.txt").write_text(
                f"{query_sha}  {query_arg}\n"
            )
            control_sequence = control["sequence"]
            control_length = str(len(control_sequence))
            control_values = [
                AGG.CONTROL, "YP_009121238.1", "YP_009121238.1",
                "ref|YP_009121238.1|", "100", control_length,
                control_length, control_length, "1", control_length, "1",
                control_length, "0", "500", "100", "44586",
                "Panax ginseng", "mock L2", control_sequence, control_sequence,
            ]
            hit_rows = [dict(zip(AGG.HIT_FIELDS, control_values))]
            reports = []
            for query_id, query_len in lengths.items():
                hits = []
                if query_id == AGG.CONTROL:
                    hits = [{
                        "description": [{
                            "id": "ref|YP_009121238.1|",
                            "accession": "YP_009121238.1",
                            "title": "mock L2", "taxid": 44586,
                        }],
                        "len": query_len,
                        "hsps": [{
                            "query_from": 1, "query_to": query_len,
                            "hit_from": 1, "hit_to": query_len,
                            "align_len": query_len, "identity": query_len,
                            "evalue": 0.0, "bit_score": 500.0,
                            "qseq": control_sequence,
                            "hseq": control_sequence,
                        }],
                    }]
                search = {
                    "query_title": query_id,
                    "query_len": query_len,
                    "hits": hits,
                    "stat": {
                        "db_num": 1000, "db_len": 500000,
                        "kappa": 0.041, "lambda": 0.267,
                        "entropy": 0.14,
                    },
                }
                if not hits:
                    search["message"] = "No hits found"
                reports.append({"report": {
                    "program": "blastp",
                    "version": "BLASTP 2.17.0+",
                    "reference": "mock BLAST reference",
                    "search_target": {"db": "nr"},
                    "params": {
                        "matrix": "BLOSUM62", "expect": 1e-5,
                        "gap_open": 11, "gap_extend": 1,
                        "filter": "L;", "cbs": 2,
                    },
                    "results": {"search": search},
                }})
            (split / "RESULTS.json").write_text(
                json.dumps({"BlastOutput2": reports})
            )
            (split / "HITS.tsv").write_text("\t".join(control_values) + "\n")
            per_query = AGG._summarize_hits(
                hit_rows, lengths, [], f"fixture:{candidate}"
            )
            control_results = AGG._recompute_control_results(
                hit_rows, {AGG.CONTROL: AGG.CONTROL_SPEC}
            )
            attempt = {
                "attempt": "1", "start_utc": "2026-08-21T00:00:00Z",
                "end_utc": "2026-08-21T00:00:01Z",
                "backoff_before_seconds": "0", "attempt_timeout_seconds": "60",
                "blast_rc": "0", "json_formatter_rc": "0",
                "tsv_formatter_rc": "0", "validator_rc": "0",
                "failure_stage": "none", "failure_class": "none",
                "retryable": "0",
                "result_archive_bytes": str((split / "RESULTS.asn").stat().st_size),
                "result_archive_sha256": hashlib.sha256(
                    (split / "RESULTS.asn").read_bytes()
                ).hexdigest(),
            }
            (split / "REMOTE_ATTEMPTS.tsv").write_text(
                "\t".join(attempt) + "\n" + "\t".join(attempt.values()) + "\n"
            )
            (split / "ATTEMPT_COUNT.txt").write_text("1\n")
            (split / "SUCCESS_ATTEMPT.txt").write_text("1\n")
            (split / "SEARCH_SUCCESS.txt").write_text("1\n")
            (split / "TERMINATION_REASON.txt").write_text("success\n")
            status = {
                "mode": "protein_nonviral", "database": "nr",
                "query_file": query_arg, "query_sha256": query_sha,
                "query_count": 2, "query_ids": [candidate, AGG.CONTROL],
                "validation_control_ids": [AGG.CONTROL],
                "expected_query_lengths": lengths,
                "request_strategy": "protein_nonviral_candidate_control_split_v1",
                "split_candidate_id": candidate,
                "technical_complete": True,
                "command_completed_successfully": True,
                "result_archive_valid": True,
                "attempt_count": 1,
                "attempt_history": [attempt],
                "termination_reason": "success",
                "per_query": per_query,
                "validation_control_results": control_results,
            }
            (split / "SEARCH_STATUS.json").write_text(
                json.dumps(status, indent=2) + "\n"
            )
            write_checksums(split)
        return out, candidate_fasta, control_fasta

    def aggregate(self, out, candidate_fasta, control_fasta):
        return AGG.aggregate(
            out, candidate_fasta, control_fasta, self.query_prefix
        )

    def test_success_aggregates_three_candidates_and_deduplicates_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, candidates, control = self.build_fixture(Path(tmp))
            status = self.aggregate(out, candidates, control)
            self.assertTrue(status["technical_complete"])
            self.assertEqual(set(status["split_results"]), set(AGG.CANDIDATES))
            self.assertEqual(status["per_query"][AGG.CONTROL]["hit_count"], 1)
            hit_lines = (out / "HITS.tsv").read_text().splitlines()
            self.assertEqual(len(hit_lines), 1)
            self.assertTrue(hit_lines[0].startswith(AGG.CONTROL + "\t"))

    def test_success_satisfies_independent_finalizer_split_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            self.assertEqual(
                FINALIZER.validate_protein_nonviral_split_contract(
                    collected, status
                ),
                [],
            )

    def test_missing_and_invalid_split_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, candidates, control = self.build_fixture(
                Path(tmp), omit="PNX_Picorna_B"
            )
            status = self.aggregate(out, candidates, control)
            self.assertFalse(status["technical_complete"])
            self.assertIn("split_missing:PNX_Picorna_B", status["split_validation_failures"])
            self.assertEqual((out / "HITS.tsv").read_text(), "")

        with tempfile.TemporaryDirectory() as tmp:
            out, candidates, control = self.build_fixture(Path(tmp))
            split = out / "SPLITS" / "PNX_Picorna_A2"
            path = split / "EXPECTED_QUERIES.json"
            payload = json.loads(path.read_text())
            payload["queries"][0]["sequence_sha256"] = "0" * 64
            path.write_text(json.dumps(payload, indent=2) + "\n")
            write_checksums(split)
            status = self.aggregate(out, candidates, control)
            self.assertFalse(status["technical_complete"])
            self.assertIn(
                "split_expected_queries_mismatch:PNX_Picorna_A2",
                status["split_validation_failures"],
            )

    def test_control_failure_in_any_split_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, candidates, control = self.build_fixture(Path(tmp))
            split = out / "SPLITS" / "PNX_Picorna_B"
            path = split / "SEARCH_STATUS.json"
            payload = json.loads(path.read_text())
            payload["validation_control_results"][AGG.CONTROL]["validated"] = False
            path.write_text(json.dumps(payload, indent=2) + "\n")
            write_checksums(split)
            status = self.aggregate(out, candidates, control)
            self.assertFalse(status["technical_complete"])
            self.assertIn(
                "split_control_status_failed:PNX_Picorna_B:PNX_Panax_L2_control",
                status["split_validation_failures"],
            )

    def test_cross_split_database_signature_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, candidates, control = self.build_fixture(Path(tmp))
            split = out / "SPLITS" / "PNX_Picorna_A1"
            path = split / "RESULTS.json"
            payload = json.loads(path.read_text())
            for report in payload["BlastOutput2"]:
                report["report"]["results"]["search"]["stat"]["db_len"] += 1
            path.write_text(json.dumps(payload))
            write_checksums(split)
            status = self.aggregate(out, candidates, control)
            self.assertFalse(status["technical_complete"])
            self.assertIn(
                "cross_split_database_or_request_signature_mismatch",
                status["split_validation_failures"],
            )

    def test_final_archive_replacement_and_success_marker_tamper_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            split = out / "SPLITS" / "PNX_Picorna_A1"
            (split / "RESULTS.asn").write_bytes(b"replacement archive\n")
            write_checksums(split)
            failed = self.aggregate(out, candidates, control)
            self.assertFalse(failed["technical_complete"])
            self.assertIn(
                "split_success_archive_mismatch:PNX_Picorna_A1",
                failed["split_validation_failures"],
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            rebind_root_archive_metadata(
                remote, status, "PNX_Picorna_A1", b"replacement archive\n"
            )
            failures = FINALIZER.validate_protein_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_success_archive_mismatch:protein_nonviral:"
                "PNX_Picorna_A1",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            split = remote / "SPLITS" / "PNX_Picorna_A1"
            (split / "SUCCESS_ATTEMPT.txt").write_text("0\n")
            failures = FINALIZER.validate_protein_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_success_marker_mismatch:protein_nonviral:"
                "PNX_Picorna_A1",
                failures,
            )

    def test_finalizer_recomputes_raw_summaries_controls_and_signature(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            candidate = "PNX_Picorna_A1"
            child_path = remote / "SPLITS" / candidate / "SEARCH_STATUS.json"
            child = json.loads(child_path.read_text())
            forged = {
                "query_length": 12, "hit_count": 1,
                "near_identical_qcov80_pident90_count": 1,
                "near_identical_qcov80_pident95_count": 1,
                "top_hit": {"saccver": "FORGED.1", "bitscore": "999"},
            }
            child["per_query"][candidate] = forged
            status["per_query"][candidate] = forged
            sync_child_status_metadata(remote, status, candidate, child)
            failures = FINALIZER.validate_protein_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_per_query_summary_mismatch:PNX_Picorna_A1",
                failures,
            )
            self.assertIn(
                "remote_split_candidate_raw_aggregation_mismatch:"
                "PNX_Picorna_A1",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            (remote / "SPLITS" / "PNX_Picorna_A1" / "HITS.tsv").write_text("")
            failures = FINALIZER.validate_protein_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_archive_structural_failure:protein_nonviral:"
                "PNX_Picorna_A1",
                failures,
            )
            self.assertIn(
                "remote_split_control_summary_mismatch:PNX_Picorna_A1",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, control = self.build_fixture(root)
            status = self.aggregate(out, candidates, control)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_partial_orfs.faa").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-protein_nonviral"
            out.rename(remote)
            candidate = "PNX_Picorna_A1"
            status["split_results"][candidate]["database_signature_sha256"] = (
                "f" * 64
            )
            child = json.loads(
                (remote / "SPLITS" / candidate / "SEARCH_STATUS.json").read_text()
            )
            sync_child_status_metadata(remote, status, candidate, child)
            failures = FINALIZER.validate_protein_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_signature_invalid:PNX_Picorna_A1",
                failures,
            )


class NtNonviralSplitAggregationTests(unittest.TestCase):
    query_prefix = "remote-nt_nonviral"

    def build_fixture(self, root):
        source = Path(AGG.__file__).parent
        candidate_fasta = root / "candidates.fna"
        candidate_fasta.write_text(
            ">PNX_Picorna_A1 source=mock\nACGTACGTACGT\n"
            ">PNX_Picorna_A2 source=mock\nCCCCAAAATTTT\n"
            ">PNX_Picorna_B source=mock\nGGGGTTTTAAAA\n"
        )
        control_fastas = [
            source / "remote_partition_controls.fna",
            source / "remote_nonpanax_control.fna",
        ]
        out = root / "out"
        out.mkdir()
        (out / "SPLITS").mkdir()
        (out / "CANDIDATE_QUERIES.fna").write_bytes(candidate_fasta.read_bytes())
        for control_path in control_fastas:
            (out / control_path.name).write_bytes(control_path.read_bytes())
        (out / "SEARCH_QUERIES.fna").write_bytes(
            candidate_fasta.read_bytes()
            + b"".join(path.read_bytes() for path in control_fastas)
        )
        candidates = AGG._fasta_records(candidate_fasta)
        controls = {}
        for path in control_fastas:
            controls.update(AGG._fasta_records(path))
        config = AGG.MODE_CONFIGS["nt_nonviral"]
        control_ids = sorted(config["controls"])
        for candidate in AGG.CANDIDATES:
            split = out / "SPLITS" / candidate
            split.mkdir()
            query_bytes = candidates[candidate]["payload"] + b"".join(
                controls[control_id]["payload"] for control_id in controls
            )
            query_arg = f"{self.query_prefix}/SPLITS/{candidate}/SEARCH_QUERIES.fna"
            query_sha = hashlib.sha256(query_bytes).hexdigest()
            query_specs = [AGG._query_spec(candidate, candidates[candidate])] + [
                AGG._query_spec(control_id, controls[control_id])
                for control_id in controls
            ]
            lengths = {row["id"]: row["length"] for row in query_specs}
            (split / "SEARCH_QUERIES.fna").write_bytes(query_bytes)
            (split / "RESULTS.asn").write_bytes(f"nt:{candidate}\n".encode())
            (split / "ATTEMPT_ARCHIVES").mkdir()
            (split / "ATTEMPT_ARCHIVES" / "attempt1.asn").write_bytes(
                (split / "RESULTS.asn").read_bytes()
            )
            expected = {
                "query_file": query_arg,
                "query_file_sha256": query_sha,
                "candidate_ids": [candidate],
                "validation_control_ids": control_ids,
                "validation_controls": [
                    {"id": control_id, **config["controls"][control_id]}
                    for control_id in control_ids
                ],
                "queries": query_specs,
                "split_candidate_id": candidate,
            }
            (split / "EXPECTED_QUERIES.json").write_text(
                json.dumps(expected, indent=2) + "\n"
            )
            (split / "COMMAND.txt").write_text(
                shlex.join(AGG._expected_argv(query_arg, "nt_nonviral")) + "\n"
            )
            (split / "QUERY_SHA256.txt").write_text(
                f"{query_sha}  {query_arg}\n"
            )
            hit_lines = []
            hit_rows = []
            json_hits = {}
            for control_id in control_ids:
                sequence = controls[control_id]["sequence"]
                length = str(len(sequence))
                accession = config["controls"][control_id]["expected_accession"]
                values = [
                    control_id, accession, accession, f"ref|{accession}|",
                    "100", length, length, length, "1", length, "1", length,
                    "0", "500", "100", "1", "mock", "mock control",
                    sequence, sequence,
                ]
                hit_lines.append("\t".join(values))
                hit_rows.append(dict(zip(AGG.HIT_FIELDS, values)))
                json_hits[control_id] = [{
                    "description": [{
                        "id": f"ref|{accession}|", "accession": accession,
                        "title": "mock control", "taxid": 1,
                    }],
                    "len": len(sequence),
                    "hsps": [{
                        "query_from": 1, "query_to": len(sequence),
                        "hit_from": 1, "hit_to": len(sequence),
                        "align_len": len(sequence), "identity": len(sequence),
                        "evalue": 0.0, "bit_score": 500.0,
                        "qseq": sequence, "hseq": sequence,
                    }],
                }]
            reports = []
            for query_id, query_len in lengths.items():
                hits = json_hits.get(query_id, [])
                search = {
                    "query_title": query_id, "query_len": query_len,
                    "hits": hits,
                    "stat": {
                        "db_num": 2000, "db_len": 900000,
                        "kappa": 0.7, "lambda": 1.2, "entropy": 1.3,
                    },
                }
                if not hits:
                    search["message"] = "No hits found"
                reports.append({"report": {
                    "program": "blastn", "version": "BLASTN 2.17.0+",
                    "reference": "mock BLAST reference",
                    "search_target": {"db": "nt"},
                    "params": {
                        "expect": 1e-5, "sc_match": 2, "sc_mismatch": -3,
                        "gap_open": 5, "gap_extend": 2, "filter": "L;m;",
                    },
                    "results": {"search": search},
                }})
            (split / "RESULTS.json").write_text(
                json.dumps({"BlastOutput2": reports})
            )
            (split / "HITS.tsv").write_text("\n".join(hit_lines) + "\n")
            per_query = AGG._summarize_hits(
                hit_rows, lengths, [], f"fixture:{candidate}"
            )
            control_results = AGG._recompute_control_results(
                hit_rows, config["controls"]
            )
            attempt_header = [
                "attempt", "start_utc", "end_utc", "backoff_before_seconds",
                "attempt_timeout_seconds", "blast_rc", "json_formatter_rc",
                "tsv_formatter_rc", "validator_rc", "failure_stage",
                "failure_class", "retryable", "result_archive_bytes",
                "result_archive_sha256",
            ]
            attempt_values = [
                "1", "2026-08-21T00:00:00Z", "2026-08-21T00:00:01Z",
                "0", "60", "0", "0", "0", "0", "none", "none", "0",
                str((split / "RESULTS.asn").stat().st_size),
                hashlib.sha256((split / "RESULTS.asn").read_bytes()).hexdigest(),
            ]
            (split / "REMOTE_ATTEMPTS.tsv").write_text(
                "\t".join(attempt_header) + "\n" + "\t".join(attempt_values) + "\n"
            )
            attempt = dict(zip(attempt_header, attempt_values))
            (split / "ATTEMPT_COUNT.txt").write_text("1\n")
            (split / "SUCCESS_ATTEMPT.txt").write_text("1\n")
            (split / "SEARCH_SUCCESS.txt").write_text("1\n")
            (split / "TERMINATION_REASON.txt").write_text("success\n")
            status = {
                "mode": "nt_nonviral", "database": "nt",
                "query_file": query_arg, "query_sha256": query_sha,
                "query_count": 3,
                "query_ids": [candidate, *list(controls)],
                "validation_control_ids": control_ids,
                "expected_query_lengths": lengths,
                "request_strategy": config["child_strategy"],
                "split_candidate_id": candidate,
                "technical_complete": True,
                "command_completed_successfully": True,
                "result_archive_valid": True,
                "attempt_count": 1,
                "attempt_history": [attempt],
                "termination_reason": "success",
                "per_query": per_query,
                "validation_control_results": control_results,
            }
            (split / "SEARCH_STATUS.json").write_text(
                json.dumps(status, indent=2) + "\n"
            )
            write_checksums(split)
        return out, candidate_fasta, control_fastas

    def run_aggregate(self, out, candidates, controls):
        return AGG.aggregate(
            out, candidates, controls, self.query_prefix, "nt_nonviral"
        )

    def test_success_requires_both_controls_and_satisfies_finalizer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            self.assertTrue(status["technical_complete"])
            self.assertEqual(
                set(status["validation_control_results"]),
                set(AGG.MODE_CONFIGS["nt_nonviral"]["controls"]),
            )
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            out.rename(collected / "panax-remote-nt_nonviral")
            self.assertEqual(
                FINALIZER.validate_nt_nonviral_split_contract(collected, status),
                [],
            )

    def test_one_of_two_controls_failing_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            split = out / "SPLITS" / "PNX_Picorna_A2"
            path = split / "SEARCH_STATUS.json"
            payload = json.loads(path.read_text())
            payload["validation_control_results"][
                "PNX_NonPanax_mtDNA_control"
            ]["validated"] = False
            path.write_text(json.dumps(payload, indent=2) + "\n")
            write_checksums(split)
            status = self.run_aggregate(out, candidates, controls)
            self.assertFalse(status["technical_complete"])
            self.assertIn(
                "split_control_status_failed:PNX_Picorna_A2:"
                "PNX_NonPanax_mtDNA_control",
                status["split_validation_failures"],
            )

    def test_finalizer_rejects_cross_candidate_child_hit_injection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-nt_nonviral"
            out.rename(remote)
            split_hits = remote / "SPLITS" / "PNX_Picorna_A1" / "HITS.tsv"
            injected = split_hits.read_text().splitlines()[0].split("\t")
            injected[0] = "PNX_Picorna_A2"
            split_hits.write_text(
                split_hits.read_text() + "\t".join(injected) + "\n"
            )
            aggregate_hits = remote / "HITS.tsv"
            aggregate_hits.write_text(
                aggregate_hits.read_text() + "\t".join(injected) + "\n"
            )
            failures = FINALIZER.validate_nt_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_hit_query_mismatch:nt_nonviral:"
                "PNX_Picorna_A1:PNX_Picorna_A2",
                failures,
            )
            self.assertIn(
                "remote_split_hit_aggregation_mismatch:nt_nonviral",
                failures,
            )

    def test_nt_finalizer_recomputes_candidate_and_both_raw_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-nt_nonviral"
            out.rename(remote)
            candidate = "PNX_Picorna_A2"
            child_path = remote / "SPLITS" / candidate / "SEARCH_STATUS.json"
            child = json.loads(child_path.read_text())
            forged = {
                "query_length": 12, "hit_count": 1,
                "near_identical_qcov80_pident90_count": 1,
                "near_identical_qcov80_pident95_count": 1,
                "top_hit": {"saccver": "FORGED.1", "bitscore": "999"},
            }
            child["per_query"][candidate] = forged
            status["per_query"][candidate] = forged
            sync_child_status_metadata(remote, status, candidate, child)
            failures = FINALIZER.validate_nt_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_per_query_summary_mismatch:nt_nonviral:"
                "PNX_Picorna_A2",
                failures,
            )
            self.assertIn(
                "remote_split_candidate_raw_aggregation_mismatch:nt_nonviral:"
                "PNX_Picorna_A2",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-nt_nonviral"
            out.rename(remote)
            candidate = "PNX_Picorna_A2"
            hits = remote / "SPLITS" / candidate / "HITS.tsv"
            kept = [
                line for line in hits.read_text().splitlines()
                if not line.startswith("PNX_NonPanax_mtDNA_control\t")
            ]
            hits.write_text("\n".join(kept) + "\n")
            failures = FINALIZER.validate_nt_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_control_hit_failed:nt_nonviral:PNX_Picorna_A2:"
                "PNX_NonPanax_mtDNA_control",
                failures,
            )
            self.assertIn(
                "remote_split_archive_structural_failure:nt_nonviral:"
                "PNX_Picorna_A2",
                failures,
            )

    def test_empty_final_archive_and_success_row_tamper_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            split = out / "SPLITS" / "PNX_Picorna_A2"
            (split / "RESULTS.asn").write_bytes(b"")
            write_checksums(split)
            failed = self.run_aggregate(out, candidates, controls)
            self.assertFalse(failed["technical_complete"])
            self.assertIn(
                "split_archive_missing:PNX_Picorna_A2",
                failed["split_validation_failures"],
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-nt_nonviral"
            out.rename(remote)
            rebind_root_archive_metadata(remote, status, "PNX_Picorna_A2", b"")
            failures = FINALIZER.validate_nt_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_final_archive_invalid:nt_nonviral:"
                "PNX_Picorna_A2",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, candidates, controls = self.build_fixture(root)
            status = self.run_aggregate(out, candidates, controls)
            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / "panax_three_contigs.fna").write_bytes(
                candidates.read_bytes()
            )
            remote = collected / "panax-remote-nt_nonviral"
            out.rename(remote)
            attempts = (
                remote / "SPLITS" / "PNX_Picorna_A2" / "REMOTE_ATTEMPTS.tsv"
            )
            lines = attempts.read_text().splitlines()
            values = lines[1].split("\t")
            values[-1] = "0" * 64
            attempts.write_text(lines[0] + "\n" + "\t".join(values) + "\n")
            failures = FINALIZER.validate_nt_nonviral_split_contract(
                collected, status
            )
            self.assertIn(
                "remote_split_success_archive_mismatch:nt_nonviral:"
                "PNX_Picorna_A2",
                failures,
            )


if __name__ == "__main__":
    unittest.main()
