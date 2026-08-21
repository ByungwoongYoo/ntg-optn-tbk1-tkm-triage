import hashlib
import importlib.util
import json
import shlex
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "finalize_panax_sequence_gate.py"
SPEC = importlib.util.spec_from_file_location("panax_sequence_finalizer", MODULE_PATH)
FINALIZER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(FINALIZER)


def fasta_rows(payload):
    records = []
    name = None
    chunks = []
    for raw in payload.decode().splitlines():
        if raw.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = raw[1:].split()[0]
            chunks = []
        else:
            chunks.append(raw.strip())
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return [
        {
            "id": name,
            "length": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
        for name, sequence in records
    ]


class FinalizerRemoteContractTests(unittest.TestCase):
    def build_nt_panax_fixture(self, root):
        collected = root / "collected"
        preflight = collected / "panax-query-preflight"
        remote = collected / "panax-remote-nt_panax"
        preflight.mkdir(parents=True)
        remote.mkdir(parents=True)

        candidate_fasta = "".join(
            f">{candidate}\n{'ACGT' * (index + 2)}\n"
            for index, candidate in enumerate(FINALIZER.CANDIDATES)
        ).encode()
        (preflight / "panax_three_contigs.fna").write_bytes(candidate_fasta)
        source = MODULE_PATH.parent
        query_bytes = b"".join([
            candidate_fasta,
            (source / "remote_partition_controls.fna").read_bytes(),
        ])
        (remote / "SEARCH_QUERIES.fna").write_bytes(query_bytes)

        mode = "nt_panax"
        argv = FINALIZER.expected_remote_argv(mode)
        (remote / "COMMAND.txt").write_text(shlex.join(argv) + "\n")
        query_sha = hashlib.sha256(query_bytes).hexdigest()
        query_argument = "remote-nt_panax/SEARCH_QUERIES.fna"
        (remote / "QUERY_SHA256.txt").write_text(
            f"{query_sha}  {query_argument}\n"
        )
        rows = fasta_rows(query_bytes)
        specs = FINALIZER.MODE_CONTROL_SPECS[mode]
        control_ids = sorted(specs)
        expected = {
            "query_file": query_argument,
            "query_file_sha256": query_sha,
            "candidate_ids": sorted(FINALIZER.CANDIDATES),
            "validation_control_ids": control_ids,
            "validation_controls": [
                {"id": control, **specs[control]} for control in control_ids
            ],
            "queries": rows,
        }
        (remote / "EXPECTED_QUERIES.json").write_text(
            json.dumps(expected, indent=2) + "\n"
        )
        status = {
            "mode": mode,
            "database": "nt",
            "query_file": query_argument,
            "query_sha256": query_sha,
            "query_count": len(rows),
            "query_ids": [row["id"] for row in rows],
            "validation_control_ids": control_ids,
            "expected_query_lengths": {
                row["id"]: row["length"] for row in rows
            },
            "validation_control_results": {
                control: {
                    "expected_accession": spec.get("expected_accession"),
                    "min_query_coverage": spec.get("min_query_coverage"),
                    "min_identity": spec.get("min_identity"),
                    "validated_accessions": [spec["expected_accession"]],
                    "validated": True,
                }
                for control, spec in specs.items()
            },
        }
        return collected, remote, status

    def test_exact_command_expected_and_preflight_fasta_contract_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            collected, _, status = self.build_nt_panax_fixture(Path(tmp))
            self.assertEqual(
                FINALIZER.validate_remote_contract(
                    collected, "nt_panax", status
                ),
                [],
            )

    def test_entrez_deletion_and_replacement_are_rejected(self):
        for mutation in ("delete", "replace"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                collected, remote, status = self.build_nt_panax_fixture(Path(tmp))
                argv = FINALIZER.expected_remote_argv("nt_panax")
                index = argv.index("-entrez_query")
                if mutation == "delete":
                    del argv[index:index + 2]
                else:
                    argv[index + 1] = "all[filter] NOT txid10239[ORGN]"
                (remote / "COMMAND.txt").write_text(shlex.join(argv) + "\n")
                failures = FINALIZER.validate_remote_contract(
                    collected, "nt_panax", status
                )
                self.assertIn(
                    "remote_command_contract_mismatch:nt_panax", failures
                )

    def test_expected_sequence_hash_and_id_mutations_are_rejected(self):
        for mutation in ("hash", "id"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                collected, remote, status = self.build_nt_panax_fixture(Path(tmp))
                path = remote / "EXPECTED_QUERIES.json"
                expected = json.loads(path.read_text())
                if mutation == "hash":
                    expected["queries"][0]["sequence_sha256"] = "0" * 64
                else:
                    expected["queries"][0]["id"] = "PNX_Picorna_MUTATED"
                path.write_text(json.dumps(expected, indent=2) + "\n")
                failures = FINALIZER.validate_remote_contract(
                    collected, "nt_panax", status
                )
                self.assertIn(
                    "remote_expected_query_contract_mismatch:nt_panax",
                    failures,
                )

    def test_expected_control_spec_mutation_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            collected, remote, status = self.build_nt_panax_fixture(Path(tmp))
            path = remote / "EXPECTED_QUERIES.json"
            expected = json.loads(path.read_text())
            expected["validation_controls"][0]["min_identity"] = 50.0
            path.write_text(json.dumps(expected, indent=2) + "\n")
            failures = FINALIZER.validate_remote_contract(
                collected, "nt_panax", status
            )
            self.assertIn(
                "remote_expected_query_contract_mismatch:nt_panax", failures
            )

    def test_status_and_preflight_fasta_mutations_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            collected, _, status = self.build_nt_panax_fixture(Path(tmp))
            status["query_sha256"] = "0" * 64
            failures = FINALIZER.validate_remote_contract(
                collected, "nt_panax", status
            )
            self.assertIn(
                "remote_status_contract_mismatch:nt_panax:query_sha256",
                failures,
            )

        with tempfile.TemporaryDirectory() as tmp:
            collected, _, status = self.build_nt_panax_fixture(Path(tmp))
            preflight_query = (
                collected / "panax-query-preflight" / "panax_three_contigs.fna"
            )
            preflight_query.write_text(
                preflight_query.read_text().replace("ACGT", "TGCA", 1)
            )
            failures = FINALIZER.validate_remote_contract(
                collected, "nt_panax", status
            )
            self.assertIn("remote_constructed_query_mismatch:nt_panax", failures)
            self.assertIn(
                "remote_expected_query_contract_mismatch:nt_panax", failures
            )


if __name__ == "__main__":
    unittest.main()
