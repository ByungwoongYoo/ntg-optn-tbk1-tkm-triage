import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
RUNNER = ROOT / "run_panax_remote_search.sh"
CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")


FAKE_BLASTP = r'''#!/usr/bin/env python3
import pathlib,sys
if "-version" in sys.argv:
    print("blastp: 2.17.0+")
    raise SystemExit(0)
out = pathlib.Path(sys.argv[sys.argv.index("-out") + 1])
out.write_bytes(b"mock ASN archive\n")
'''


FAKE_FORMATTER = r'''#!/usr/bin/env python3
import json,os,pathlib,sys
if "-version" in sys.argv:
    print("blast_formatter: 2.17.0+")
    raise SystemExit(0)
out = pathlib.Path(sys.argv[sys.argv.index("-out") + 1])
outfmt = sys.argv[sys.argv.index("-outfmt") + 1]
if outfmt == "15":
    reports=[]
    for qid in ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B"):
        reports.append({"report": {
            "program": "blastp",
            "version": "BLASTP 2.17.0+",
            "reference": "mock BLAST reference",
            "search_target": {"db": "tsa_nr"},
            "params": {
                "matrix": "BLOSUM62", "expect": 1e-5,
                "gap_open": 11, "gap_extend": 1, "filter": "L;", "cbs": 2,
            },
            "results": {"search": {
                "query_title": qid, "query_len": 12,
                "hits": [], "message": "No hits found",
                "stat": {"db_num": 1000, "db_len": 500000,
                         "kappa": 0.041, "lambda": 0.267, "entropy": 0.14},
            }},
        }})
    out.write_text(json.dumps({"BlastOutput2": reports}))
else:
    fields = [
        "PNX_Picorna_A1", "REP00001.1", "REP00001", "ref|REP00001.1|",
        "100", "12", "12", "12", "1", "12", "1", "12", "0", "500",
        "100", "1", "N/A", "mock protein", "AAAAAAAAAAAA", "AAAAAAAAAAAA",
    ]
    case=os.environ["FAKE_TSV_CASE"]
    if case == "empty":
        out.write_text("")
        raise SystemExit(0)
    elif case == "nineteen_fields":
        fields.pop()
    elif case == "nonnumeric_bitscore":
        fields[13]="not-a-number"
    out.write_text("\t".join(fields) + "\n")
'''


class RemoteStatusArtifactTests(unittest.TestCase):
    def test_valid_no_hit_archive_completes_and_is_checksummed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            queries = root / "queries"
            out = root / "out"
            tools.mkdir()
            queries.mkdir()
            for name, content in (("blastp", FAKE_BLASTP), ("blast_formatter", FAKE_FORMATTER)):
                path = tools / name
                path.write_text(content)
                path.chmod(0o755)
            (queries / "panax_three_partial_orfs.faa").write_text(
                "".join(f">{name}\nAAAAAAAAAAAA\n" for name in CANDIDATES)
            )
            env = dict(os.environ)
            env.update({
                "PATH": f"{tools}:{env['PATH']}",
                "PANAX_REMOTE_MAX_ATTEMPTS": "1",
                "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS": "60",
                "PANAX_REMOTE_SEARCH_BUDGET_SECONDS": "300",
                "FAKE_TSV_CASE": "empty",
            })
            completed = subprocess.run(
                ["bash", str(RUNNER), "protein_tsa", str(queries), str(out)],
                env=env, text=True, capture_output=True, timeout=30,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            status = json.loads((out / "SEARCH_STATUS.json").read_text())
            self.assertTrue(status["technical_complete"])
            self.assertTrue(status["command_completed_successfully"])
            self.assertEqual(status["result_row_count"], 0)
            verified = subprocess.run(
                ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=out,
                text=True, capture_output=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)

    def run_malformed_case(self, case):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            queries = root / "queries"
            out = root / "out"
            tools.mkdir()
            queries.mkdir()
            blastp = tools / "blastp"
            formatter = tools / "blast_formatter"
            blastp.write_text(FAKE_BLASTP)
            formatter.write_text(FAKE_FORMATTER)
            blastp.chmod(0o755)
            formatter.chmod(0o755)
            (queries / "panax_three_partial_orfs.faa").write_text(
                "".join(f">{name}\nAAAAAAAAAAAA\n" for name in CANDIDATES)
            )
            env = dict(os.environ)
            env.update({
                "PATH": f"{tools}:{env['PATH']}",
                "PANAX_REMOTE_MAX_ATTEMPTS": "1",
                "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS": "60",
                "PANAX_REMOTE_SEARCH_BUDGET_SECONDS": "300",
                "FAKE_TSV_CASE": case,
            })
            completed = subprocess.run(
                ["bash", str(RUNNER), "protein_tsa", str(queries), str(out)],
                env=env, text=True, capture_output=True, timeout=30,
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            status = json.loads((out / "SEARCH_STATUS.json").read_text())
            self.assertFalse(status["technical_complete"])
            self.assertTrue(status["tsv_validation_errors"])
            self.assertEqual(status["result_row_count"], 0)
            self.assertTrue(all(
                row["hit_count"] == 0 and row["top_hit"] is None
                for row in status["per_query"].values()
            ))
            self.assertIn("no hit or control annotation", status["annotation_validation"])
            self.assertTrue((out / "FINISHED_UTC.txt").is_file())
            self.assertTrue((out / "SHA256SUMS.txt").is_file())
            verified = subprocess.run(
                ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=out,
                text=True, capture_output=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)

    def test_malformed_field_count_still_finalizes_failure_artifact(self):
        self.run_malformed_case("nineteen_fields")

    def test_nonnumeric_row_still_finalizes_failure_artifact(self):
        self.run_malformed_case("nonnumeric_bitscore")

    def test_missing_formatter_still_finalizes_preflight_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            queries = root / "queries"
            out = root / "out"
            tools.mkdir()
            queries.mkdir()
            blastp = tools / "blastp"
            blastp.write_text(FAKE_BLASTP)
            blastp.chmod(0o755)
            formatter = tools / "blast_formatter"
            formatter.write_text("#!/usr/bin/env bash\nexit 127\n")
            formatter.chmod(0o755)
            (queries / "panax_three_partial_orfs.faa").write_text(
                "".join(f">{name}\nAAAAAAAAAAAA\n" for name in CANDIDATES)
            )
            env = dict(os.environ)
            env["PATH"] = f"{tools}:{env['PATH']}"
            completed = subprocess.run(
                ["bash", str(RUNNER), "protein_tsa", str(queries), str(out)],
                env=env, text=True, capture_output=True, timeout=30,
            )
            self.assertNotEqual(completed.returncode, 0)
            status = json.loads((out / "SEARCH_STATUS.json").read_text())
            self.assertFalse(status["technical_complete"])
            self.assertEqual(status["failure_stage"], "preflight_or_internal")
            self.assertTrue((out / "FINISHED_UTC.txt").is_file())
            verified = subprocess.run(
                ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=out,
                text=True, capture_output=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)

    def test_overflow_sized_timeout_still_finalizes_preflight_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            queries = root / "queries"
            out = root / "out"
            tools.mkdir()
            queries.mkdir()
            for name, content in (("blastp", FAKE_BLASTP), ("blast_formatter", FAKE_FORMATTER)):
                path = tools / name
                path.write_text(content)
                path.chmod(0o755)
            (queries / "panax_three_partial_orfs.faa").write_text(
                "".join(f">{name}\nAAAAAAAAAAAA\n" for name in CANDIDATES)
            )
            env = dict(os.environ)
            env.update({
                "PATH": f"{tools}:{env['PATH']}",
                "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS": "18446744073709552216",
                "FAKE_TSV_CASE": "nineteen_fields",
            })
            completed = subprocess.run(
                ["bash", str(RUNNER), "protein_tsa", str(queries), str(out)],
                env=env, text=True, capture_output=True, timeout=30,
            )
            self.assertEqual(completed.returncode, 2)
            status = json.loads((out / "SEARCH_STATUS.json").read_text())
            self.assertFalse(status["technical_complete"])
            self.assertEqual(status["runner_exit_code"], 2)
            self.assertTrue((out / "SHA256SUMS.txt").is_file())


if __name__ == "__main__":
    unittest.main()
