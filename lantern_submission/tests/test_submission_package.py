#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOADER = ROOT / "lantern_submission/scripts/load_explicit_mapping.py"
TOY = ROOT / "lantern_submission/configs/toy_official_mapping.tsv"
ACTIVE = ROOT / "lantern_submission/configs/active_official_mapping.tsv"


def run_loader(mapping: Path, out: Path, individuals: int, timepoints: int, sample_ids: str):
    return subprocess.run([
        sys.executable, str(LOADER),
        "--mapping", str(mapping),
        "--expected-individuals", str(individuals),
        "--expected-timepoints", str(timepoints),
        "--expected-sample-ids", sample_ids,
        "--require-consecutive-timepoints",
        "--out", str(out),
    ], capture_output=True, text=True)


class MappingTests(unittest.TestCase):
    def test_toy_exact_mapping_and_determinism(self):
        expected = [["0", "1"], ["2", "7"], ["3", "5"], ["4", "8"], ["6", "19"],
                    ["9", "12"], ["10", "15"], ["11", "13"], ["14", "16"], ["17", "18"]]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            first = root / "a"
            second = root / "b"
            ids = ",".join(str(i) for i in range(20))
            self.assertEqual(run_loader(TOY, first, 10, 2, ids).returncode, 0)
            self.assertEqual(run_loader(TOY, second, 10, 2, ids).returncode, 0)
            a = json.loads((first / "MAPPING_FREEZE.json").read_text())
            self.assertEqual([group["samples"] for group in a["groups"]], expected)
            self.assertFalse(a["sequential_pairing_used"])
            self.assertFalse(a["similarity_pairing_used"])
            self.assertEqual(
                (first / "MAPPING_FREEZE.json").read_bytes(),
                (second / "MAPPING_FREEZE.json").read_bytes(),
            )

    def test_active_exact_mapping(self):
        expected = [["0", "2", "5", "19"], ["1", "4", "11", "30"], ["3", "8", "20", "24"],
                    ["6", "7", "9", "31"], ["10", "28", "29", "32"], ["12", "14", "18", "34"],
                    ["13", "22", "25", "26"], ["15", "16", "27", "33"], ["17", "21", "23", "35"]]
        with tempfile.TemporaryDirectory() as td:
            result = run_loader(ACTIVE, Path(td), 9, 4, ",".join(str(i) for i in range(36)))
            self.assertEqual(result.returncode, 0, result.stderr)
            obj = json.loads((Path(td) / "MAPPING_FREEZE.json").read_text())
            self.assertEqual([group["samples"] for group in obj["groups"]], expected)

    def test_duplicate_sample_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mapping = root / "bad.tsv"
            mapping.write_text(
                "individual_id\tsample_id\ttimepoint\n1\t0\t1\n1\t0\t2\n",
                encoding="utf-8",
            )
            result = run_loader(mapping, root / "out", 1, 2, "0")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("duplicate sample_id", result.stderr)

    def test_duplicate_individual_timepoint_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mapping = root / "bad.tsv"
            mapping.write_text(
                "individual_id\tsample_id\ttimepoint\n1\t0\t1\n1\t1\t1\n",
                encoding="utf-8",
            )
            result = run_loader(mapping, root / "out", 1, 2, "0,1")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("duplicate individual_id/timepoint", result.stderr)

    def test_missing_and_extra_sample_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mapping = root / "bad.tsv"
            mapping.write_text(
                "individual_id\tsample_id\ttimepoint\n1\t0\t1\n1\t2\t2\n",
                encoding="utf-8",
            )
            result = run_loader(mapping, root / "out", 1, 2, "0,1")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("sample set mismatch", result.stderr)

    def test_nonconsecutive_timepoint_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mapping = root / "bad.tsv"
            mapping.write_text(
                "individual_id\tsample_id\ttimepoint\n1\t0\t1\n1\t1\t3\n",
                encoding="utf-8",
            )
            result = run_loader(mapping, root / "out", 1, 2, "0,1")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("expected [1, 2]", result.stderr)


if __name__ == "__main__":
    unittest.main()
