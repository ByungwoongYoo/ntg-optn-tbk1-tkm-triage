import tempfile
import unittest
from pathlib import Path

from virus_hunt.finalization.summarize_panax_fragment_depth import (
    augment,
    longest_zero_run,
    read_depth,
    read_fasta_lengths,
    summarize,
)


class FragmentDepthTests(unittest.TestCase):
    def test_summary_and_internal_zero_run(self):
        result = summarize([0, 2, 0, 0, 5, 0], "proper")
        self.assertEqual(result["proper_breadth_1x"], "0.333333")
        self.assertEqual(result["proper_breadth_5x"], "0.166667")
        self.assertEqual(result["proper_mean_depth"], "1.166667")
        self.assertEqual(result["proper_max_zero_run"], 2)
        self.assertEqual(result["proper_max_internal_zero_run"], 2)
        self.assertEqual(longest_zero_run([]), 0)

    def test_full_augmentation_and_coordinate_validation(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            fasta = root / "refs.fna"
            fasta.write_text(">A\nAAAA\n>B note\nCCCCC\n")
            depth_dir = root / "depth"
            depth_dir.mkdir()
            expected = {"A": [1, 2, 3, 4], "B": [0, 0, 5, 5, 0]}
            for reference, values in expected.items():
                for label in ("preduplicate", "nonduplicate"):
                    (depth_dir / f"{reference}.{label}.depth.tsv").write_text(
                        "".join(
                            f"{reference}\t{position}\t{value}\n"
                            for position, value in enumerate(values, 1)
                        )
                    )
            rows = [
                {"run": "R", "reference": "A", "count": "1"},
                {"run": "R", "reference": "B", "count": "2"},
            ]
            output = augment(rows, read_fasta_lengths(fasta), depth_dir)
            self.assertEqual(output[0]["proper_nonduplicate_breadth_1x"], "1.000000")
            self.assertEqual(output[1]["proper_preduplicate_max_internal_zero_run"], 0)

            bad = depth_dir / "A.preduplicate.depth.tsv"
            bad.write_text("A\t2\t1\n")
            with self.assertRaisesRegex(ValueError, "non-contiguous"):
                read_depth(bad, "A", 4)

    def test_rejects_reference_set_mismatch(self):
        with self.assertRaisesRegex(ValueError, "set mismatch"):
            augment([{"reference": "A"}], {"A": 4, "B": 5}, Path("unused"))


if __name__ == "__main__":
    unittest.main()
