import unittest

from virus_hunt.finalization.finalize_panax_sequence_gate import (
    completed_remote_near_identical,
)


class CompletedRemoteEvidenceTests(unittest.TestCase):
    def test_incomplete_modes_cannot_drive_near_identical_flags(self):
        status = {
            "protein_viral": {
                "per_query": {
                    "candidate": {"near_identical_qcov80_pident90_count": 3}
                }
            },
            "nt_viral": {
                "per_query": {
                    "candidate": {"near_identical_qcov80_pident95_count": 2}
                }
            },
        }
        observed = completed_remote_near_identical(
            "candidate", status, {"protein_viral": False, "nt_viral": False}
        )
        self.assertEqual(observed, (False, False, False))

    def test_only_complete_modes_contribute(self):
        status = {
            "protein_viral": {
                "per_query": {
                    "candidate": {"near_identical_qcov80_pident90_count": 1}
                }
            },
            "nt_viral": {
                "per_query": {
                    "candidate": {"near_identical_qcov80_pident95_count": 1}
                }
            },
        }
        observed = completed_remote_near_identical(
            "candidate", status, {"protein_viral": True, "nt_viral": False}
        )
        self.assertEqual(observed, (True, False, False))
        self.assertEqual(
            completed_remote_near_identical(
                "candidate", status, {"protein_viral": True, "nt_viral": True}
            ),
            (True, True, True),
        )


if __name__ == "__main__":
    unittest.main()
