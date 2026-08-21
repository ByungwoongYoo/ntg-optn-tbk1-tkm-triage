import csv
import hashlib
import json
import unittest
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).parents[1]
EXPECTED_FIELDS = (
    "accession", "context_group", "expected_title", "expected_length",
    "sequence_sha256", "expected_queries", "distinct_rank",
)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if header is not None:
                records[header.split()[0]] = (header, "".join(chunks))
            header = raw[1:]
            chunks = []
        else:
            chunks.append(raw.strip())
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks))
    return records


class CurrentNrPanelContractTests(unittest.TestCase):
    def test_exact_six_panel_matches_fasta_and_curated_manifest(self):
        panel = read_tsv(ROOT / "current_nr_top_hit_proteins.tsv")
        curated = read_tsv(ROOT / "panax_rdrp_curated_references.tsv")
        fasta = read_fasta(ROOT / "current_nr_top_hit_proteins.faa")
        panel_by_accession = {row["accession"]: row for row in panel}
        curated_current = {
            row["accession"]: row for row in curated
            if row["role"] == "current_nr_top_hit_context"
        }
        self.assertEqual(len(panel), 6)
        self.assertEqual(tuple(panel[0]), EXPECTED_FIELDS)
        self.assertEqual(len(panel_by_accession), 6)
        self.assertEqual(set(panel_by_accession), set(curated_current))
        self.assertEqual(set(panel_by_accession), set(fasta))
        ranks = defaultdict(list)
        for accession, row in panel_by_accession.items():
            header, sequence = fasta[accession]
            curated_row = curated_current[accession]
            self.assertEqual(row["expected_title"], header.split(maxsplit=1)[1])
            self.assertEqual(len(sequence), int(row["expected_length"]))
            self.assertEqual(
                hashlib.sha256(sequence.encode()).hexdigest(),
                row["sequence_sha256"],
            )
            self.assertEqual(row["context_group"], curated_row["context_group"])
            self.assertEqual(row["sequence_sha256"], curated_row["expected_sequence_sha256"])
            self.assertEqual(int(row["expected_length"]), int(curated_row["expected_aa_length"]))
            for query in row["expected_queries"].split(";"):
                ranks[query].append(int(row["distinct_rank"]))
            self.assertIn(
                row["expected_queries"],
                {"PNX_Picorna_A1;PNX_Picorna_A2", "PNX_Picorna_B"},
            )
        self.assertEqual(
            {query: sorted(values) for query, values in ranks.items()},
            {
                "PNX_Picorna_A1": [1, 2],
                "PNX_Picorna_A2": [1, 2],
                "PNX_Picorna_B": [1, 2, 3, 4],
            },
        )

    def test_panel_checksum_manifest_is_current(self):
        expected = {}
        for raw in (ROOT / "current_nr_top_hit_proteins.sha256").read_text().splitlines():
            if not raw or raw.startswith("#"):
                continue
            digest, name = raw.split(maxsplit=1)
            expected[name] = digest
        self.assertEqual(
            expected,
            {
                "current_nr_top_hit_proteins.faa": hashlib.sha256(
                    (ROOT / "current_nr_top_hit_proteins.faa").read_bytes()
                ).hexdigest(),
                "current_nr_top_hit_proteins.tsv": hashlib.sha256(
                    (ROOT / "current_nr_top_hit_proteins.tsv").read_bytes()
                ).hexdigest(),
            },
        )

    def test_remote_partition_control_files_are_hash_bound(self):
        manifest = json.loads((ROOT / "remote_partition_controls.json").read_text())
        controls = manifest["controls"]
        self.assertEqual(len({row["control"] for row in controls}), len(controls))
        for row in controls:
            path = ROOT / row["fasta_file"]
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                row["fasta_file_sha256"],
            )
            records = read_fasta(path)
            self.assertEqual(set(records), {row["control"]})
            sequence = records[row["control"]][1]
            self.assertEqual(len(sequence), int(row["length"]))
            self.assertEqual(
                hashlib.sha256(sequence.encode()).hexdigest(),
                row["sequence_sha256"],
            )
        by_mode = {
            mode: {row["control"] for row in controls if mode in row["required_modes"]}
            for mode in ("nt_nonviral", "nt_panax")
        }
        self.assertEqual(
            by_mode["nt_nonviral"] - by_mode["nt_panax"],
            {"PNX_NonPanax_mtDNA_control"},
        )


if __name__ == "__main__":
    unittest.main()
