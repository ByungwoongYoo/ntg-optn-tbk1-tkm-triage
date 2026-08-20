#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "baseline_preserving_augment.py"


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, sequence in records:
            f.write(f">{name}\n{sequence}\n")


def run(cmd):
    subprocess.run(cmd, check=True)


def main():
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        b1 = "A" * 1200
        b2 = "C" * 1300
        temporal = "G" * 1400
        long_only = "T" * 1500
        consensus = "ACGT" * 400
        redundant = "A" * 1190 + "G" * 10
        exact = b2
        write_fasta(d / "backbone.fa", [("B1", b1), ("B2", b2)])
        write_fasta(
            d / "candidates.fa",
            [("TEMP", temporal), ("LONG", long_only), ("CONS", consensus), ("RED", redundant), ("EXACT", exact)],
        )
        with open(d / "meta.tsv", "w", newline="") as f:
            fields = ["representative_id", "assembler_count", "source_count", "scopes", "assemblers"]
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(
                [
                    {"representative_id": "TEMP", "assembler_count": 1, "source_count": 1, "scopes": "longitudinal", "assemblers": "flye"},
                    {"representative_id": "LONG", "assembler_count": 1, "source_count": 1, "scopes": "longitudinal", "assemblers": "flye"},
                    {"representative_id": "CONS", "assembler_count": 2, "source_count": 2, "scopes": "longitudinal", "assemblers": "flye,megahit"},
                    {"representative_id": "RED", "assembler_count": 2, "source_count": 2, "scopes": "longitudinal", "assemblers": "flye,megahit"},
                    {"representative_id": "EXACT", "assembler_count": 2, "source_count": 2, "scopes": "longitudinal", "assemblers": "flye,megahit"},
                ]
            )
        fields = ["contig_id", "sample_id", "short_breadth", "short_depth", "long_breadth", "long_depth", "long_spanning_reads"]
        with open(d / "evidence.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for cid in ["TEMP", "LONG", "CONS", "RED", "EXACT"]:
                for sample in ["t0", "t1"]:
                    row = {"contig_id": cid, "sample_id": sample, "short_breadth": 0.8, "short_depth": 2, "long_breadth": 0, "long_depth": 0, "long_spanning_reads": 0}
                    if cid == "LONG":
                        row.update({"short_breadth": 0.2, "short_depth": 0.2, "long_breadth": 0.9, "long_depth": 2, "long_spanning_reads": 2 if sample == "t0" else 0})
                    w.writerow(row)
        with open(d / "c2b.paf", "w") as f:
            f.write("RED\t1200\t0\t1200\t+\tB1\t1200\t0\t1200\t1190\t1200\t60\n")
        cfg = {
            "version": "test",
            "minimum_contig_length": 1000,
            "maximum_n_fraction": 0.05,
            "support_breadth": 0.5,
            "support_depth": 0.75,
            "minimum_unique_long_spans": 2,
            "minimum_rescue_timepoints": 2,
            "minimum_assembler_sources_for_consensus": 2,
            "weights": {"length": 0.12, "source_consensus": 0.23, "temporal_recurrence": 0.25, "short_breadth": 0.12, "long_breadth": 0.10, "long_spanning": 0.18},
            "selection_score_minimum": 0.28,
            "cluster_identity": 0.97,
            "cluster_shorter_coverage": 0.85,
        }
        (d / "cfg.json").write_text(json.dumps(cfg))

        def execute(ablation):
            o = d / ablation
            run(
                [
                    "python",
                    str(SCRIPT),
                    "--backbone",
                    str(d / "backbone.fa"),
                    "--candidates",
                    str(d / "candidates.fa"),
                    "--metadata",
                    str(d / "meta.tsv"),
                    "--evidence",
                    str(d / "evidence.tsv"),
                    "--candidate-to-backbone-paf",
                    str(d / "c2b.paf"),
                    "--config",
                    str(d / "cfg.json"),
                    "--out",
                    str(o),
                    "--ablation",
                    ablation,
                ]
            )
            rows = {r["contig_id"]: r for r in csv.DictReader(open(o / "AUGMENTATION_DECISIONS.tsv"), delimiter="\t")}
            records = []
            name = None
            parts = []
            for line in (o / "LANTERN_BACKBONE_AUGMENTED.fasta").read_text().splitlines():
                if line.startswith(">"):
                    if name is not None:
                        records.append((name, "".join(parts)))
                    name = line[1:]
                    parts = []
                else:
                    parts.append(line)
            if name is not None:
                records.append((name, "".join(parts)))
            assert records[:2] == [("B1", b1), ("B2", b2)]
            return rows, json.loads((o / "AUGMENTATION_SUMMARY.json").read_text())

        full, summary = execute("full")
        assert full["TEMP"]["selected"] == "true"
        assert full["LONG"]["selected"] == "true"
        assert full["CONS"]["selected"] == "true"
        assert full["RED"]["selected"] == "false" and full["RED"]["reason"] == "alignment_redundant_to_backbone"
        assert full["EXACT"]["selected"] == "false" and full["EXACT"]["reason"] == "exact_backbone_or_prior_duplicate"
        assert summary["backbone_sequences_preserved_exactly"]

        no_longitudinal, _ = execute("no_longitudinal")
        assert no_longitudinal["TEMP"]["selected"] == "false"
        no_long, _ = execute("no_long")
        assert no_long["LONG"]["selected"] == "false"
        no_consensus, _ = execute("no_consensus")
        assert no_consensus["CONS"]["consensus_gate"] == "false"

    print("BASELINE_PRESERVING_AUGMENTATION_TEST_PASS")


if __name__ == "__main__":
    main()
