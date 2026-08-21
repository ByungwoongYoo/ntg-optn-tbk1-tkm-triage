#!/usr/bin/env python3
"""Fail-closed QC for unrooted homologous PF00680-core phylogenies."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
CURRENT_PANEL_FIELDS = (
    "accession", "context_group", "expected_title", "expected_length",
    "sequence_sha256", "expected_queries", "distinct_rank",
)


def canonical_panel_sha(rows: list[dict[str, str]]) -> str:
    payload = [
        {field: row[field] for field in CURRENT_PANEL_FIELDS}
        for row in sorted(rows, key=lambda item: item["accession"])
    ]
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise SystemExit(f"duplicate alignment ID in {path}: {name}")
            records[name] = ""
        elif name is None:
            raise SystemExit(f"sequence before FASTA header in {path}")
        else:
            records[name] += line.upper()
    if not records:
        raise SystemExit(f"empty FASTA: {path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise SystemExit(f"unaligned FASTA supplied as alignment: {path} {sorted(lengths)}")
    return records


def read_unaligned_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            if not name or name in records:
                raise SystemExit(f"missing or duplicate FASTA ID in {path}: {name}")
            records[name] = ""
        elif name is None:
            raise SystemExit(f"sequence before FASTA header in {path}")
        else:
            records[name] += line.upper()
    if not records or any(not sequence for sequence in records.values()):
        raise SystemExit(f"empty FASTA or record in {path}")
    return records


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def parse(self) -> Node:
        root = self.node()
        self.skip_space()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self.skip_space()
        if self.i != len(self.text):
            raise ValueError(f"unparsed Newick suffix at {self.i}: {self.text[self.i:self.i+40]}")
        return root

    def skip_space(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def label(self) -> str:
        self.skip_space()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",():;":
            self.i += 1
        return self.text[start:self.i].strip().strip("'\"")

    def branch(self) -> float:
        self.skip_space()
        if self.i >= len(self.text) or self.text[self.i] != ":":
            return 0.0
        self.i += 1
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        return float(self.text[start:self.i].strip())

    def node(self) -> Node:
        self.skip_space()
        if self.i < len(self.text) and self.text[self.i] == "(":
            self.i += 1
            children = [self.node()]
            while True:
                self.skip_space()
                if self.i < len(self.text) and self.text[self.i] == ",":
                    self.i += 1; children.append(self.node()); continue
                if self.i >= len(self.text) or self.text[self.i] != ")":
                    raise ValueError(f"expected ')' at Newick position {self.i}")
                self.i += 1; break
            node = Node(name=self.label(), children=children)
            node.length = self.branch()
            for child in children:
                child.parent = node
            return node
        name = self.label()
        if not name:
            raise ValueError(f"empty leaf at Newick position {self.i}")
        return Node(name=name, length=self.branch())


def leaves(root: Node) -> dict[str, Node]:
    out: dict[str, Node] = {}
    stack = [root]
    while stack:
        node = stack.pop()
        if node.children:
            stack.extend(node.children)
        else:
            if node.name in out:
                raise SystemExit(f"duplicate tree tip: {node.name}")
            out[node.name] = node
    return out


def ancestors(node: Node) -> dict[Node, float]:
    out = {node: 0.0}
    distance = 0.0
    while node.parent is not None:
        distance += max(0.0, node.length)
        node = node.parent
        out[node] = distance
    return out


def distance(left: Node, right: Node) -> float:
    a = ancestors(left)
    distance_right = 0.0
    node = right
    while node not in a:
        if node.parent is None:
            raise SystemExit("tree is disconnected")
        distance_right += max(0.0, node.length)
        node = node.parent
    return a[node] + distance_right


def alignment_stats(path: Path, expected: set[str]) -> tuple[dict[str, str], dict[str, object]]:
    records = read_fasta(path)
    if set(records) != expected:
        raise SystemExit(f"alignment tip mismatch in {path}: missing={sorted(expected-set(records))}, unexpected={sorted(set(records)-expected)}")
    length = len(next(iter(records.values())))
    occupancies = {name: sum(char not in "-?X" for char in seq) / length for name, seq in records.items()}
    informative = sum(
        len({seq[column] for seq in records.values() if seq[column] not in "-?X"}) >= 2
        for column in range(length)
    )
    stats = {
        "alignment": path.name, "sequence_count": len(records), "alignment_length": length,
        "overall_occupancy": sum(occupancies.values()) / len(occupancies),
        "minimum_occupancy": min(occupancies.values()),
        "minimum_candidate_occupancy": min(occupancies[name] for name in CANDIDATES),
        "variable_or_informative_columns": informative,
    }
    return records, stats


def tree_stats(path: Path, expected: set[str], reference_context: dict[str, str]) -> tuple[dict[str, object], dict[str, str]]:
    root = NewickParser(path.read_text()).parse()
    tips = leaves(root)
    if set(tips) != expected:
        raise SystemExit(f"tree tip mismatch in {path}: missing={sorted(expected-set(tips))}, unexpected={sorted(set(tips)-expected)}")
    reference_names = sorted(expected - set(CANDIDATES))
    ref_terminal = [max(0.0, tips[name].length) for name in reference_names]
    median_ref_terminal = statistics.median(ref_terminal) if ref_terminal else 0.0
    nearest_context: dict[str, str] = {}
    nearest_rows = []
    for candidate in CANDIDATES:
        ranked = sorted((distance(tips[candidate], tips[reference]), reference) for reference in reference_names)
        nearest_distance, nearest = ranked[0]
        nearest_context[candidate] = reference_context[nearest]
        ratio = None if median_ref_terminal == 0 else max(0.0, tips[candidate].length) / median_ref_terminal
        nearest_rows.append({
            "tree": path.name, "candidate": candidate, "nearest_reference": nearest,
            "nearest_reference_context": reference_context[nearest], "patristic_distance": f"{nearest_distance:.8f}",
            "candidate_terminal_branch": f"{max(0.0,tips[candidate].length):.8f}",
            "median_reference_terminal_branch": f"{median_ref_terminal:.8f}",
            "terminal_branch_ratio": "" if ratio is None else f"{ratio:.6f}",
            "long_branch_warning": str(ratio is not None and ratio > 5).lower(),
        })
    pairwise = {}
    for i, left in enumerate(CANDIDATES):
        for right in CANDIDATES[i + 1:]:
            pairwise[f"{left}|{right}"] = distance(tips[left], tips[right])
    return {"nearest_rows": nearest_rows, "candidate_pairwise_patristic": pairwise}, nearest_context


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cores", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--current-panel-manifest", type=Path, required=True)
    parser.add_argument("--untrimmed-alignment", type=Path, required=True)
    parser.add_argument("--trimmed-alignment", type=Path, required=True)
    parser.add_argument("--untrimmed-tree", type=Path, required=True)
    parser.add_argument("--trimmed-tree", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    with args.current_panel_manifest.open(newline="") as handle:
        panel_reader = csv.DictReader(handle, delimiter="\t")
        current_panel = list(panel_reader)
    if tuple(panel_reader.fieldnames or ()) != CURRENT_PANEL_FIELDS:
        raise SystemExit(
            "expanded phylogeny current-nr panel schema mismatch: "
            f"{panel_reader.fieldnames}"
        )
    if any(
        not row.get(field, "").strip()
        for row in current_panel for field in CURRENT_PANEL_FIELDS
    ):
        raise SystemExit("expanded phylogeny current-nr panel has blank contract fields")
    reference_context = {row["accession"]: row["context_group"] for row in manifest}
    if len(reference_context) != len(manifest) or not reference_context:
        raise SystemExit("empty or duplicated reference manifest")
    if len(manifest) != 25:
        raise SystemExit(f"expanded phylogeny requires exactly 25 references: {len(manifest)}")
    if sum(row.get("role") == "rooting_sensitivity_reference" for row in manifest) != 2:
        raise SystemExit("expanded phylogeny requires exactly two rooting-sensitivity references")
    panel_by_accession = {row.get("accession", ""): row for row in current_panel}
    manifest_current = {
        row["accession"]: row for row in manifest
        if row.get("role") == "current_nr_top_hit_context"
    }
    if len(current_panel) != 6 or len(panel_by_accession) != 6 or set(panel_by_accession) != set(manifest_current):
        raise SystemExit(
            "expanded phylogeny current-nr contract mismatch: "
            f"panel={sorted(panel_by_accession)}, manifest={sorted(manifest_current)}"
        )
    for accession, panel_row in panel_by_accession.items():
        manifest_row = manifest_current[accession]
        if panel_row.get("context_group") != manifest_row.get("context_group"):
            raise SystemExit(f"current-nr context mismatch in phylogeny manifest: {accession}")
        if panel_row.get("sequence_sha256") != manifest_row.get("expected_sequence_sha256"):
            raise SystemExit(f"current-nr full-sequence hash mismatch in phylogeny manifest: {accession}")
        if int(panel_row.get("expected_length", 0)) != int(manifest_row.get("expected_aa_length", 0)):
            raise SystemExit(f"current-nr length mismatch in phylogeny manifest: {accession}")
    observed_ranks = {candidate: [] for candidate in CANDIDATES}
    for accession, panel_row in panel_by_accession.items():
        query_contract = panel_row["expected_queries"]
        if query_contract not in {
            "PNX_Picorna_A1;PNX_Picorna_A2", "PNX_Picorna_B",
        }:
            raise SystemExit(f"current-nr query contract mismatch: {accession}")
        try:
            rank = int(panel_row["distinct_rank"])
        except ValueError as exc:
            raise SystemExit(f"current-nr rank is not an integer: {accession}") from exc
        for query in query_contract.split(";"):
            observed_ranks[query].append(rank)
    expected_ranks = {
        "PNX_Picorna_A1": [1, 2],
        "PNX_Picorna_A2": [1, 2],
        "PNX_Picorna_B": [1, 2, 3, 4],
    }
    if {key: sorted(value) for key, value in observed_ranks.items()} != expected_ranks:
        raise SystemExit(f"current-nr rank contract mismatch: {observed_ranks}")
    expected = set(CANDIDATES) | set(reference_context)
    if len(expected) != 28:
        raise SystemExit(f"expanded phylogeny requires exactly 28 tips: {len(expected)}")

    raw_cores = read_unaligned_fasta(args.cores)
    if set(raw_cores) != expected:
        raise SystemExit(
            f"raw core tip mismatch: missing={sorted(expected-set(raw_cores))}, "
            f"unexpected={sorted(set(raw_cores)-expected)}"
        )

    untrimmed_records, untrimmed = alignment_stats(args.untrimmed_alignment, expected)
    trimmed_records, trimmed = alignment_stats(args.trimmed_alignment, expected)
    for sequence_id, core in raw_cores.items():
        if untrimmed_records[sequence_id].replace("-", "") != core:
            raise SystemExit(
                f"untrimmed alignment is not sequence-bound to raw core: {sequence_id}"
            )
    tip_order = sorted(expected)
    untrimmed_columns = [
        tuple(untrimmed_records[name][column] for name in tip_order)
        for column in range(len(next(iter(untrimmed_records.values()))))
    ]
    trimmed_columns = [
        tuple(trimmed_records[name][column] for name in tip_order)
        for column in range(len(next(iter(trimmed_records.values()))))
    ]
    cursor = 0
    for trimmed_column in trimmed_columns:
        while cursor < len(untrimmed_columns) and untrimmed_columns[cursor] != trimmed_column:
            cursor += 1
        if cursor == len(untrimmed_columns):
            raise SystemExit(
                "trimmed alignment is not a common ordered column subset of the untrimmed alignment"
            )
        cursor += 1
    tree_untrimmed, contexts_untrimmed = tree_stats(args.untrimmed_tree, expected, reference_context)
    tree_trimmed, contexts_trimmed = tree_stats(args.trimmed_tree, expected, reference_context)
    nearest_rows = tree_untrimmed["nearest_rows"] + tree_trimmed["nearest_rows"]
    with (args.out / "NEAREST_REFERENCE.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(nearest_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(nearest_rows)
    alignment_rows = []
    for row in (untrimmed, trimmed):
        alignment_rows.append({key: f"{value:.6f}" if isinstance(value, float) else value for key, value in row.items()})
    with (args.out / "ALIGNMENT_QC.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(alignment_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(alignment_rows)
    sensitivity = [{
        "candidate": candidate, "untrimmed_nearest_context": contexts_untrimmed[candidate],
        "trimmed_nearest_context": contexts_trimmed[candidate],
        "nearest_context_stable": str(contexts_untrimmed[candidate] == contexts_trimmed[candidate]).lower(),
        "interpretation": "contextual sensitivity only; the broad Picornavirales trees are unrooted",
    } for candidate in CANDIDATES]
    with (args.out / "TREE_SENSITIVITY.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sensitivity[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(sensitivity)

    failures = []
    if int(untrimmed["alignment_length"]) < 250: failures.append("untrimmed_alignment_too_short")
    if int(trimmed["alignment_length"]) < 120: failures.append("trimmed_alignment_too_short")
    if float(trimmed["overall_occupancy"]) < 0.50: failures.append("trimmed_overall_occupancy_below_0.50")
    if float(trimmed["minimum_candidate_occupancy"]) < 0.50: failures.append("trimmed_candidate_occupancy_below_0.50")
    if int(trimmed["variable_or_informative_columns"]) < 50: failures.append("too_few_variable_or_informative_columns")
    status = {
        "generated_utc": datetime.now(timezone.utc).isoformat(), "technical_complete": not failures,
        "failures": failures, "primary_tree_interpretation": "unrooted broad Picornavirales PF00680-core context",
        "reference_count": len(reference_context), "expected_tip_count": len(expected),
        "untrimmed_bound_to_raw_cores": True,
        "trimmed_is_ordered_column_subset": True,
        "current_nr_reference_count": len(panel_by_accession),
        "current_nr_panel_contract_sha256": canonical_panel_sha(current_panel),
        "current_nr_reference_contract": {
            accession: {
                "role": manifest_current[accession]["role"],
                "context_group": manifest_current[accession]["context_group"],
                "sequence_sha256": panel_by_accession[accession]["sequence_sha256"],
                "expected_length": int(panel_by_accession[accession]["expected_length"]),
                "PF00680_core_sha256": hashlib.sha256(
                    raw_cores[accession].encode()
                ).hexdigest(),
                "present_in_both_alignments_and_trees": True,
            }
            for accession in sorted(panel_by_accession)
        },
        "rooting_note": "Picornaviridae references are retained as sensitivity/context references, not asserted as an external outgroup to Picornavirales.",
        "untrimmed_alignment": untrimmed, "trimmed_alignment": trimmed,
        "untrimmed_tree": tree_untrimmed, "trimmed_tree": tree_trimmed,
        "nearest_context_stability": sensitivity,
        "claim_boundary": "Phylogenetic placement is contextual and cannot by itself delimit a species or establish host, replication, or disease association.",
    }
    (args.out / "TREE_QC.json").write_text(json.dumps(status, indent=2) + "\n")
    report = [
        "# Panax A1/A2/B homologous RdRP-core phylogeny", "",
        f"Technical gate: **{'PASS' if not failures else 'FAIL'}**", "",
        "The primary analysis is an unrooted broad Picornavirales context tree. Candidate and reference sequences are restricted to locally detected PF00680 homologous cores. Both untrimmed and trimAl-filtered alignments were analyzed independently with model selection, ultrafast bootstrap, and SH-aLRT support.", "",
        "Picornaviridae references are not treated as a proven external outgroup for the order-wide tree. Placement remains contextual, not a taxonomic species decision.",
    ]
    (args.out / "PHYLOGENY_REPORT.md").write_text("\n".join(report) + "\n")
    sums = []
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            sums.append(f"{sha(path)}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n")
    print(json.dumps(status, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
