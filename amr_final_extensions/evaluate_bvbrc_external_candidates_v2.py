#!/usr/bin/env python3
"""Exact-ID wrapper and feasibility audit for BV-BRC external validation.

The underlying frozen-candidate evaluator is unchanged statistically. This wrapper:
1. preserves dot-containing BV-BRC genome IDs as strings when reading tab-delimited
   matrices, preventing identifiers such as 573.46140 from being coerced to floats;
2. requires exact alignment to the full QC-passed cohort; and
3. audits whether the prespecified source/country random-effects gate is actually
   estimable from the frozen cohort before interpreting zero survivors as a negative
   replication result.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from amr_final_extensions import evaluate_bvbrc_external_candidates as base


_ORIGINAL_READ_CSV = pd.read_csv


def _exact_id_read_csv(filepath_or_buffer, *args, **kwargs):
    """Preserve the first field of tab-delimited indexed matrices as exact text."""
    if kwargs.get("index_col") == 0 and kwargs.get("sep") == "\t":
        dtype = kwargs.get("dtype")
        if dtype is None:
            kwargs["dtype"] = {0: str}
        elif isinstance(dtype, dict) and 0 not in dtype:
            kwargs["dtype"] = {**dtype, 0: str}
    return _ORIGINAL_READ_CSV(filepath_or_buffer, *args, **kwargs)


def _group_table(groups: pd.Series, phenotype: pd.Series) -> tuple[list[dict], int]:
    frame = pd.DataFrame(
        {
            "group": groups.map(base.clean_group),
            "phenotype": phenotype.astype(str),
        }
    )
    frame = frame[frame["group"].ne("")]
    if frame.empty:
        return [], 0
    table = pd.crosstab(frame["group"], frame["phenotype"])
    for label in ["R", "S"]:
        if label not in table.columns:
            table[label] = 0
    table = table[["R", "S"]].astype(int)
    records = [
        {"group": str(group), "R": int(row.R), "S": int(row.S)}
        for group, row in table.iterrows()
    ]
    mixed = int(((table["R"] > 0) & (table["S"] > 0)).sum())
    return records, mixed


def _rewrite_hashes(out: Path) -> None:
    rows = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            rows.append(
                f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}"
            )
    (out / "SHA256SUMS.txt").write_text("\n".join(rows) + "\n")


def main() -> None:
    args = base.parse_args()

    # base.pd and this module's pd reference the same pandas module object.
    base.pd.read_csv = _exact_id_read_csv
    base.main()

    out = Path(args.out)
    summary_path = out / "BVBRC_EXTERNAL_VALIDATION_SUMMARY.json"
    summary = json.loads(summary_path.read_text())
    cohort = _ORIGINAL_READ_CSV(args.cohort, dtype=str).fillna("")
    cohort = cohort.drop_duplicates("genome_id").copy()

    if int(summary.get("n_external_samples", -1)) != len(cohort):
        raise RuntimeError(
            "Exact-ID alignment failed: evaluator used "
            f"{summary.get('n_external_samples')} of {len(cohort)} QC-passed genomes"
        )

    source_groups = base.group_series(
        cohort, ["bioproject_accession", "source_group", "pmid_text"]
    )
    country_groups = base.group_series(
        cohort, ["country", "geographic_location"]
    )
    source_table, source_mixed = _group_table(source_groups, cohort["phenotype"])
    country_table, country_mixed = _group_table(country_groups, cohort["phenotype"])

    min_groups = int(args.min_groups)
    feasible = bool(source_mixed >= min_groups and country_mixed >= min_groups)
    summary["exact_id_alignment"] = {
        "n_qc_passed_cohort": int(len(cohort)),
        "n_evaluated": int(summary["n_external_samples"]),
        "complete": True,
    }
    summary["external_group_gate_feasibility"] = {
        "required_mixed_groups_per_axis": min_groups,
        "source_mixed_phenotype_groups": source_mixed,
        "country_mixed_phenotype_groups": country_mixed,
        "source_group_table": source_table,
        "country_group_table": country_table,
        "feasible": feasible,
    }

    if not feasible:
        summary["raw_candidate_gate_status"] = summary.get("status")
        summary["status"] = "EXTERNAL_GATE_INFEASIBLE_SOURCE_COUNTRY_CONFOUNDING"
        summary["boundary"] = (
            "All QC-passed laboratory-method genomes were evaluated with exact string IDs. "
            f"However, the frozen cohort contains only {source_mixed} source groups and "
            f"{country_mixed} nonblank country groups with both resistant and susceptible "
            f"isolates, below the prespecified requirement of {min_groups} for each axis. "
            "Therefore zero complete-gate survivors cannot be interpreted as evidence that "
            "the frozen candidates failed independent external replication. Apparent signals "
            "remain inseparable from source/lineage structure in this cohort; a better balanced "
            "external phenotype-linked cohort is required."
        )

    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    report = [
        "# Laboratory-only BV-BRC external validation of frozen candidates",
        "",
        f"- QC-passed and evaluated genomes: **{len(cohort):,}**",
        f"- Phenotypes: **{cohort['phenotype'].value_counts().to_dict()}**",
        f"- Frozen candidates: **{summary['n_frozen_candidates']:,}**",
        f"- Represented / analysable: **{summary['n_represented_candidates']:,} / {summary['n_analysable_candidates']:,}**",
        f"- Complete external gate: **{summary['n_strict_external_replicates']:,}**",
        f"- Mixed source groups available / required: **{source_mixed} / {min_groups}**",
        f"- Mixed nonblank country groups available / required: **{country_mixed} / {min_groups}**",
        f"- Gate feasible: **{feasible}**",
        "",
        "## Claim boundary",
        "",
        summary["boundary"],
    ]
    (out / "BVBRC_EXTERNAL_VALIDATION_REPORT.md").write_text("\n".join(report) + "\n")
    _rewrite_hashes(out)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
