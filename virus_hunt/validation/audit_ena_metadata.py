#!/usr/bin/env python3
"""Fetch ENA metadata for the Panax candidate runs and audit independence.

The audit distinguishes archived run, experiment, and sample accessions.  It does
not infer independent plants or biological replication unless the submitted
metadata explicitly supports that conclusion.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RICH_FIELDS = [
    "run_accession", "study_accession", "secondary_study_accession",
    "sample_accession", "secondary_sample_accession", "experiment_accession",
    "scientific_name", "tax_id", "sample_title", "sample_description",
    "library_name", "library_strategy", "library_source", "library_selection",
    "library_layout", "instrument_platform", "instrument_model", "center_name",
    "first_public", "last_updated", "country", "collection_date",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
]
FALLBACK_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "experiment_accession",
    "scientific_name", "sample_title", "library_name", "library_strategy",
    "library_source", "library_selection", "library_layout", "instrument_platform",
    "instrument_model", "first_public", "read_count", "base_count",
    "fastq_ftp", "fastq_md5", "fastq_bytes",
]


def fetch_one(session: requests.Session, run: str) -> tuple[dict[str, str], dict[str, Any]]:
    errors: list[str] = []
    for mode, fields in (("rich", RICH_FIELDS), ("fallback", FALLBACK_FIELDS)):
        params = {
            "accession": run,
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
        }
        for attempt in range(5):
            try:
                response = session.get(ENA, params=params, timeout=120)
                if response.status_code >= 400:
                    raise RuntimeError(f"HTTP {response.status_code}: {response.text[:300]}")
                rows = list(csv.DictReader(io.StringIO(response.text), delimiter="\t"))
                if not rows:
                    raise RuntimeError("ENA returned no data rows")
                row = {k: (v or "") for k, v in rows[0].items()}
                row["metadata_field_mode"] = mode
                return row, {
                    "run": run,
                    "mode": mode,
                    "request_url": response.url,
                    "status": "ok",
                }
            except Exception as exc:
                errors.append(f"{mode}/attempt{attempt+1}: {exc!r}")
                time.sleep(min(20, 2 ** attempt))
    return {"run_accession": run, "metadata_field_mode": "failed"}, {
        "run": run,
        "status": "failed",
        "errors": errors,
    }


def unique(rows: list[dict[str, str]], *keys: str) -> list[str]:
    values: set[str] = set()
    for row in rows:
        for key in keys:
            value = (row.get(key) or "").strip()
            if value:
                values.add(value)
    return sorted(values)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Panax-RdRP-validation/0.2 (independent public-data audit)",
        "Accept": "text/tab-separated-values",
    })

    rows: list[dict[str, str]] = []
    requests_log: list[dict[str, Any]] = []
    for run in args.runs:
        row, log = fetch_one(session, run)
        rows.append(row)
        requests_log.append(log)
        time.sleep(0.4)

    fieldnames = sorted({k for row in rows for k in row})
    with (args.out_dir / "ENA_RUN_METADATA.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    (args.out_dir / "ENA_RUN_METADATA.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (args.out_dir / "ENA_REQUEST_LOG.json").write_text(
        json.dumps(requests_log, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    run_ids = unique(rows, "run_accession")
    experiment_ids = unique(rows, "experiment_accession")
    sample_ids = unique(rows, "secondary_sample_accession", "sample_accession")
    study_ids = unique(rows, "secondary_study_accession", "study_accession")
    library_names = unique(rows, "library_name")
    sample_titles = unique(rows, "sample_title")
    failed_runs = [r.get("run_accession", "") for r in rows if r.get("metadata_field_mode") == "failed"]

    if failed_runs:
        grade = "metadata_incomplete"
    elif len(sample_ids) >= 2 and len(experiment_ids) >= 2:
        grade = "distinct_archived_samples_and_experiments"
    elif len(experiment_ids) >= 2:
        grade = "distinct_experiments_sample_independence_unresolved"
    elif len(run_ids) >= 2:
        grade = "distinct_runs_only"
    else:
        grade = "single_archived_run"

    audit = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "requested_runs": args.runs,
        "metadata_complete": not bool(failed_runs),
        "failed_runs": failed_runs,
        "run_accessions": run_ids,
        "experiment_accessions": experiment_ids,
        "sample_accessions": sample_ids,
        "study_accessions": study_ids,
        "library_names": library_names,
        "sample_titles": sample_titles,
        "distinct_run_count": len(run_ids),
        "distinct_experiment_count": len(experiment_ids),
        "distinct_archived_sample_count": len(sample_ids),
        "independence_grade": grade,
        "interpretation_boundary": (
            "Distinct archive sample or experiment accessions do not by themselves prove "
            "independent plants, independent infections, or the true biological host."
        ),
    }
    (args.out_dir / "INDEPENDENCE_AUDIT.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# ENA metadata and independence audit",
        "",
        f"Generated (UTC): {audit['generated_utc']}",
        "",
        f"- Independence grade: **`{grade}`**",
        f"- Distinct runs: **{len(run_ids)}**",
        f"- Distinct experiments: **{len(experiment_ids)}**",
        f"- Distinct archived samples: **{len(sample_ids)}**",
        f"- Study accessions: {', '.join(study_ids) or 'not returned'}",
        "",
        "This audit does not infer independent plants or a Panax host unless the submitted metadata explicitly establishes those facts.",
        "",
        "| Run | Experiment | Archived sample | Sample title | Library | First public |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        sample = row.get("secondary_sample_accession") or row.get("sample_accession") or ""
        lines.append(
            f"| `{row.get('run_accession','')}` | `{row.get('experiment_accession','')}` | "
            f"`{sample}` | {row.get('sample_title','')[:80]} | {row.get('library_name','')[:60]} | "
            f"{row.get('first_public','')} |"
        )
    (args.out_dir / "INDEPENDENCE_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0 if not failed_runs else 1


if __name__ == "__main__":
    raise SystemExit(main())
