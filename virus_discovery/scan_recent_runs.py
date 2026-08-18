#!/usr/bin/env python3
"""Scan recent public RNA-seq runs from medicinal plants and medicinal fungi.

The script deliberately focuses on runs first released after the original 2020
Serratus sweep.  It does not call any sequence a virus.  Its sole purpose is to
identify tractable, non-virus-targeted public projects for a second-stage RdRP
screen and raw-read validation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import requests

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
CUTOFF_DATE = "2021-01-01"

# Each label may have multiple current or historical taxonomy names in ENA.
SPECIES_GROUPS = [
    ("Panax ginseng", ["Panax ginseng"]),
    ("Panax notoginseng", ["Panax notoginseng"]),
    ("Panax quinquefolius", ["Panax quinquefolius"]),
    ("Rehmannia glutinosa", ["Rehmannia glutinosa"]),
    ("Coptis chinensis", ["Coptis chinensis"]),
    ("Scutellaria baicalensis", ["Scutellaria baicalensis"]),
    ("Astragalus membranaceus", ["Astragalus membranaceus", "Astragalus mongholicus"]),
    ("Glycyrrhiza uralensis", ["Glycyrrhiza uralensis"]),
    ("Angelica sinensis", ["Angelica sinensis"]),
    ("Paeonia lactiflora", ["Paeonia lactiflora"]),
    ("Artemisia annua", ["Artemisia annua"]),
    ("Ephedra sinica", ["Ephedra sinica"]),
    ("Atractylodes macrocephala", ["Atractylodes macrocephala"]),
    ("Bupleurum chinense", ["Bupleurum chinense"]),
    ("Pueraria lobata", ["Pueraria lobata", "Pueraria montana var. lobata"]),
    ("Schisandra chinensis", ["Schisandra chinensis"]),
    ("Cornus officinalis", ["Cornus officinalis"]),
    ("Gardenia jasminoides", ["Gardenia jasminoides"]),
    ("Ziziphus jujuba", ["Ziziphus jujuba"]),
    ("Platycodon grandiflorus", ["Platycodon grandiflorus"]),
    ("Lonicera japonica", ["Lonicera japonica"]),
    ("Forsythia suspensa", ["Forsythia suspensa"]),
    ("Magnolia officinalis", ["Magnolia officinalis"]),
    ("Salvia miltiorrhiza", ["Salvia miltiorrhiza"]),
    ("Curcuma longa", ["Curcuma longa"]),
    ("Zingiber officinale", ["Zingiber officinale"]),
    ("Houttuynia cordata", ["Houttuynia cordata"]),
    ("Dendrobium officinale", ["Dendrobium officinale"]),
    ("Eucommia ulmoides", ["Eucommia ulmoides"]),
    ("Cinnamomum cassia", ["Cinnamomum cassia", "Cinnamomum aromaticum"]),
    ("Morus alba", ["Morus alba"]),
    ("Ganoderma lingzhi/lucidum", ["Ganoderma lingzhi", "Ganoderma lucidum"]),
    ("Cordyceps militaris", ["Cordyceps militaris"]),
    ("Ophiocordyceps sinensis", ["Ophiocordyceps sinensis", "Cordyceps sinensis"]),
    ("Sanghuangporus/Phellinus", ["Sanghuangporus sanghuang", "Phellinus linteus"]),
    ("Wolfiporia cocos", ["Wolfiporia cocos", "Poria cocos"]),
    ("Hericium erinaceus", ["Hericium erinaceus"]),
    ("Inonotus obliquus", ["Inonotus obliquus"]),
]

RICH_FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "experiment_accession",
    "run_accession",
    "sample_accession",
    "secondary_sample_accession",
    "scientific_name",
    "tax_id",
    "first_public",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "sample_title",
    "sample_description",
    "experiment_title",
    "study_title",
    "center_name",
    "country",
]

# Known to be accepted for read_run by ENA and therefore used as a fallback.
CORE_FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "experiment_accession",
    "run_accession",
    "secondary_sample_accession",
    "scientific_name",
    "tax_id",
    "first_public",
    "library_strategy",
    "library_source",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "sample_title",
    "sample_description",
    "library_name",
]

VIRUS_TARGET_TERMS = re.compile(
    r"\b(virus|viral|virome|virome|viridae|infection|infected|pathogen|pathogenic|"
    r"phytoplasma|mosaic disease|yellowing disease|disease resistance|challenge inoculation|"
    r"inoculated|mycovirus)\b",
    re.IGNORECASE,
)

LOW_VALUE_TERMS = re.compile(
    r"\b(single[- ]cell|scrna|spatial transcript|nanopore direct rna|ribosome profiling|ribo-seq)\b",
    re.IGNORECASE,
)


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def safe_int(value: Any) -> int:
    try:
        if value is None or value == "":
            return 0
        return int(float(str(value)))
    except (ValueError, TypeError):
        return 0


def split_numeric_sum(value: Any) -> int:
    if value is None:
        return 0
    total = 0
    for token in re.split(r"[;,]", str(value)):
        token = token.strip()
        if token:
            total += safe_int(token)
    return total


def compact_text(*values: Any) -> str:
    text = " | ".join(str(v).strip() for v in values if v not in (None, ""))
    return re.sub(r"\s+", " ", text)[:4000]


def api_get_json(session: requests.Session, params: dict[str, Any], *, attempts: int = 5) -> list[dict[str, Any]]:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(ENA_SEARCH, params=params, timeout=120)
            if response.status_code == 400:
                raise ValueError(f"ENA 400: {response.text[:500]}")
            response.raise_for_status()
            payload = response.json()
            if isinstance(payload, dict) and payload.get("message"):
                raise RuntimeError(f"ENA error: {payload['message']}")
            if not isinstance(payload, list):
                raise TypeError(f"Unexpected ENA payload type: {type(payload)!r}")
            return payload
        except Exception as exc:  # network and API failures are retried
            last_error = exc
            if attempt < attempts:
                time.sleep(min(20, 2 ** attempt))
    assert last_error is not None
    raise last_error


def query_alias(session: requests.Session, alias: str, *, limit: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    query = (
        f'scientific_name="{alias}" AND library_strategy="RNA-Seq" '
        f"AND first_public>={CUTOFF_DATE}"
    )
    common = {
        "result": "read_run",
        "query": query,
        "format": "json",
        "limit": str(limit),
    }
    rich_params = {**common, "fields": ",".join(RICH_FIELDS)}
    try:
        rows = api_get_json(session, rich_params)
        mode = "rich"
    except Exception as rich_error:
        core_params = {**common, "fields": ",".join(CORE_FIELDS)}
        rows = api_get_json(session, core_params)
        mode = "core_fallback"
        return rows, {
            "alias": alias,
            "query": query,
            "field_mode": mode,
            "rich_error": repr(rich_error),
            "returned": len(rows),
        }
    return rows, {
        "alias": alias,
        "query": query,
        "field_mode": mode,
        "returned": len(rows),
    }


def project_accession(row: dict[str, Any]) -> str:
    return str(row.get("secondary_study_accession") or row.get("study_accession") or "UNKNOWN")


def project_score(project: dict[str, Any]) -> tuple[float, list[str]]:
    n_runs = project["run_count"]
    size_gb = project["total_fastq_bytes"] / 1_000_000_000
    paired_fraction = project["paired_runs"] / max(1, n_runs)
    illumina_fraction = project["illumina_runs"] / max(1, n_runs)
    avg_read_len = project["total_base_count"] / max(1, project["total_read_count"])
    latest_year = project["latest_year"]
    text = project["search_text"]

    score = 0.0
    reasons: list[str] = []

    if 2 <= n_runs <= 16:
        score += 18
        reasons.append("replicated tractable run count")
    elif n_runs == 1:
        score += 5
        reasons.append("single run only")
    elif n_runs <= 40:
        score += 8
    else:
        score -= 10
        reasons.append("many runs")

    if 0.25 <= size_gb <= 12:
        score += 22
        reasons.append("download size suitable")
    elif 12 < size_gb <= 25:
        score += 8
        reasons.append("large but feasible")
    elif size_gb > 25:
        score -= min(25, (size_gb - 25) / 2)
        reasons.append("oversized project")
    else:
        score += 4
        reasons.append("very small project")

    score += 10 * paired_fraction
    if paired_fraction >= 0.8:
        reasons.append("mostly paired-end")
    score += 5 * illumina_fraction

    if avg_read_len >= 100:
        score += 10
        reasons.append("long Illumina reads")
    elif avg_read_len >= 75:
        score += 6
    elif avg_read_len and avg_read_len < 50:
        score -= 8
        reasons.append("short-read/possibly small-RNA data")

    if latest_year >= 2024:
        score += 12
        reasons.append("recent release")
    elif latest_year >= 2022:
        score += 7

    if VIRUS_TARGET_TERMS.search(text):
        score -= 30
        reasons.append("already virus/pathogen targeted")
    else:
        score += 15
        reasons.append("not explicitly virus-targeted")

    if LOW_VALUE_TERMS.search(text):
        score -= 18
        reasons.append("library type less suitable for assembly")

    if project["fastq_available_runs"] == n_runs:
        score += 5
        reasons.append("all FASTQ URLs available")
    elif project["fastq_available_runs"] == 0:
        score -= 20
        reasons.append("no FASTQ URLs")

    # A mild preference for projects with multiple sample titles, suggesting
    # biological replication rather than technical lane splitting.
    unique_titles = project["unique_sample_titles"]
    if unique_titles >= 2:
        score += min(8, math.log2(unique_titles + 1) * 2)
        reasons.append("multiple biological/sample labels")

    return round(score, 3), reasons


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="virus_discovery/scan_results")
    parser.add_argument("--limit-per-alias", type=int, default=1000)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "ByungwoongYoo-public-virus-discovery/0.1 (research audit; contact via GitHub)",
            "Accept": "application/json",
        }
    )

    all_rows_by_run: dict[str, dict[str, Any]] = {}
    query_log: list[dict[str, Any]] = []

    for label, aliases in SPECIES_GROUPS:
        for alias in aliases:
            print(f"[scan] {label} <- {alias}", flush=True)
            try:
                rows, log_entry = query_alias(session, alias, limit=args.limit_per_alias)
                log_entry["label"] = label
                log_entry["timestamp_utc"] = now_iso()
                query_log.append(log_entry)
            except Exception as exc:
                query_log.append(
                    {
                        "label": label,
                        "alias": alias,
                        "timestamp_utc": now_iso(),
                        "error": repr(exc),
                    }
                )
                print(f"[warn] ENA query failed for {alias}: {exc!r}", file=sys.stderr, flush=True)
                continue

            for row in rows:
                run = str(row.get("run_accession") or "")
                if not run:
                    continue
                row = dict(row)
                row["curated_species_label"] = label
                row["queried_alias"] = alias
                previous = all_rows_by_run.get(run)
                if previous is None or len(json.dumps(row)) > len(json.dumps(previous)):
                    all_rows_by_run[run] = row
            time.sleep(0.2)

    all_rows = sorted(all_rows_by_run.values(), key=lambda r: str(r.get("run_accession", "")))
    if not all_rows:
        (out / "FATAL_NO_RESULTS.txt").write_text(
            "No ENA records were returned. Inspect query_log.json for API errors.\n",
            encoding="utf-8",
        )
        (out / "query_log.json").write_text(json.dumps(query_log, indent=2), encoding="utf-8")
        return 2

    projects: dict[str, dict[str, Any]] = {}
    for row in all_rows:
        project_id = project_accession(row)
        entry = projects.setdefault(
            project_id,
            {
                "project_accession": project_id,
                "species_labels": set(),
                "scientific_names": set(),
                "runs": [],
                "sample_titles": set(),
                "search_fragments": [],
                "total_fastq_bytes": 0,
                "total_read_count": 0,
                "total_base_count": 0,
                "paired_runs": 0,
                "illumina_runs": 0,
                "fastq_available_runs": 0,
                "years": [],
                "centers": set(),
                "countries": set(),
            },
        )
        entry["species_labels"].add(str(row.get("curated_species_label", "")))
        entry["scientific_names"].add(str(row.get("scientific_name", "")))
        entry["runs"].append(str(row.get("run_accession", "")))
        if row.get("sample_title"):
            entry["sample_titles"].add(str(row["sample_title"]))
        entry["search_fragments"].append(
            compact_text(
                row.get("study_title"),
                row.get("experiment_title"),
                row.get("sample_title"),
                row.get("sample_description"),
                row.get("library_name"),
            )
        )
        entry["total_fastq_bytes"] += split_numeric_sum(row.get("fastq_bytes"))
        entry["total_read_count"] += safe_int(row.get("read_count"))
        entry["total_base_count"] += safe_int(row.get("base_count"))
        if str(row.get("library_layout", "")).upper() == "PAIRED":
            entry["paired_runs"] += 1
        if str(row.get("instrument_platform", "")).upper() == "ILLUMINA":
            entry["illumina_runs"] += 1
        if row.get("fastq_ftp"):
            entry["fastq_available_runs"] += 1
        match = re.match(r"(\d{4})", str(row.get("first_public", "")))
        if match:
            entry["years"].append(int(match.group(1)))
        if row.get("center_name"):
            entry["centers"].add(str(row["center_name"]))
        if row.get("country"):
            entry["countries"].add(str(row["country"]))

    ranked: list[dict[str, Any]] = []
    for project_id, entry in projects.items():
        run_count = len(set(entry["runs"]))
        years = entry["years"] or [0]
        search_text = compact_text(*entry["search_fragments"])
        normalized = {
            "project_accession": project_id,
            "species_labels": "; ".join(sorted(x for x in entry["species_labels"] if x)),
            "scientific_names": "; ".join(sorted(x for x in entry["scientific_names"] if x)),
            "run_count": run_count,
            "run_accessions": ";".join(sorted(set(entry["runs"]))),
            "total_fastq_bytes": entry["total_fastq_bytes"],
            "total_fastq_gb": round(entry["total_fastq_bytes"] / 1_000_000_000, 3),
            "total_read_count": entry["total_read_count"],
            "total_base_count": entry["total_base_count"],
            "estimated_mean_read_length": round(
                entry["total_base_count"] / max(1, entry["total_read_count"]), 2
            ),
            "paired_runs": entry["paired_runs"],
            "illumina_runs": entry["illumina_runs"],
            "fastq_available_runs": entry["fastq_available_runs"],
            "earliest_year": min(years),
            "latest_year": max(years),
            "unique_sample_titles": len(entry["sample_titles"]),
            "sample_title_examples": " || ".join(sorted(entry["sample_titles"])[:8]),
            "centers": "; ".join(sorted(entry["centers"])),
            "countries": "; ".join(sorted(entry["countries"])),
            "search_text": search_text,
        }
        score, reasons = project_score(normalized)
        normalized["screen_priority_score"] = score
        normalized["priority_reasons"] = "; ".join(reasons)
        normalized["explicit_virus_target"] = bool(VIRUS_TARGET_TERMS.search(search_text))
        ranked.append(normalized)

    ranked.sort(
        key=lambda r: (
            -float(r["screen_priority_score"]),
            float(r["total_fastq_gb"]),
            str(r["project_accession"]),
        )
    )
    for index, row in enumerate(ranked, start=1):
        row["rank"] = index

    run_fields = sorted({key for row in all_rows for key in row.keys()})
    project_fields = [
        "rank",
        "screen_priority_score",
        "project_accession",
        "species_labels",
        "scientific_names",
        "run_count",
        "total_fastq_gb",
        "total_fastq_bytes",
        "total_read_count",
        "total_base_count",
        "estimated_mean_read_length",
        "paired_runs",
        "illumina_runs",
        "fastq_available_runs",
        "earliest_year",
        "latest_year",
        "unique_sample_titles",
        "explicit_virus_target",
        "priority_reasons",
        "sample_title_examples",
        "centers",
        "countries",
        "run_accessions",
        "search_text",
    ]

    write_tsv(out / "all_runs.tsv", all_rows, run_fields)
    write_tsv(out / "ranked_projects.tsv", ranked, project_fields)
    write_tsv(out / "top30_projects.tsv", ranked[:30], project_fields)

    species_counts = Counter(row.get("curated_species_label", "") for row in all_rows)
    species_rows = [
        {"curated_species_label": label, "run_count": count}
        for label, count in sorted(species_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    write_tsv(out / "species_summary.tsv", species_rows, ["curated_species_label", "run_count"])

    (out / "query_log.json").write_text(
        json.dumps(query_log, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (out / "top30_projects.json").write_text(
        json.dumps(ranked[:30], ensure_ascii=False, indent=2), encoding="utf-8"
    )

    md_lines = [
        "# Recent medicinal-species RNA-seq candidate scan",
        "",
        f"Generated: {now_iso()}",
        f"Cutoff: first_public >= {CUTOFF_DATE}",
        f"Unique runs: {len(all_rows):,}",
        f"Projects: {len(ranked):,}",
        "",
        "This table is a metadata triage only. No project or sequence is called viral at this stage.",
        "",
        "| Rank | Score | Project | Species | Runs | FASTQ GB | Mean read nt | Latest | Virus-targeted? |",
        "|---:|---:|---|---|---:|---:|---:|---:|---|",
    ]
    for row in ranked[:30]:
        md_lines.append(
            f"| {row['rank']} | {row['screen_priority_score']} | {row['project_accession']} | "
            f"{row['species_labels'][:55]} | {row['run_count']} | {row['total_fastq_gb']} | "
            f"{row['estimated_mean_read_length']} | {row['latest_year']} | "
            f"{row['explicit_virus_target']} |"
        )
    (out / "SCAN_REPORT.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    manifest_lines = []
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_lines.append(f"{digest}  {path.name}")
    (out / "SHA256SUMS.txt").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print(f"[done] {len(all_rows)} unique runs across {len(ranked)} projects")
    print(f"[done] results: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
