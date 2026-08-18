#!/usr/bin/env python3
"""Select recent public RNA-seq runs from medicinal-plant taxa.

The script queries NCBI E-utilities and writes a compact GitHub Actions
matrix. It deliberately favors recent, modest-sized, replicated projects.
No biological claim is made from metadata alone.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import time
from collections import defaultdict
from datetime import datetime
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "public-virome-hunt/0.1 (independent research; contact via GitHub)"

TARGETS = [
    ("panax", '("Panax ginseng"[Organism] OR "Panax quinquefolius"[Organism] OR "Panax notoginseng"[Organism])'),
    ("platycodon", '"Platycodon grandiflorus"[Organism]'),
    ("scutellaria", '"Scutellaria baicalensis"[Organism]'),
    ("glycyrrhiza", '("Glycyrrhiza uralensis"[Organism] OR "Glycyrrhiza glabra"[Organism])'),
    ("rehmannia", '"Rehmannia glutinosa"[Organism]'),
    ("schisandra", '"Schisandra chinensis"[Organism]'),
    ("astragalus", '("Astragalus mongholicus"[Organism] OR "Astragalus membranaceus"[Organism])'),
    ("angelica", '("Angelica gigas"[Organism] OR "Angelica sinensis"[Organism])'),
]

# Controls/fallbacks. SRR9968562 is the public PalmID tutorial example with an
# RdRP-containing microassembly. The Panax runs provide a Korean-medicine-linked
# legacy cohort even if recent-query metadata retrieval fails.
FALLBACK = [
    {"label": "positive_control", "accession": "SRR9968562", "bioproject": "unknown", "organism": "Waxsystermes-associated sample", "release_date": "legacy", "size_mb": "", "bases": "", "reason": "PalmID tutorial positive control"},
    {"label": "panax_legacy", "accession": "SRR4835279", "bioproject": "unknown", "organism": "Panax ginseng", "release_date": "2017", "size_mb": "", "bases": "", "reason": "legacy medicinal-plant probe"},
    {"label": "panax_legacy", "accession": "SRR4835278", "bioproject": "unknown", "organism": "Panax ginseng", "release_date": "2017", "size_mb": "", "bases": "", "reason": "legacy medicinal-plant replicate"},
]


def get(url: str, params: dict[str, Any], retries: int = 4) -> bytes:
    full = url + "?" + urlencode(params)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(full, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
            with urlopen(req, timeout=120) as response:
                return response.read()
        except Exception as exc:  # network/transient API failures
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET failed after {retries} attempts: {full}: {last}")


def esearch(term: str, retmax: int = 200) -> list[str]:
    raw = get(
        f"{EUTILS}/esearch.fcgi",
        {"db": "sra", "term": term, "retmax": retmax, "retmode": "json", "sort": "date"},
    )
    obj = json.loads(raw)
    return obj.get("esearchresult", {}).get("idlist", [])


def efetch_runinfo(ids: list[str]) -> list[dict[str, str]]:
    if not ids:
        return []
    raw = get(
        f"{EUTILS}/efetch.fcgi",
        {"db": "sra", "id": ",".join(ids), "rettype": "runinfo", "retmode": "text"},
    )
    text = raw.decode("utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(text)))


def as_float(value: str | None, default: float = 1e18) -> float:
    try:
        return float(value or "")
    except ValueError:
        return default


def select_for_target(label: str, organism_clause: str, max_runs: int = 3) -> list[dict[str, str]]:
    # PDAT is intentionally broad; final filtering is done from returned RunInfo.
    term = (
        f"({organism_clause}) AND (\"rna seq\"[Strategy] OR \"transcriptomic\"[Source]) "
        'AND ("2024/01/01"[PDAT] : "2025/12/31"[PDAT])'
    )
    ids = esearch(term)
    rows = efetch_runinfo(ids)
    clean: list[dict[str, str]] = []
    for row in rows:
        accession = row.get("Run", "").strip()
        if not accession.startswith(("SRR", "ERR", "DRR")):
            continue
        strategy = (row.get("LibraryStrategy") or "").upper()
        source = (row.get("LibrarySource") or "").upper()
        if strategy not in {"RNA-SEQ", "OTHER"} and source != "TRANSCRIPTOMIC":
            continue
        size = as_float(row.get("size_MB"))
        # Avoid exceptionally large assemblies in the first exhaustive pass.
        if size > 4500:
            continue
        clean.append(row)

    # Prefer a project with multiple runs, then the latest and smallest runs.
    by_project: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in clean:
        by_project[row.get("BioProject") or "unknown"].append(row)
    projects = sorted(
        by_project.items(),
        key=lambda kv: (
            -min(len(kv[1]), 10),
            max((r.get("ReleaseDate") or "") for r in kv[1]),
        ),
        reverse=True,
    )
    chosen: list[dict[str, str]] = []
    if projects:
        project, candidates = projects[0]
        candidates.sort(key=lambda r: (as_float(r.get("size_MB")), r.get("Run") or ""))
        for row in candidates[:max_runs]:
            chosen.append(
                {
                    "label": label,
                    "accession": row.get("Run", ""),
                    "bioproject": project,
                    "organism": row.get("ScientificName", ""),
                    "release_date": row.get("ReleaseDate", ""),
                    "size_mb": row.get("size_MB", ""),
                    "bases": row.get("bases", ""),
                    "reason": "recent replicated medicinal-plant RNA-seq",
                }
            )
    return chosen


def main() -> int:
    selected: list[dict[str, str]] = []
    audit: list[dict[str, Any]] = []
    for label, clause in TARGETS:
        try:
            rows = select_for_target(label, clause)
            selected.extend(rows)
            audit.append({"label": label, "status": "ok", "selected": len(rows)})
        except Exception as exc:
            audit.append({"label": label, "status": "error", "error": str(exc)})
        time.sleep(0.35)

    # Always retain the explicit positive control and at least one Panax probe.
    seen = {x["accession"] for x in selected}
    for row in FALLBACK:
        if row["accession"] not in seen:
            selected.append(row)
            seen.add(row["accession"])

    # Hard cap to keep the first run within GitHub-hosted runner limits.
    # Keep every target represented before filling remaining slots.
    first_by_label: list[dict[str, str]] = []
    rest: list[dict[str, str]] = []
    label_seen: set[str] = set()
    for row in selected:
        if row["label"] not in label_seen:
            first_by_label.append(row)
            label_seen.add(row["label"])
        else:
            rest.append(row)
    selected = (first_by_label + rest)[:18]

    with open("selected_runs.tsv", "w", encoding="utf-8", newline="") as handle:
        fields = ["label", "accession", "bioproject", "organism", "release_date", "size_mb", "bases", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    with open("selection_audit.json", "w", encoding="utf-8") as handle:
        json.dump(
            {"generated_utc": datetime.utcnow().isoformat() + "Z", "queries": audit, "selected": selected},
            handle,
            indent=2,
            ensure_ascii=False,
        )

    matrix = {"include": selected}
    with open("matrix.json", "w", encoding="utf-8") as handle:
        json.dump(matrix, handle, separators=(",", ":"), ensure_ascii=False)
    print(json.dumps(matrix, separators=(",", ":"), ensure_ascii=False))
    return 0 if selected else 2


if __name__ == "__main__":
    raise SystemExit(main())
