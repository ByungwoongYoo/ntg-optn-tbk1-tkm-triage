#!/usr/bin/env python3
"""Independent HMMER-web audit of four Panax-associated RdRP candidates.

Runs phmmer against the current UniProtKB target database and hmmscan against
Pfam. Raw API responses are preserved. The script does not convert a no-hit
result into proof of novelty and does not infer a biological host.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import requests
from Bio import SeqIO

BASE = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def post_job(session: requests.Session, algo: str, database: str, fasta: str) -> dict[str, Any]:
    url = f"{BASE}/search/{algo}"
    body = {"database": database, "input": fasta}
    last: Exception | None = None
    for attempt in range(1, 6):
        try:
            r = session.post(url, json=body, timeout=120)
            r.raise_for_status()
            data = r.json()
            if not data.get("id"):
                raise RuntimeError(f"missing job id: {data}")
            return data
        except Exception as exc:
            last = exc
            if attempt < 5:
                time.sleep(min(60, 3 * 2 ** attempt))
    assert last is not None
    raise last


def poll_job(session: requests.Session, job_id: str, max_minutes: int = 45) -> dict[str, Any]:
    url = f"{BASE}/result/{job_id}"
    deadline = time.time() + max_minutes * 60
    last_payload: dict[str, Any] = {}
    while time.time() < deadline:
        r = session.get(url, headers={"Accept": "application/json"}, timeout=120)
        if r.status_code in (429, 502, 503, 504):
            time.sleep(15)
            continue
        r.raise_for_status()
        payload = r.json()
        last_payload = payload
        status = str(payload.get("status", "")).upper()
        if status == "SUCCESS":
            return payload
        if status in {"FAILURE", "ERROR", "CANCELLED"}:
            raise RuntimeError(f"HMMER job {job_id} ended with {status}: {payload}")
        time.sleep(8)
    raise TimeoutError(f"HMMER job {job_id} not complete; last payload={last_payload}")


def collect_hits(node: Any) -> list[dict[str, Any]]:
    """Find hit-like dictionaries recursively without relying on one API schema."""
    hits: list[dict[str, Any]] = []
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "hits" and isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        hits.append(item)
            else:
                hits.extend(collect_hits(value))
    elif isinstance(node, list):
        for value in node:
            hits.extend(collect_hits(value))
    return hits


def first_value(hit: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in hit and hit[key] not in (None, ""):
            return hit[key]
    return ""


def summarize(query: str, algo: str, database: str, payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for hit in collect_hits(payload):
        accession = str(first_value(hit, "acc", "accession", "acc2"))
        name = str(first_value(hit, "name", "id", "target", "hit_id"))
        description = str(first_value(hit, "desc", "description", "title"))
        evalue = first_value(hit, "evalue", "E-value", "eval")
        score = first_value(hit, "score", "bitscore", "bit_score")
        species = str(first_value(hit, "species", "taxname", "kingdom"))
        key = (accession, name, str(evalue))
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "query": query,
            "algorithm": algo,
            "database": database,
            "accession": accession,
            "name": name,
            "description": description,
            "species_or_taxon": species,
            "evalue": evalue,
            "score": score,
            "raw_hit_json": json.dumps(hit, ensure_ascii=False, sort_keys=True),
        })
    def number(x: Any, default: float = 1e999) -> float:
        try: return float(x)
        except Exception: return default
    rows.sort(key=lambda r: (number(r["evalue"]), -number(r["score"], -1e999)))
    return rows[:100]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "ByungwoongYoo-Panax-RdRP-audit/0.1",
        "Accept": "application/json",
    })

    all_rows: list[dict[str, Any]] = []
    status: list[dict[str, Any]] = []
    records = list(SeqIO.parse(str(args.fasta), "fasta"))
    if not records:
        raise SystemExit("No FASTA records")

    searches = [("phmmer", "uniprotkb"), ("hmmscan", "pfam")]
    for record in records:
        fasta = f">{record.id}\n{str(record.seq)}\n"
        for algo, database in searches:
            stem = f"{record.id}.{algo}.{database}"
            entry: dict[str, Any] = {
                "query": record.id,
                "algorithm": algo,
                "database": database,
                "submitted_utc": utc(),
            }
            try:
                submission = post_job(session, algo, database, fasta)
                entry["job_id"] = submission["id"]
                (args.out / f"{stem}.submission.json").write_text(
                    json.dumps(submission, indent=2) + "\n", encoding="utf-8"
                )
                payload = poll_job(session, submission["id"])
                (args.out / f"{stem}.result.json").write_text(
                    json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                rows = summarize(record.id, algo, database, payload)
                all_rows.extend(rows)
                entry.update({
                    "success": True,
                    "completed_utc": utc(),
                    "summarized_hit_count": len(rows),
                    "reported_status": payload.get("status"),
                })
            except Exception as exc:
                entry.update({"success": False, "completed_utc": utc(), "error": repr(exc)})
            status.append(entry)

    fields = [
        "query", "algorithm", "database", "accession", "name", "description",
        "species_or_taxon", "evalue", "score", "raw_hit_json",
    ]
    with (args.out / "HMMER_TOP_HITS.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(all_rows)
    (args.out / "HMMER_JOB_STATUS.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    report = [
        "# Independent HMMER-web audit",
        "",
        f"Generated: {utc()}",
        "",
        "This is a current UniProtKB/Pfam homology and domain audit. Empty results are not proof of novelty.",
        "",
        "| Query | Search | Success | Summarized hits | Job ID |",
        "|---|---|---|---:|---|",
    ]
    for x in status:
        report.append(
            f"| `{x['query']}` | {x['algorithm']} / {x['database']} | {x['success']} | "
            f"{x.get('summarized_hit_count','NA')} | `{x.get('job_id','')}` |"
        )
    (args.out / "HMMER_AUDIT_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = []
    for path in sorted(args.out.glob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            manifest.append(f"{sha256(path)}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print((args.out / "HMMER_AUDIT_REPORT.md").read_text())
    return 0 if all(x.get("success") for x in status) else 1


if __name__ == "__main__":
    raise SystemExit(main())
