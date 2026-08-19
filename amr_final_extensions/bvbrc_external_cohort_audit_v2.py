#!/usr/bin/env python3
"""Compatibility and evidence-boundary wrapper for the BV-BRC cohort audit.

BV-BRC collections do not share a universal sort field and the public API caps a
response at about 25,000 rows. This wrapper applies collection-specific sorting,
paginates deterministically, records every page, and de-duplicates stable record IDs.

Crucially, after the base audit finishes, this wrapper separately freezes a balanced
BioProject-disjoint cohort restricted to records whose BV-BRC evidence field is
`Laboratory Method`. Records labelled `Computational Method` are retained only in the
availability audit and are never promoted as independent phenotype validation data.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amr_final_extensions import bvbrc_external_cohort_audit as base


def _request_page(
    collection: str,
    expression: str,
    fields: list[str],
    *,
    sort_field: str,
    page_size: int,
    start: int,
    timeout: int,
    attempts: int,
) -> tuple[list[dict], dict, bytes]:
    suffix = (
        f"{expression}&select({','.join(fields)})&sort(%2B{sort_field})"
        f"&limit({page_size},{start})"
    )
    url = f"{base.API}/{collection}/?{suffix}"
    headers = {"Accept": "application/json", "User-Agent": base.UA}
    errors: list[str] = []
    for attempt in range(1, attempts + 1):
        try:
            response = requests.get(url, headers=headers, timeout=timeout)
            if response.status_code == 200:
                payload = response.json()
                if not isinstance(payload, list):
                    raise RuntimeError(f"Unexpected {collection} payload type: {type(payload)}")
                meta = {
                    "method": "GET",
                    "url": url,
                    "sort_field": sort_field,
                    "status": response.status_code,
                    "content_range": response.headers.get("Content-Range"),
                    "etag": response.headers.get("ETag"),
                    "date": response.headers.get("Date"),
                    "sha256": base.sha256_bytes(response.content),
                    "bytes": len(response.content),
                    "n_records": len(payload),
                    "start": start,
                    "requested_page_size": page_size,
                }
                return payload, meta, response.content
            errors.append(f"GET {response.status_code}: {response.text[:500]}")
        except Exception as exc:
            errors.append(f"GET exception: {exc!r}")

        try:
            body = (
                f"{expression}&select({','.join(fields)})&sort(+{sort_field})"
                f"&limit({page_size},{start})"
            )
            response = requests.post(
                f"{base.API}/{collection}/",
                data=body.encode("utf-8"),
                headers={
                    **headers,
                    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
                },
                timeout=timeout,
            )
            if response.status_code == 200:
                payload = response.json()
                if not isinstance(payload, list):
                    raise RuntimeError(f"Unexpected {collection} POST payload type: {type(payload)}")
                meta = {
                    "method": "POST",
                    "url": f"{base.API}/{collection}/",
                    "query": body,
                    "sort_field": sort_field,
                    "status": response.status_code,
                    "content_range": response.headers.get("Content-Range"),
                    "etag": response.headers.get("ETag"),
                    "date": response.headers.get("Date"),
                    "sha256": base.sha256_bytes(response.content),
                    "bytes": len(response.content),
                    "n_records": len(payload),
                    "start": start,
                    "requested_page_size": page_size,
                }
                return payload, meta, response.content
            errors.append(f"POST {response.status_code}: {response.text[:500]}")
        except Exception as exc:
            errors.append(f"POST exception: {exc!r}")
        if attempt < attempts:
            time.sleep(5 * attempt)
    raise RuntimeError("BV-BRC page query failed after retries:\n" + "\n".join(errors))


def _total_from_range(value: str | None) -> int | None:
    if not value:
        return None
    match = re.search(r"/(\d+)\s*$", value)
    return int(match.group(1)) if match else None


def request_rql(
    collection: str,
    expression: str,
    fields: list[str],
    *,
    limit: int = 250000,
    start: int = 0,
    timeout: int = 180,
    attempts: int = 4,
) -> tuple[list[dict], dict]:
    sort_field = "genome_id" if collection == "genome" else "id"
    api_cap = 25000
    requested_total = max(0, int(limit))
    cursor = max(0, int(start))
    records: list[dict] = []
    pages: list[dict] = []
    page_hash_material: list[str] = []
    known_total: int | None = None

    while len(records) < requested_total:
        remaining = requested_total - len(records)
        page_size = min(api_cap, remaining)
        payload, meta, _raw = _request_page(
            collection,
            expression,
            fields,
            sort_field=sort_field,
            page_size=page_size,
            start=cursor,
            timeout=timeout,
            attempts=attempts,
        )
        pages.append(meta)
        page_hash_material.append(meta["sha256"])
        if known_total is None:
            known_total = _total_from_range(meta.get("content_range"))
        if not payload:
            break
        records.extend(payload)
        cursor += len(payload)
        if len(payload) < page_size:
            break
        if known_total is not None and cursor >= known_total:
            break

    dedup_key = "id" if records and "id" in records[0] else "genome_id"
    seen: set[str] = set()
    deduped: list[dict] = []
    for record in records:
        key = str(record.get(dedup_key, ""))
        if key and key in seen:
            continue
        if key:
            seen.add(key)
        deduped.append(record)
    records = deduped

    combined_hash = hashlib.sha256("\n".join(page_hash_material).encode()).hexdigest()
    return records, {
        "method": "PAGINATED_GET_POST_FALLBACK",
        "collection": collection,
        "expression": expression,
        "fields": fields,
        "sort_field": sort_field,
        "requested_limit": requested_total,
        "requested_start": start,
        "api_page_cap": api_cap,
        "reported_total": known_total,
        "n_pages": len(pages),
        "n_records": len(records),
        "sha256": combined_hash,
        "pages": pages,
    }


def freeze_laboratory_only() -> None:
    args = base.parse_args()
    out = Path(args.out)
    strict_path = out / "BVBRC_BIOPROJECT_DISJOINT.csv"
    if not strict_path.exists():
        raise FileNotFoundError(strict_path)
    strict = pd.read_csv(strict_path, dtype=str).fillna("")
    laboratory = strict[strict["evidence"].astype(str).eq("Laboratory Method")].copy()
    laboratory.to_csv(out / "BVBRC_BIOPROJECT_DISJOINT_LABORATORY_METHOD.csv", index=False)
    counts = laboratory["phenotype"].value_counts().to_dict()
    balanced_n = min(int(args.per_class), int(counts.get("R", 0)), int(counts.get("S", 0)))
    frozen = base.select_diverse_balanced(laboratory, balanced_n, args.seed)
    frozen.to_csv(out / "BVBRC_EXTERNAL_LABORATORY_METHOD_FROZEN_COHORT.csv", index=False)
    (out / "BVBRC_EXTERNAL_LABORATORY_METHOD_GENOME_IDS.txt").write_text(
        "\n".join(frozen.get("genome_id", pd.Series(dtype=str)).astype(str))
        + ("\n" if len(frozen) else "")
    )

    summary_path = out / "BVBRC_EXTERNAL_COHORT_SUMMARY.json"
    summary = json.loads(summary_path.read_text())
    frozen_counts = frozen["phenotype"].value_counts().to_dict()
    summary.update({
        "bioproject_disjoint_laboratory_method": int(len(laboratory)),
        "bioproject_disjoint_laboratory_method_counts": counts,
        "laboratory_method_balanced_n_per_class": int(balanced_n),
        "laboratory_method_frozen": int(len(frozen)),
        "laboratory_method_frozen_counts": frozen_counts,
        "laboratory_method_external_validation_feasible": bool(
            frozen_counts.get("R", 0) >= 50 and frozen_counts.get("S", 0) >= 50
        ),
        "computational_method_excluded_from_external_validation": True,
        "external_validation_boundary": (
            "Only records explicitly labelled Laboratory Method are eligible for independent "
            "phenotype validation. Computational Method records are excluded from candidate testing."
        ),
    })
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    report_path = out / "BVBRC_EXTERNAL_COHORT_REPORT.md"
    with report_path.open("a", encoding="utf-8") as handle:
        handle.write("\n## Laboratory-evidence gate\n\n")
        handle.write(f"- BioProject-disjoint laboratory-method genomes: **{len(laboratory):,}** {counts}\n")
        handle.write(f"- Candidate-blind balanced frozen laboratory cohort: **{len(frozen):,}** {frozen_counts}\n")
        handle.write(f"- At least 50 per class: **{summary['laboratory_method_external_validation_feasible']}**\n\n")
        handle.write("Records labelled Computational Method are retained in the availability audit but excluded from external phenotype validation.\n")

    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps({
        "laboratory_method_records": len(laboratory),
        "laboratory_method_counts": counts,
        "frozen_records": len(frozen),
        "frozen_counts": frozen_counts,
    }, indent=2))


base.request_rql = request_rql

if __name__ == "__main__":
    base.main()
    freeze_laboratory_only()
