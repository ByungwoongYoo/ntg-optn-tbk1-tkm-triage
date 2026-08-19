#!/usr/bin/env python3
"""Compatibility wrapper for the candidate-blind BV-BRC cohort audit.

BV-BRC collections do not share a universal `id` sort field. The genome_amr
collection supports `id`, whereas the genome collection must be sorted by
`genome_id`. In addition, the public API caps a response at about 25,000 rows even
when a larger RQL limit is requested. This wrapper applies collection-specific sorting,
paginates deterministically, records every page, de-duplicates exact record IDs, and
then runs the frozen audit unchanged.
"""
from __future__ import annotations

import hashlib
import re
import sys
import time
from pathlib import Path

import requests

# When invoked as `python amr_final_extensions/<file>.py`, Python places only the
# script directory on sys.path. Add the repository root explicitly so the sibling
# package can be imported reproducibly on GitHub Actions and local clean rooms.
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

    # Stable de-duplication is a safety guard for API pagination boundary behavior.
    if records and "id" in records[0]:
        seen: set[str] = set()
        deduped: list[dict] = []
        for record in records:
            key = str(record.get("id", ""))
            if key and key in seen:
                continue
            if key:
                seen.add(key)
            deduped.append(record)
        records = deduped
    elif records and "genome_id" in records[0]:
        seen = set()
        deduped = []
        for record in records:
            key = str(record.get("genome_id", ""))
            if key and key in seen:
                continue
            if key:
                seen.add(key)
            deduped.append(record)
        records = deduped

    combined_hash = hashlib.sha256("\n".join(page_hash_material).encode()).hexdigest()
    meta = {
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
    return records, meta


base.request_rql = request_rql

if __name__ == "__main__":
    base.main()
