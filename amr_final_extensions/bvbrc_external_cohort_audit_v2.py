#!/usr/bin/env python3
"""Compatibility wrapper for the candidate-blind BV-BRC cohort audit.

BV-BRC collections do not share a universal `id` sort field. The genome_amr
collection supports `id`, whereas the genome collection must be sorted by
`genome_id`. This wrapper replaces only the HTTP request helper and then runs the
frozen audit unchanged.
"""
from __future__ import annotations

import time

import requests

from amr_final_extensions import bvbrc_external_cohort_audit as base


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
    suffix = (
        f"{expression}&select({','.join(fields)})&sort(%2B{sort_field})"
        f"&limit({limit},{start})"
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
                return payload, {
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
                }
            errors.append(f"GET {response.status_code}: {response.text[:500]}")
        except Exception as exc:
            errors.append(f"GET exception: {exc!r}")

        try:
            body = (
                f"{expression}&select({','.join(fields)})&sort(+{sort_field})"
                f"&limit({limit},{start})"
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
                return payload, {
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
                }
            errors.append(f"POST {response.status_code}: {response.text[:500]}")
        except Exception as exc:
            errors.append(f"POST exception: {exc!r}")
        if attempt < attempts:
            time.sleep(5 * attempt)
    raise RuntimeError("BV-BRC query failed after retries:\n" + "\n".join(errors))


base.request_rql = request_rql

if __name__ == "__main__":
    base.main()
