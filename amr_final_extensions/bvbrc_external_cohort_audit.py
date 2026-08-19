#!/usr/bin/env python3
"""Freeze an independent BV-BRC K. pneumoniae-colistin phenotype cohort.

This script queries the public BV-BRC Data API before any candidate sequence is tested.
It audits measured colistin phenotypes, joins genome metadata, removes overlap with the
AMR Portal discovery corpus at assembly/BioSample/BioProject levels, and freezes
source-diverse R/S cohorts when the public data permit it.

The result is a cohort-availability and independence audit. It does not test candidate
markers and cannot establish novelty, causality, diagnostic validity, or treatment use.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

API = "https://www.bv-brc.org/api"
UA = "ByungwoongYoo-AMR-audit/20260819 (public reproducibility study)"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--portal-labels", required=True)
    p.add_argument("--portal-manifest", required=False)
    p.add_argument("--out", required=True)
    p.add_argument("--per-class", type=int, default=250)
    p.add_argument("--min-per-class", type=int, default=20)
    p.add_argument("--seed", type=int, default=20260819)
    p.add_argument("--timeout", type=int, default=180)
    return p.parse_args()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def normalize_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


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
    """Query BV-BRC with GET then POST fallback, preserving the literal RQL syntax."""
    suffix = (
        f"{expression}&select({','.join(fields)})&sort(%2Bid)"
        f"&limit({limit},{start})"
    )
    url = f"{API}/{collection}/?{suffix}"
    headers = {"Accept": "application/json", "User-Agent": UA}
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
                    "status": response.status_code,
                    "content_range": response.headers.get("Content-Range"),
                    "etag": response.headers.get("ETag"),
                    "date": response.headers.get("Date"),
                    "sha256": sha256_bytes(response.content),
                    "bytes": len(response.content),
                    "n_records": len(payload),
                }
                return payload, meta
            errors.append(f"GET {response.status_code}: {response.text[:500]}")
        except Exception as exc:  # network/transient parsing errors
            errors.append(f"GET exception: {exc!r}")

        try:
            body = f"{expression}&select({','.join(fields)})&sort(+id)&limit({limit},{start})"
            response = requests.post(
                f"{API}/{collection}/",
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
                    "url": f"{API}/{collection}/",
                    "query": body,
                    "status": response.status_code,
                    "content_range": response.headers.get("Content-Range"),
                    "etag": response.headers.get("ETag"),
                    "date": response.headers.get("Date"),
                    "sha256": sha256_bytes(response.content),
                    "bytes": len(response.content),
                    "n_records": len(payload),
                }
                return payload, meta
            errors.append(f"POST {response.status_code}: {response.text[:500]}")
        except Exception as exc:
            errors.append(f"POST exception: {exc!r}")
        if attempt < attempts:
            time.sleep(5 * attempt)
    raise RuntimeError("BV-BRC query failed after retries:\n" + "\n".join(errors))


def query_genomes(genome_ids: list[str], timeout: int) -> tuple[pd.DataFrame, list[dict]]:
    fields = [
        "genome_id", "genome_name", "taxon_id", "assembly_accession",
        "biosample_accession", "bioproject_accession", "country",
        "geographic_location", "collection_date", "host_name", "host_gender",
        "host_age", "isolation_source", "body_sample_site", "body_sample_subsite",
        "sequencing_status", "genome_status", "public", "owner", "date_inserted",
        "date_modified",
    ]
    records: list[dict] = []
    provenance: list[dict] = []
    for i in range(0, len(genome_ids), 100):
        batch = genome_ids[i : i + 100]
        values = ",".join(batch)
        payload, meta = request_rql(
            "genome", f"in(genome_id,({values}))", fields,
            limit=max(200, len(batch) + 10), timeout=timeout,
        )
        records.extend(payload)
        meta["batch_start"] = i
        meta["batch_n_requested"] = len(batch)
        provenance.append(meta)
    return pd.DataFrame(records), provenance


def canonical_phenotype(value: object) -> str | None:
    text = normalize_text(value).lower()
    if text in {"resistant", "non-susceptible", "nonsusceptible", "non susceptible"}:
        return "R"
    if text == "susceptible":
        return "S"
    return None


def first_nonempty(series: pd.Series) -> str:
    for value in series:
        text = normalize_text(value)
        if text:
            return text
    return ""


def stable_hash(value: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode()).hexdigest()


def select_diverse_balanced(df: pd.DataFrame, per_class: int, seed: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    work = df.copy()
    for c in ["bioproject_accession", "country", "pmid_text", "laboratory_typing_method"]:
        if c not in work:
            work[c] = ""
        work[c] = work[c].fillna("").astype(str)
    work["source_group"] = work["bioproject_accession"].where(
        work["bioproject_accession"].str.len() > 0,
        work["pmid_text"].where(work["pmid_text"].str.len() > 0, "UNKNOWN_SOURCE"),
    )
    work["diversity_key"] = (
        work["source_group"] + "|" + work["country"] + "|" +
        work["laboratory_typing_method"]
    )
    work["stable_hash"] = work["genome_id"].map(lambda x: stable_hash(str(x), seed))
    selected: list[pd.DataFrame] = []
    for phenotype in ["R", "S"]:
        sub = work[work["phenotype"].eq(phenotype)].copy()
        if sub.empty:
            continue
        # Round-robin across source/country/method groups; no genomic feature is inspected.
        sub = sub.sort_values(["diversity_key", "stable_hash", "genome_id"])
        sub["within_group_rank"] = sub.groupby("diversity_key").cumcount()
        sub = sub.sort_values(["within_group_rank", "stable_hash", "diversity_key"])
        selected.append(sub.head(per_class))
    return pd.concat(selected, ignore_index=True) if selected else work.iloc[0:0].copy()


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    prov = out / "provenance"
    prov.mkdir(exist_ok=True)

    portal_labels = pd.read_csv(args.portal_labels, dtype=str).fillna("")
    portal_manifest = (
        pd.read_csv(args.portal_manifest, dtype=str).fillna("")
        if args.portal_manifest and Path(args.portal_manifest).exists()
        else pd.DataFrame()
    )
    portal_assemblies = set(portal_labels.get("assembly_ID", pd.Series(dtype=str)).astype(str))
    portal_assembly_bases = {x.split(".")[0] for x in portal_assemblies if x}
    portal_biosamples = set(portal_labels.get("BioSample_ID", pd.Series(dtype=str)).astype(str))
    portal_bioprojects: set[str] = set()
    if not portal_manifest.empty:
        for col in ["BioProject", "BioProjects_all"]:
            if col in portal_manifest:
                for value in portal_manifest[col].astype(str):
                    portal_bioprojects.update(x for x in re.split(r"[;,| ]+", value) if x.startswith("PRJ"))

    amr_fields = [
        "id", "genome_id", "genome_name", "taxon_id", "antibiotic", "evidence",
        "resistant_phenotype", "measurement", "measurement_sign", "measurement_value",
        "measurement_unit", "laboratory_typing_method", "laboratory_typing_method_version",
        "laboratory_typing_platform", "testing_standard", "testing_standard_year", "vendor",
        "source", "pmid", "public", "date_inserted", "date_modified",
    ]
    amr_payload, amr_meta = request_rql(
        "genome_amr", "eq(antibiotic,colistin)", amr_fields,
        limit=250000, timeout=args.timeout,
    )
    (prov / "BVBRC_GENOME_AMR_QUERY.json").write_text(
        json.dumps(amr_meta, indent=2, ensure_ascii=False) + "\n"
    )
    raw_amr = pd.DataFrame(amr_payload)
    raw_amr.to_csv(out / "BVBRC_RAW_COLISTIN_AMR_RECORDS.csv", index=False)
    if raw_amr.empty:
        raise RuntimeError("BV-BRC returned no colistin AMR records")

    raw_amr["phenotype"] = raw_amr.get("resistant_phenotype", "").map(canonical_phenotype)
    raw_amr["genome_name_norm"] = raw_amr.get("genome_name", "").fillna("").astype(str).str.lower()
    raw_amr["taxon_id_text"] = raw_amr.get("taxon_id", "").fillna("").astype(str)
    kp = raw_amr[
        raw_amr["phenotype"].isin(["R", "S"])
        & (
            raw_amr["genome_name_norm"].str.contains(r"klebsiella\s+pneumoniae", regex=True)
            | raw_amr["taxon_id_text"].isin(["573", "72407", "1195464", "1463165", "1463164"])
        )
    ].copy()
    if kp.empty:
        raise RuntimeError("No K. pneumoniae-complex R/S colistin records after filtering")

    genome_ids = sorted(set(kp["genome_id"].dropna().astype(str)))
    genomes, genome_provenance = query_genomes(genome_ids, args.timeout)
    (prov / "BVBRC_GENOME_QUERY_BATCHES.json").write_text(
        json.dumps(genome_provenance, indent=2, ensure_ascii=False) + "\n"
    )
    genomes.to_csv(out / "BVBRC_MATCHED_GENOME_METADATA.csv", index=False)
    merged = kp.merge(genomes, on="genome_id", how="left", suffixes=("_amr", "_genome"))

    # Harmonize identifiers and provenance fields.
    def choose_col(frame: pd.DataFrame, names: list[str]) -> pd.Series:
        for name in names:
            if name in frame:
                return frame[name].fillna("").astype(str)
        return pd.Series([""] * len(frame), index=frame.index, dtype=str)

    merged["assembly_accession"] = choose_col(merged, ["assembly_accession", "assembly_accession_genome"])
    merged["biosample_accession"] = choose_col(merged, ["biosample_accession", "biosample_accession_genome"])
    merged["bioproject_accession"] = choose_col(merged, ["bioproject_accession", "bioproject_accession_genome"])
    merged["country"] = choose_col(merged, ["country", "country_genome", "geographic_location"])
    merged["collection_date"] = choose_col(merged, ["collection_date", "collection_date_genome"])
    merged["laboratory_typing_method"] = choose_col(merged, ["laboratory_typing_method"])
    merged["laboratory_typing_platform"] = choose_col(merged, ["laboratory_typing_platform"])
    merged["measurement_value"] = choose_col(merged, ["measurement_value"])
    merged["measurement_sign"] = choose_col(merged, ["measurement_sign"])
    merged["measurement_unit"] = choose_col(merged, ["measurement_unit"])
    merged["pmid_text"] = choose_col(merged, ["pmid"])
    merged["assembly_base"] = merged["assembly_accession"].str.split(".").str[0]
    merged["method_text"] = (
        merged["laboratory_typing_method"] + " " + merged["laboratory_typing_platform"]
    ).str.lower()
    merged["numeric_mic"] = pd.to_numeric(merged["measurement_value"], errors="coerce")
    merged["mic_like"] = (
        merged["numeric_mic"].notna()
        & merged["measurement_unit"].str.lower().str.contains(r"mg/l|ug/ml|µg/ml|μg/ml", regex=True)
    )
    merged["broth_like"] = merged["method_text"].str.contains(
        r"broth|microdilution|sensititre|mic", regex=True
    )
    merged["reference_method_like"] = merged["broth_like"] & merged["mic_like"]

    # Resolve duplicate tests at genome level without hiding R/S conflicts.
    conflict_rows: list[dict] = []
    resolved_rows: list[dict] = []
    for genome_id, group in merged.groupby("genome_id", dropna=False):
        states = sorted(set(group["phenotype"].dropna().astype(str)))
        if len(states) != 1:
            conflict_rows.append({
                "genome_id": genome_id,
                "phenotypes": ";".join(states),
                "n_records": len(group),
                "record_ids": ";".join(group.get("id", pd.Series(dtype=str)).astype(str)),
            })
            continue
        preferred = group.sort_values(
            ["reference_method_like", "mic_like", "numeric_mic"],
            ascending=[False, False, False],
        ).iloc[0]
        row = preferred.to_dict()
        row["n_amr_records"] = int(len(group))
        row["all_methods"] = ";".join(sorted(set(group["laboratory_typing_method"].astype(str))))
        row["all_platforms"] = ";".join(sorted(set(group["laboratory_typing_platform"].astype(str))))
        row["all_pmids"] = ";".join(sorted(set(group["pmid_text"].astype(str))))
        row["any_reference_method_like"] = bool(group["reference_method_like"].any())
        resolved_rows.append(row)
    conflicts = pd.DataFrame(conflict_rows)
    resolved = pd.DataFrame(resolved_rows)
    conflicts.to_csv(out / "BVBRC_GENOME_PHENOTYPE_CONFLICTS.csv", index=False)
    resolved.to_csv(out / "BVBRC_ALL_KP_COLISTIN_NONCONFLICT.csv", index=False)

    if resolved.empty:
        raise RuntimeError("No non-conflicting BV-BRC K. pneumoniae colistin genomes")
    resolved["overlap_assembly"] = resolved["assembly_base"].isin(portal_assembly_bases)
    resolved["overlap_biosample"] = resolved["biosample_accession"].isin(portal_biosamples) & resolved["biosample_accession"].ne("")
    resolved["overlap_bioproject"] = resolved["bioproject_accession"].isin(portal_bioprojects) & resolved["bioproject_accession"].ne("")
    resolved["assembly_biosample_disjoint"] = ~(resolved["overlap_assembly"] | resolved["overlap_biosample"])
    resolved["strict_bioproject_disjoint"] = resolved["assembly_biosample_disjoint"] & ~resolved["overlap_bioproject"]

    ab_disjoint = resolved[resolved["assembly_biosample_disjoint"]].copy()
    strict_disjoint = resolved[resolved["strict_bioproject_disjoint"]].copy()
    strict_method = strict_disjoint[strict_disjoint["any_reference_method_like"].eq(True)].copy()
    ab_disjoint.to_csv(out / "BVBRC_ASSEMBLY_BIOSAMPLE_DISJOINT.csv", index=False)
    strict_disjoint.to_csv(out / "BVBRC_BIOPROJECT_DISJOINT.csv", index=False)
    strict_method.to_csv(out / "BVBRC_BIOPROJECT_DISJOINT_REFERENCE_METHOD.csv", index=False)

    frozen = select_diverse_balanced(ab_disjoint, args.per_class, args.seed)
    strict_frozen = select_diverse_balanced(strict_disjoint, args.per_class, args.seed)
    strict_method_frozen = select_diverse_balanced(strict_method, args.per_class, args.seed)
    frozen.to_csv(out / "BVBRC_EXTERNAL_FROZEN_COHORT.csv", index=False)
    strict_frozen.to_csv(out / "BVBRC_EXTERNAL_STRICT_FROZEN_COHORT.csv", index=False)
    strict_method_frozen.to_csv(out / "BVBRC_EXTERNAL_REFERENCE_METHOD_FROZEN_COHORT.csv", index=False)
    for name, frame in [
        ("BVBRC_EXTERNAL_GENOME_IDS.txt", frozen),
        ("BVBRC_EXTERNAL_STRICT_GENOME_IDS.txt", strict_frozen),
        ("BVBRC_EXTERNAL_REFERENCE_METHOD_GENOME_IDS.txt", strict_method_frozen),
    ]:
        (out / name).write_text("\n".join(frame.get("genome_id", pd.Series(dtype=str)).astype(str)) + ("\n" if len(frame) else ""))

    def counts(frame: pd.DataFrame) -> dict:
        return frame.get("phenotype", pd.Series(dtype=str)).value_counts().sort_index().to_dict()

    summary = {
        "api": API,
        "amr_query": amr_meta,
        "raw_colistin_records": int(len(raw_amr)),
        "kp_rs_records": int(len(kp)),
        "unique_genome_ids_queried": int(len(genome_ids)),
        "genome_metadata_records": int(len(genomes)),
        "phenotype_conflict_genomes": int(len(conflicts)),
        "nonconflict_genomes": int(len(resolved)),
        "nonconflict_counts": counts(resolved),
        "with_ncbi_assembly": int(resolved["assembly_accession"].ne("").sum()),
        "assembly_biosample_disjoint": int(len(ab_disjoint)),
        "assembly_biosample_disjoint_counts": counts(ab_disjoint),
        "bioproject_disjoint": int(len(strict_disjoint)),
        "bioproject_disjoint_counts": counts(strict_disjoint),
        "bioproject_disjoint_reference_method": int(len(strict_method)),
        "bioproject_disjoint_reference_method_counts": counts(strict_method),
        "frozen": int(len(frozen)),
        "frozen_counts": counts(frozen),
        "strict_frozen": int(len(strict_frozen)),
        "strict_frozen_counts": counts(strict_frozen),
        "strict_method_frozen": int(len(strict_method_frozen)),
        "strict_method_frozen_counts": counts(strict_method_frozen),
        "external_validation_feasible": bool(
            counts(strict_frozen).get("R", 0) >= args.min_per_class
            and counts(strict_frozen).get("S", 0) >= args.min_per_class
        ),
        "reference_method_external_validation_feasible": bool(
            counts(strict_method_frozen).get("R", 0) >= args.min_per_class
            and counts(strict_method_frozen).get("S", 0) >= args.min_per_class
        ),
        "boundary": (
            "BV-BRC phenotypes are curated public records of heterogeneous provenance. "
            "Cohorts were formed before candidate testing. Assembly/BioSample and BioProject "
            "disjointness are reported separately; method metadata are incomplete."
        ),
    }
    (out / "BVBRC_EXTERNAL_COHORT_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str) + "\n"
    )
    report = [
        "# Candidate-blind BV-BRC external colistin cohort audit",
        "",
        f"- Raw public colistin AMR records: **{len(raw_amr):,}**",
        f"- K. pneumoniae-complex R/S records: **{len(kp):,}**",
        f"- Nonconflicting genome-level phenotypes: **{len(resolved):,}** {counts(resolved)}",
        f"- Assembly/BioSample-disjoint: **{len(ab_disjoint):,}** {counts(ab_disjoint)}",
        f"- BioProject-disjoint: **{len(strict_disjoint):,}** {counts(strict_disjoint)}",
        f"- BioProject-disjoint with broth/MIC-like metadata: **{len(strict_method):,}** {counts(strict_method)}",
        f"- External validation feasible at >= {args.min_per_class} per class: **{summary['external_validation_feasible']}**",
        f"- Reference-method subset feasible: **{summary['reference_method_external_validation_feasible']}**",
        "",
        "## Boundary",
        "",
        summary["boundary"],
    ]
    (out / "BVBRC_EXTERNAL_COHORT_REPORT.md").write_text("\n".join(report) + "\n")

    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
