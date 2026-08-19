#!/usr/bin/env python3
"""Freeze an external NCBI Pathogen Detection K. pneumoniae-colistin cohort.

The current official Klebsiella metadata TSV is downloaded and hashed. Submitter-supplied
colistin R/S calls are extracted without inspecting genomic candidates. Assemblies and
BioSamples already present in the EMBL-EBI AMR Portal cohort are removed. A stricter
manifest additionally removes overlapping BioProjects.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from collections import deque
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests

INDEXES = [
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Metadata/",
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_kmer/Metadata/",
]
ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")
BIOPROJECT_RE = re.compile(r"PRJ(?:NA|EB|DB)\d+", re.I)
BIOSAMPLE_RE = re.compile(r"(?:SAMN|SAMEA|SAMD)\d+", re.I)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--portal-labels", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-per-class", type=int, default=250)
    p.add_argument("--strict-max-per-class", type=int, default=150)
    return p.parse_args()


def get_text(url: str, timeout: int = 60) -> str:
    response = requests.get(url, timeout=timeout, headers={"User-Agent": "AMR-external-cohort-audit/1.0"})
    response.raise_for_status()
    return response.text


def discover_metadata() -> tuple[dict, list[dict]]:
    found: list[dict] = []
    for index in INDEXES:
        try:
            text = get_text(index)
        except Exception as exc:
            found.append({"index": index, "error": repr(exc)})
            continue
        for href in re.findall(r"href=[\"']([^\"']+\.metadata\.tsv)[\"']", text, re.I):
            name = html.unescape(href)
            match = re.search(r"PDG\d+\.(\d+)\.metadata\.tsv$", name)
            found.append({
                "index": index,
                "url": urljoin(index, name),
                "name": Path(name).name,
                "version": int(match.group(1)) if match else -1,
            })
    viable = [item for item in found if item.get("url")]
    if not viable:
        raise RuntimeError(f"No metadata TSV discovered: {found}")
    viable.sort(key=lambda item: (item["version"], "latest_snps" in item["index"]), reverse=True)
    return viable[0], found


def download(url: str, path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with requests.get(
        url,
        stream=True,
        timeout=(60, 600),
        headers={"User-Agent": "AMR-external-cohort-audit/1.0"},
    ) as response:
        response.raise_for_status()
        with open(path, "wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
    return digest.hexdigest(), size


def normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def choose_column(columns: list[str], *aliases: str) -> str | None:
    normalized = {normalize(column): column for column in columns}
    for alias in aliases:
        if normalize(alias) in normalized:
            return normalized[normalize(alias)]
    return None


def extract_label(value: str) -> str | None:
    text = str(value)
    resistant = bool(re.search(r"(?i)(?:^|[^a-z0-9])colistin\s*=\s*R(?:$|[^a-z0-9])", text))
    susceptible = bool(re.search(r"(?i)(?:^|[^a-z0-9])colistin\s*=\s*S(?:$|[^a-z0-9])", text))
    if resistant and susceptible:
        return "CONFLICT"
    if resistant:
        return "R"
    if susceptible:
        return "S"
    return None


def extract_tokens(value: str, pattern: re.Pattern) -> list[str]:
    return sorted(set(match.upper() for match in pattern.findall(str(value))))


def round_robin(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    group_key = "BioProject" if frame.BioProject.ne("").any() else "country"
    queues = []
    for _, group in frame.sort_values([group_key, "assembly_ID"]).groupby(group_key, dropna=False, sort=True):
        queues.append(deque(group.to_dict("records")))
    output = []
    while queues and len(output) < n:
        next_queues = []
        for queue in queues:
            if queue and len(output) < n:
                output.append(queue.popleft())
            if queue:
                next_queues.append(queue)
        queues = next_queues
    return pd.DataFrame(output)


def freeze_balanced(frame: pd.DataFrame, cap: int) -> pd.DataFrame:
    counts = frame.phenotype.value_counts()
    n = min(cap, int(counts.min())) if {"R", "S"}.issubset(set(counts.index)) else 0
    if n == 0:
        return pd.DataFrame(columns=frame.columns)
    return pd.concat(
        [round_robin(frame[frame.phenotype.eq(phenotype)], n) for phenotype in ["R", "S"]],
        ignore_index=True,
    )


def main() -> None:
    a = parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    selected, discovered = discover_metadata()
    (out / "DISCOVERED_METADATA_CANDIDATES.json").write_text(json.dumps(discovered, indent=2) + "\n")
    raw = out / "ncbi_klebsiella_metadata.tsv"
    metadata_sha, metadata_size = download(selected["url"], raw)
    with open(raw, "r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    (out / "METADATA_HEADER.json").write_text(json.dumps(header, indent=2) + "\n")

    assembly_column = choose_column(header, "asm_acc", "assembly", "assembly accession")
    biosample_column = choose_column(header, "biosample_acc", "biosample")
    bioproject_column = choose_column(header, "bioproject_acc", "bioproject")
    scientific_name_column = choose_column(header, "scientific_name", "organism", "species")
    ast_column = choose_column(header, "AST_phenotypes", "ast phenotypes")
    genotype_column = choose_column(header, "AMR_genotypes", "amr genotypes")
    country_column = choose_column(header, "geo_loc_name", "country", "location")
    collection_column = choose_column(header, "collection_date", "collection date")
    source_column = choose_column(header, "isolation_source", "isolation source")

    required = {
        "assembly": assembly_column,
        "scientific_name": scientific_name_column,
        "AST_phenotypes": ast_column,
    }
    if any(value is None for value in required.values()):
        raise RuntimeError(f"Missing required columns {required}; header={header}")
    use_columns = [
        value
        for value in [
            assembly_column,
            biosample_column,
            bioproject_column,
            scientific_name_column,
            ast_column,
            genotype_column,
            country_column,
            collection_column,
            source_column,
        ]
        if value is not None
    ]

    records = []
    parser_stats = {
        "rows_read": 0,
        "kp_rows": 0,
        "colistin_rows": 0,
        "conflicts": 0,
        "missing_assembly": 0,
    }
    for chunk in pd.read_csv(
        raw,
        sep="\t",
        dtype=str,
        usecols=use_columns,
        chunksize=200000,
        keep_default_na=False,
        quoting=csv.QUOTE_MINIMAL,
        low_memory=False,
    ):
        parser_stats["rows_read"] += len(chunk)
        scientific_name = chunk[scientific_name_column].astype(str)
        kp_mask = (
            scientific_name.str.match(r"(?i)^Klebsiella pneumoniae(?:$|\s)")
            & ~scientific_name.str.contains("quasipneumoniae|variicola", case=False, regex=True)
        )
        subset = chunk[kp_mask].copy()
        parser_stats["kp_rows"] += len(subset)
        subset["phenotype"] = subset[ast_column].map(extract_label)
        parser_stats["conflicts"] += int(subset.phenotype.eq("CONFLICT").sum())
        subset = subset[subset.phenotype.isin(["R", "S"])].copy()
        parser_stats["colistin_rows"] += len(subset)
        for _, row in subset.iterrows():
            accession_match = ACC_RE.search(str(row[assembly_column]))
            if not accession_match:
                parser_stats["missing_assembly"] += 1
                continue
            biosamples = extract_tokens(row[biosample_column], BIOSAMPLE_RE) if biosample_column else []
            bioprojects = extract_tokens(row[bioproject_column], BIOPROJECT_RE) if bioproject_column else []
            records.append({
                "assembly_ID": accession_match.group(0),
                "BioSample_ID": biosamples[0] if biosamples else "",
                "BioProjects_all": ";".join(bioprojects),
                "BioProject": bioprojects[0] if bioprojects else "",
                "scientific_name": row[scientific_name_column],
                "phenotype": row["phenotype"],
                "country": row[country_column] if country_column else "",
                "collection_date": row[collection_column] if collection_column else "",
                "isolation_source": row[source_column] if source_column else "",
                "AMR_genotypes": row[genotype_column] if genotype_column else "",
                "AST_phenotypes": row[ast_column],
            })

    all_candidates = pd.DataFrame(records)
    if all_candidates.empty:
        raise RuntimeError("No K. pneumoniae colistin R/S rows parsed")
    all_candidates = all_candidates.drop_duplicates()
    conflict_counts = all_candidates.groupby("assembly_ID").phenotype.nunique()
    conflict_accessions = set(conflict_counts[conflict_counts > 1].index)
    all_candidates[all_candidates.assembly_ID.isin(conflict_accessions)].to_csv(
        out / "NCBI_ASSEMBLY_PHENOTYPE_CONFLICTS.csv", index=False
    )
    all_candidates = all_candidates[~all_candidates.assembly_ID.isin(conflict_accessions)].drop_duplicates("assembly_ID")

    portal = pd.read_csv(a.portal_labels, dtype=str, keep_default_na=False)
    portal_assemblies = set(portal.assembly_ID.astype(str))
    portal_biosamples = set(portal.get("BioSample_ID", pd.Series(dtype=str)).astype(str))
    portal_projects: set[str] = set()
    for value in portal.get("BioProject", pd.Series(dtype=str)).astype(str):
        portal_projects.update(extract_tokens(value, BIOPROJECT_RE))
    for value in portal.get("BioProjects_all", pd.Series(dtype=str)).astype(str):
        portal_projects.update(extract_tokens(value, BIOPROJECT_RE))

    all_candidates["assembly_overlap"] = all_candidates.assembly_ID.isin(portal_assemblies)
    all_candidates["biosample_overlap"] = (
        all_candidates.BioSample_ID.isin(portal_biosamples) & all_candidates.BioSample_ID.ne("")
    )
    all_candidates["bioproject_overlap"] = all_candidates.BioProjects_all.map(
        lambda value: bool(set(value.split(";")) & portal_projects) if value else False
    )
    all_candidates["has_mcr_catalog_call"] = all_candidates.AMR_genotypes.str.contains(
        r"(?i)(^|[^a-z0-9])mcr[-_]?\d", regex=True, na=False
    )
    all_candidates.to_csv(out / "NCBI_ALL_COLISTIN_RS_ASSEMBLIES.csv", index=False)

    assembly_external = all_candidates[
        ~all_candidates.assembly_overlap & ~all_candidates.biosample_overlap
    ].copy()
    strict_external = assembly_external[~assembly_external.bioproject_overlap].copy()
    assembly_external.to_csv(out / "NCBI_ASSEMBLY_BIOSAMPLE_DISJOINT.csv", index=False)
    strict_external.to_csv(out / "NCBI_BIOPROJECT_DISJOINT.csv", index=False)

    frozen = freeze_balanced(assembly_external, a.max_per_class)
    strict_frozen = freeze_balanced(strict_external, a.strict_max_per_class)
    frozen.to_csv(out / "NCBI_EXTERNAL_FROZEN_COHORT.csv", index=False)
    strict_frozen.to_csv(out / "NCBI_EXTERNAL_STRICT_FROZEN_COHORT.csv", index=False)
    (out / "NCBI_EXTERNAL_ACCESSIONS.txt").write_text(
        "\n".join(frozen.assembly_ID.astype(str)) + "\n" if len(frozen) else ""
    )
    (out / "NCBI_EXTERNAL_STRICT_ACCESSIONS.txt").write_text(
        "\n".join(strict_frozen.assembly_ID.astype(str)) + "\n" if len(strict_frozen) else ""
    )

    summary = {
        "metadata_url": selected["url"],
        "metadata_name": selected["name"],
        "metadata_version": selected["version"],
        "metadata_sha256": metadata_sha,
        "metadata_size_bytes": metadata_size,
        "parser_stats": parser_stats,
        "all_unique_nonconflict": len(all_candidates),
        "all_counts": all_candidates.phenotype.value_counts().to_dict(),
        "assembly_biosample_disjoint": len(assembly_external),
        "assembly_disjoint_counts": assembly_external.phenotype.value_counts().to_dict(),
        "bioproject_disjoint": len(strict_external),
        "bioproject_disjoint_counts": strict_external.phenotype.value_counts().to_dict(),
        "frozen": len(frozen),
        "frozen_counts": frozen.phenotype.value_counts().to_dict() if len(frozen) else {},
        "strict_frozen": len(strict_frozen),
        "strict_frozen_counts": strict_frozen.phenotype.value_counts().to_dict() if len(strict_frozen) else {},
        "mcr_catalog_calls_in_frozen": int(frozen.has_mcr_catalog_call.sum()) if len(frozen) else 0,
        "external_validation_feasible": bool(
            len(frozen) >= 60 and frozen.phenotype.value_counts().min() >= 30
        ),
        "strict_external_validation_feasible": bool(
            len(strict_frozen) >= 40 and strict_frozen.phenotype.value_counts().min() >= 20
        ),
        "boundary": "NCBI AST phenotypes are submitter-supplied calls. Assembly/BioSample-disjoint and BioProject-disjoint sets are reported separately. Cohort formation is candidate-blind.",
    }
    (out / "NCBI_EXTERNAL_COHORT_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    raw.unlink()
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
