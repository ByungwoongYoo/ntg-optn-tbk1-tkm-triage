#!/usr/bin/env python3
"""Download exact public BV-BRC genome assemblies for a frozen cohort.

Each genome is queried by its already-frozen BV-BRC genome_id. The script writes one
FASTA per genome, records response hashes and sequence-quality statistics, and never
uses candidate sequence information. The resulting genomes are public database copies;
they are removed before the final workflow artifact is uploaded.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import re
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://www.bv-brc.org/api/genome_sequence/"
UA = "ByungwoongYoo-AMR-external-validation/20260819"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--attempts", type=int, default=4)
    p.add_argument("--min-genome", type=int, default=4_500_000)
    p.add_argument("--max-genome", type=int, default=7_500_000)
    p.add_argument("--max-contigs", type=int, default=1000)
    p.add_argument("--min-n50", type=int, default=5000)
    p.add_argument("--max-n-fraction", type=float, default=0.02)
    p.add_argument("--min-per-class", type=int, default=50)
    return p.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def fasta_stats(data: bytes) -> dict:
    lengths: list[int] = []
    n_bases = 0
    total = 0
    current = 0
    saw_header = False
    for raw in data.decode("utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            saw_header = True
            if current:
                lengths.append(current)
                current = 0
            continue
        seq = re.sub(r"\s+", "", line).upper()
        if not re.fullmatch(r"[ACGTRYSWKMBDHVN.-]+", seq):
            raise ValueError("Non-FASTA sequence payload received")
        current += len(seq)
        total += len(seq)
        n_bases += seq.count("N")
    if current:
        lengths.append(current)
    if not saw_header or not lengths:
        raise ValueError("Empty or invalid FASTA payload")
    lengths.sort(reverse=True)
    half = total / 2
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half:
            n50 = length
            break
    return {
        "genome_length": int(total),
        "contigs": int(len(lengths)),
        "n50": int(n50),
        "n_fraction": float(n_bases / total) if total else 1.0,
        "largest_contig": int(lengths[0]),
    }


def get_one(genome_id: str, out_dir: Path, timeout: int, attempts: int) -> dict:
    expression = f"eq(genome_id,{genome_id})"
    url = (
        f"{API}?{expression}&sort(%2Bsequence_id)&limit(25000,0)"
        "&http_accept=application/dna%2Bfasta"
    )
    headers = {"Accept": "application/dna+fasta", "User-Agent": UA}
    errors: list[str] = []
    for attempt in range(1, attempts + 1):
        try:
            response = requests.get(url, headers=headers, timeout=timeout)
            if response.status_code == 200:
                content = response.content
                stats = fasta_stats(content)
                path = out_dir / f"{safe_name(genome_id)}.fna"
                path.write_bytes(content)
                return {
                    "genome_id": genome_id,
                    "status": "downloaded",
                    "path": str(path),
                    "bytes": len(content),
                    "sha256": hashlib.sha256(content).hexdigest(),
                    "content_range": response.headers.get("Content-Range"),
                    "etag": response.headers.get("ETag"),
                    "date": response.headers.get("Date"),
                    **stats,
                }
            errors.append(f"HTTP {response.status_code}: {response.text[:300]}")
        except Exception as exc:
            errors.append(repr(exc))
        if attempt < attempts:
            time.sleep(5 * attempt)
    return {"genome_id": genome_id, "status": "failed", "path": "", "error": " | ".join(errors)[-4000:]}


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    assemblies = out / "assemblies"
    assemblies.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(args.cohort, dtype=str).fillna("")
    if not {"genome_id", "phenotype"}.issubset(cohort.columns):
        raise ValueError("Cohort must contain genome_id and phenotype")
    cohort = cohort[cohort["phenotype"].isin(["R", "S"])].drop_duplicates("genome_id")
    genome_ids = sorted(cohort["genome_id"].astype(str))
    if not genome_ids:
        raise ValueError("Frozen cohort is empty")

    results: list[dict] = []
    with cf.ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        futures = {pool.submit(get_one, gid, assemblies, args.timeout, args.attempts): gid for gid in genome_ids}
        for i, future in enumerate(cf.as_completed(futures), start=1):
            results.append(future.result())
            if i % 25 == 0 or i == len(futures):
                print(f"downloaded_or_attempted={i}/{len(futures)}", flush=True)

    table = pd.DataFrame(results)
    for column in ["genome_length", "contigs", "n50", "n_fraction", "largest_contig"]:
        if column not in table:
            table[column] = math.nan
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["qc_pass"] = (
        table["status"].eq("downloaded")
        & table["genome_length"].between(args.min_genome, args.max_genome, inclusive="both")
        & table["contigs"].le(args.max_contigs)
        & table["n50"].ge(args.min_n50)
        & table["n_fraction"].le(args.max_n_fraction)
    )
    table = cohort.merge(table, on="genome_id", how="left", validate="one_to_one")
    table.to_csv(out / "BVBRC_GENOME_DOWNLOAD_AND_QC.csv", index=False)
    passed = table[table["qc_pass"].eq(True)].copy()
    passed.to_csv(out / "BVBRC_QC_PASSED_COHORT.csv", index=False)
    (out / "BVBRC_QC_PASSED_GENOME_IDS.txt").write_text(
        "\n".join(passed["genome_id"].astype(str)) + ("\n" if len(passed) else "")
    )
    refs = [str((assemblies / f"{safe_name(gid)}.fna").resolve()) for gid in passed["genome_id"]]
    (out / "BVBRC_QC_PASSED_REFS.txt").write_text("\n".join(refs) + ("\n" if refs else ""))

    counts = passed["phenotype"].value_counts().sort_index().to_dict()
    summary = {
        "n_frozen": int(len(cohort)),
        "frozen_counts": cohort["phenotype"].value_counts().sort_index().to_dict(),
        "n_downloaded": int(table["status"].eq("downloaded").sum()),
        "n_failed": int(table["status"].eq("failed").sum()),
        "n_qc_pass": int(len(passed)),
        "qc_pass_counts": counts,
        "qc_thresholds": {
            "min_genome": args.min_genome,
            "max_genome": args.max_genome,
            "max_contigs": args.max_contigs,
            "min_n50": args.min_n50,
            "max_n_fraction": args.max_n_fraction,
        },
        "min_per_class": args.min_per_class,
        "validation_size_gate": bool(
            counts.get("R", 0) >= args.min_per_class and counts.get("S", 0) >= args.min_per_class
        ),
        "boundary": "Genome quality filters and the minimum class size were fixed before candidate sequence testing. Only the laboratory-method frozen cohort is eligible for validation.",
    }
    (out / "BVBRC_GENOME_DOWNLOAD_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt" and "assemblies" not in path.parts:
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out)}")
    (out / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if not summary["validation_size_gate"]:
        raise SystemExit(f"Fewer than {args.min_per_class} QC-passing genomes per phenotype")


if __name__ == "__main__":
    main()
