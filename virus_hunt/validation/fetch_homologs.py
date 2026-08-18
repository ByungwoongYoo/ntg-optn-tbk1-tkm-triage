#!/usr/bin/env python3
"""Fetch top current viral BLASTP homologs and construct small audit trees."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

USER_AGENT = "public-virome-hunt/0.2 (independent reproducibility audit)"


def fetch_fasta(accessions: list[str]) -> str:
    params = urllib.parse.urlencode({
        "db": "protein",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    })
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    last: Exception | None = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"NCBI EFetch failed: {last}")


def safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--blast", type=Path, required=True)
    p.add_argument("--queries", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--top", type=int, default=15)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    query_records = {r.id: r for r in SeqIO.parse(args.queries, "fasta")}
    by_query: dict[str, list[str]] = defaultdict(list)
    if args.blast.exists():
        with args.blast.open(encoding="utf-8", errors="replace") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 2:
                    continue
                q, acc = row[0], row[1]
                if acc not in by_query[q] and len(by_query[q]) < args.top:
                    by_query[q].append(acc)

    status_lines = ["query\thomologs_requested\thomologs_retrieved\tmafft\tfasttree"]
    for q, qrec in query_records.items():
        qdir = args.out_dir / safe(q)
        qdir.mkdir(exist_ok=True)
        accessions = by_query.get(q, [])
        fasta_text = ""
        if accessions:
            fasta_text = fetch_fasta(accessions)
            time.sleep(0.4)
        homolog_path = qdir / "homologs.faa"
        homolog_path.write_text(fasta_text, encoding="utf-8")
        combined = qdir / "candidate_plus_homologs.faa"
        with combined.open("w", encoding="utf-8") as out:
            out.write(f">{q}|CANDIDATE\n{str(qrec.seq)}\n")
            out.write(fasta_text)
        retrieved = sum(1 for _ in SeqIO.parse(homolog_path, "fasta"))
        mafft_status = "not_run"
        tree_status = "not_run"
        if retrieved >= 2:
            aln = qdir / "alignment.faa"
            tree = qdir / "tree.nwk"
            with aln.open("w") as out:
                proc = subprocess.run(
                    ["mafft", "--auto", "--quiet", str(combined)],
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            mafft_status = f"exit_{proc.returncode}"
            (qdir / "mafft.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
            if proc.returncode == 0:
                with tree.open("w") as out:
                    proc2 = subprocess.run(
                        ["FastTree", "-wag", str(aln)],
                        stdout=out,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                tree_status = f"exit_{proc2.returncode}"
                (qdir / "fasttree.stderr.txt").write_text(proc2.stderr or "", encoding="utf-8")
        status_lines.append(f"{q}\t{len(accessions)}\t{retrieved}\t{mafft_status}\t{tree_status}")
    (args.out_dir / "phylogeny_status.tsv").write_text("\n".join(status_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
