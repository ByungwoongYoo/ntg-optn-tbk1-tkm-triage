#!/usr/bin/env python3
"""Run one candidate against EMBL-EBI Job Dispatcher NCBI BLAST.

Protein queries use UniProtKB; nucleotide contigs use the current ENA sequence
release (em_rel). The raw submission, status, result-type, XML and text outputs
are preserved. No-hit and service-failure states are distinguished explicitly.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests
from Bio import SeqIO

BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_record(path: Path, query_id: str) -> str:
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id == query_id:
            return f">{record.id}\n{str(record.seq)}\n"
    raise KeyError(f"{query_id!r} absent from {path}")


def submit(session: requests.Session, params: dict[str, Any]) -> str:
    last: Exception | None = None
    for attempt in range(1, 6):
        try:
            r = session.post(f"{BASE}/run", data=params, timeout=120)
            r.raise_for_status()
            job = r.text.strip()
            if not job or "<" in job:
                raise RuntimeError(f"unexpected job id: {job[:200]}")
            return job
        except Exception as exc:
            last = exc
            if attempt < 5:
                time.sleep(min(60, 5 * 2 ** attempt))
    assert last is not None
    raise last


def poll(session: requests.Session, job: str, max_minutes: int) -> list[dict[str, Any]]:
    history: list[dict[str, Any]] = []
    deadline = time.time() + max_minutes * 60
    while time.time() < deadline:
        r = session.get(f"{BASE}/status/{job}", timeout=120)
        if r.status_code in (429, 502, 503, 504):
            history.append({"utc": utc(), "http": r.status_code, "status": "transient"})
            time.sleep(15)
            continue
        r.raise_for_status()
        status = r.text.strip().upper()
        history.append({"utc": utc(), "http": r.status_code, "status": status})
        if status == "FINISHED":
            return history
        if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"EBI BLAST job {job} ended as {status}")
        time.sleep(10)
    raise TimeoutError(f"EBI BLAST job {job} did not finish in {max_minutes} minutes")


def get_result(session: requests.Session, job: str, result_type: str) -> requests.Response:
    r = session.get(f"{BASE}/result/{job}/{result_type}", timeout=300)
    r.raise_for_status()
    return r


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_xml(xml_text: str, query: str, program: str, database: str) -> list[dict[str, Any]]:
    root = ET.fromstring(xml_text)
    rows: list[dict[str, Any]] = []
    rank = 0
    for hit in root.iter():
        if local_name(hit.tag) != "hit":
            continue
        rank += 1
        attrs = hit.attrib
        access = attrs.get("ac", "")
        hit_id = attrs.get("id", "")
        desc = attrs.get("description", "")
        hit_len = attrs.get("length", "")
        alignments = [x for x in hit.iter() if local_name(x.tag) == "alignment"]
        if not alignments:
            rows.append({
                "query": query, "program": program, "database": database,
                "rank": rank, "accession": access, "hit_id": hit_id,
                "description": desc, "hit_length": hit_len,
            })
            continue
        for aln_index, aln in enumerate(alignments, start=1):
            values = {local_name(x.tag): (x.text or "") for x in aln.iter() if x is not aln}
            qstart = values.get("queryStart", "")
            qend = values.get("queryEnd", "")
            qlen = values.get("queryLength", values.get("queryLen", ""))
            qcov = ""
            try:
                qcov = 100.0 * (abs(int(qend) - int(qstart)) + 1) / int(qlen)
            except Exception:
                pass
            rows.append({
                "query": query,
                "program": program,
                "database": database,
                "rank": rank,
                "alignment_index": aln_index,
                "accession": access,
                "hit_id": hit_id,
                "description": desc,
                "hit_length": hit_len,
                "score": values.get("score", ""),
                "bits": values.get("bits", ""),
                "evalue": values.get("expectation", ""),
                "identity_percent": values.get("identity", ""),
                "positives_percent": values.get("positives", ""),
                "gaps": values.get("gaps", ""),
                "query_start": qstart,
                "query_end": qend,
                "query_length": qlen,
                "query_coverage_percent": qcov,
                "match_start": values.get("matchStart", ""),
                "match_end": values.get("matchEnd", ""),
                "query_seq": values.get("querySeq", ""),
                "match_seq": values.get("matchSeq", ""),
            })
    return rows


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--query-id", required=True)
    p.add_argument("--program", choices=["blastp", "blastn"], required=True)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--database", required=True)
    p.add_argument("--email", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--max-minutes", type=int, default=55)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "ByungwoongYoo-Panax-RdRP-EBI-BLAST-audit/0.1",
        "Accept": "text/plain,application/json,application/xml",
    })
    sequence = read_record(args.fasta, args.query_id)
    params: dict[str, Any] = {
        "email": args.email,
        "program": args.program,
        "stype": "protein" if args.program == "blastp" else "dna",
        "sequence": sequence,
        "database": args.database,
        "alignments": "100",
        "scores": "100",
        "exp": "0.00001",
        "filter": "T",
        "gapalign": "true",
    }
    if args.program == "blastp":
        params["matrix"] = "BLOSUM62"
    else:
        params["match_scores"] = "2,-3"

    metadata: dict[str, Any] = {
        "query": args.query_id,
        "program": args.program,
        "database": args.database,
        "submitted_utc": utc(),
        "parameters": {k: v for k, v in params.items() if k not in {"sequence", "email"}},
        "input_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
    }
    try:
        for endpoint in ("parameters", "parameterdetails/database"):
            r = session.get(f"{BASE}/{endpoint}", headers={"Accept": "application/json"}, timeout=120)
            if r.ok:
                (args.out / endpoint.replace("/", "_")).with_suffix(".json").write_text(r.text, encoding="utf-8")
        job = submit(session, params)
        metadata["job_id"] = job
        (args.out / "JOB_ID.txt").write_text(job + "\n", encoding="utf-8")
        history = poll(session, job, args.max_minutes)
        (args.out / "STATUS_HISTORY.json").write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
        rt = session.get(f"{BASE}/resulttypes/{job}", headers={"Accept": "application/json"}, timeout=120)
        rt.raise_for_status()
        (args.out / "RESULT_TYPES.json").write_text(rt.text, encoding="utf-8")

        xml_r = get_result(session, job, "xml")
        (args.out / "result.xml").write_text(xml_r.text, encoding="utf-8")
        try:
            out_r = get_result(session, job, "out")
            (args.out / "result.out.txt").write_text(out_r.text, encoding="utf-8")
        except Exception as exc:
            (args.out / "result_out_error.txt").write_text(repr(exc) + "\n", encoding="utf-8")

        rows = parse_xml(xml_r.text, args.query_id, args.program, args.database)
        fields = [
            "query", "program", "database", "rank", "alignment_index", "accession",
            "hit_id", "description", "hit_length", "score", "bits", "evalue",
            "identity_percent", "positives_percent", "gaps", "query_start", "query_end",
            "query_length", "query_coverage_percent", "match_start", "match_end",
            "query_seq", "match_seq",
        ]
        with (args.out / "hits.tsv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        metadata.update({"success": True, "completed_utc": utc(), "parsed_hits": len(rows)})
    except Exception as exc:
        metadata.update({"success": False, "completed_utc": utc(), "error": repr(exc)})
    (args.out / "STATUS.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    manifest = []
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            manifest.append(f"{sha256(path)}  {path.name}")
    (args.out / "SHA256SUMS.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0 if metadata.get("success") else 1


if __name__ == "__main__":
    raise SystemExit(main())
