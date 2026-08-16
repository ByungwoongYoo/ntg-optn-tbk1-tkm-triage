#!/usr/bin/env python3
from __future__ import annotations

import gzip
import io
import json
import os
from pathlib import Path

import pandas as pd
import requests
from remotezip import RemoteZip

OUT = Path(os.environ.get("OUT_DIR", "artifact/cpt_sample_probe"))
OUT.mkdir(parents=True, exist_ok=True)
record = requests.get("https://zenodo.org/api/records/7954657", timeout=120).json()
results = []
for file in record.get("files", []):
    name = file["key"]
    if not name.endswith(".zip"):
        continue
    url = file.get("links", {}).get("content") or file.get("links", {}).get("self")
    rz = RemoteZip(url)
    members = [n for n in rz.namelist() if n.endswith(".csv.gz")]
    selected = []
    for target in ["ACADM_HUMAN", "BRCA1_HUMAN", "PTEN_HUMAN", "APOE_HUMAN"]:
        match = next((m for m in members if f"/{target}.csv.gz" in m), None)
        if match:
            selected.append(match)
    if not selected and members:
        selected = members[:1]
    archive_result = {"archive": name, "selected": []}
    for member in selected[:3]:
        raw = rz.read(member)
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
            frame = pd.read_csv(gz, nrows=10)
        archive_result["selected"].append(
            {
                "member": member,
                "bytes_compressed": len(raw),
                "columns": frame.columns.tolist(),
                "preview": frame.head(5).to_dict(orient="records"),
            }
        )
    rz.close()
    results.append(archive_result)
(OUT / "sample_probe.json").write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
