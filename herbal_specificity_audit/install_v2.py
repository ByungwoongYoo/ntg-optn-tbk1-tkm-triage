#!/usr/bin/env python3
"""Reconstruct the frozen v2 audit source from repository-safe compressed parts."""
from pathlib import Path
import base64
import gzip
import hashlib

root = Path(__file__).resolve().parent
parts = sorted((root / "v2_parts").glob("part*.b64"))
if len(parts) != 6:
    raise RuntimeError(f"Expected 6 source parts, found {len(parts)}")
payload = "".join(p.read_text(encoding="utf-8").strip() for p in parts)
source = gzip.decompress(base64.b64decode(payload))
out = root / "audit.py"
out.write_bytes(source)
print(f"reconstructed={out} bytes={len(source)} sha256={hashlib.sha256(source).hexdigest()}")
