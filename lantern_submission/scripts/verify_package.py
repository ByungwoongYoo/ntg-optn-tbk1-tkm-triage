#!/usr/bin/env python3
"""Verify ZIP structure and every entry listed in the internal SHA-256 manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import PurePosixPath, Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("zip_path", type=Path)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    if not a.zip_path.is_file() or a.zip_path.stat().st_size == 0:
        raise FileNotFoundError(a.zip_path)
    errors: list[str] = []
    with zipfile.ZipFile(a.zip_path) as archive:
        bad = archive.testzip()
        if bad:
            errors.append(f"corrupt member: {bad}")
        names = archive.namelist()
        if len(names) != len(set(names)):
            errors.append("duplicate ZIP entry")
        for name in names:
            path = PurePosixPath(name)
            if path.is_absolute() or ".." in path.parts:
                errors.append(f"unsafe ZIP path: {name}")
        required = {"SHA256SUMS.txt", "PACKAGE_MANIFEST.json"}
        missing_required = sorted(required - set(names))
        if missing_required:
            errors.append(f"missing required entries: {missing_required}")
        checked = 0
        if "SHA256SUMS.txt" in names:
            manifest = archive.read("SHA256SUMS.txt").decode("utf-8")
            for line_number, line in enumerate(manifest.splitlines(), 1):
                if not line.strip():
                    continue
                try:
                    expected, name = line.split("  ", 1)
                except ValueError:
                    errors.append(f"malformed SHA line {line_number}")
                    continue
                if name not in names:
                    errors.append(f"manifest entry missing from ZIP: {name}")
                    continue
                actual = hashlib.sha256(archive.read(name)).hexdigest()
                if actual != expected:
                    errors.append(f"checksum mismatch: {name}")
                checked += 1
        package_manifest = {}
        if "PACKAGE_MANIFEST.json" in names:
            package_manifest = json.loads(archive.read("PACKAGE_MANIFEST.json"))
            if package_manifest.get("active_reads_included") is not False:
                errors.append("package manifest does not exclude active reads")
            if package_manifest.get("derived_active_assemblies_included") is not False:
                errors.append("package manifest does not exclude active assemblies")
    result = {
        "status": "PASS" if not errors else "FAIL",
        "zip_path": str(a.zip_path),
        "zip_bytes": a.zip_path.stat().st_size,
        "zip_sha256": hashlib.sha256(a.zip_path.read_bytes()).hexdigest(),
        "entries": len(names),
        "manifest_entries_checked": checked,
        "errors": errors,
        "package_manifest": package_manifest,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
