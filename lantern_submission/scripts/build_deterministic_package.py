#!/usr/bin/env python3
"""Build a deterministic source package with an internal SHA-256 manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import zipfile
from pathlib import Path

FIXED_TIME = (1980, 1, 1, 0, 0, 0)
EXCLUDED_PARTS = {".git", "__pycache__", ".pytest_cache", ".mypy_cache"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--include", action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    return parser.parse_args()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def zip_info(name: str, executable: bool = False) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, FIXED_TIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    mode = 0o755 if executable else 0o644
    info.external_attr = (stat.S_IFREG | mode) << 16
    info.create_system = 3
    return info


def is_executable(path: Path) -> bool:
    return path.suffix in {".sh", ".py"} or bool(path.stat().st_mode & stat.S_IXUSR)


def main() -> None:
    args = parse_args()
    root = args.source_root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    files: dict[str, bytes] = {}
    for raw in args.include:
        item = (root / raw).resolve()
        try:
            item.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"include escapes source root: {raw}") from exc
        if not item.exists():
            raise FileNotFoundError(item)
        candidates = [item] if item.is_file() else sorted(path for path in item.rglob("*") if path.is_file())
        for path in candidates:
            relative = path.relative_to(root)
            if any(part in EXCLUDED_PARTS for part in relative.parts):
                continue
            name = relative.as_posix()
            if name in files:
                raise ValueError(f"duplicate archive path: {name}")
            files[name] = path.read_bytes()
    if not files:
        raise ValueError("no files selected")

    checksums = "".join(f"{sha256_bytes(files[name])}  {name}\n" for name in sorted(files))
    package_manifest = {
        "status": "DETERMINISTIC_SUBMISSION_PACKAGE",
        "source_commit": args.source_commit,
        "files_hashed": len(files),
        "active_reads_included": False,
        "derived_active_assemblies_included": False,
        "public_release_authorized": False,
    }
    generated = {
        "SHA256SUMS.txt": checksums.encode("utf-8"),
        "PACKAGE_MANIFEST.json": (json.dumps(package_manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9, allowZip64=True) as archive:
        for name in sorted({**files, **generated}):
            data = generated[name] if name in generated else files[name]
            source_path = root / name
            executable = source_path.exists() and is_executable(source_path)
            archive.writestr(zip_info(name, executable=executable), data)
    if args.out.stat().st_size == 0:
        raise RuntimeError("empty ZIP produced")
    print(json.dumps({
        "status": "PASS",
        "archive": str(args.out),
        "bytes": args.out.stat().st_size,
        "sha256": hashlib.sha256(args.out.read_bytes()).hexdigest(),
        "entries": len(files) + len(generated),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
