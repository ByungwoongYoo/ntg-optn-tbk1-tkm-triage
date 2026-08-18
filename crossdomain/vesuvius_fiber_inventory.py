#!/usr/bin/env python3
"""Inventory and positive-control preparation for the public Vesuvius fiber-skeleton dataset.

This script does not claim new fiber recovery. It records the exact archive layout, loads
available NRRD/TIFF/NIfTI arrays, measures label content, and exports deterministic 2-D
previews for a later CPU tracing benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import numpy as np


def sha256(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def summarize_array(a: np.ndarray, *, sample_cap: int = 2_000_000) -> dict[str, Any]:
    if a.size <= sample_cap:
        s = np.asarray(a)
    else:
        step = max(1, int(round((a.size / sample_cap) ** (1 / max(1, a.ndim)))))
        s = np.asarray(a[tuple(slice(None, None, step) for _ in range(a.ndim))])
    finite = np.isfinite(s)
    vals = s[finite]
    out: dict[str, Any] = {
        "shape": list(a.shape),
        "dtype": str(a.dtype),
        "size": int(a.size),
        "sampled_values": int(vals.size),
    }
    if vals.size:
        out.update(
            min=float(vals.min()),
            max=float(vals.max()),
            mean=float(vals.mean()),
            std=float(vals.std()),
            nonzero_fraction=float(np.mean(vals != 0)),
        )
        if np.issubdtype(a.dtype, np.integer) and vals.size <= sample_cap:
            u, c = np.unique(vals, return_counts=True)
            if len(u) <= 256:
                out["unique_values"] = [int(x) for x in u.tolist()]
                out["unique_counts"] = [int(x) for x in c.tolist()]
                out["positive_instance_count"] = int(np.sum(u > 0))
    return out


def load_array(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    low = path.name.lower()
    if low.endswith(".nrrd") or low.endswith(".nhdr"):
        import nrrd

        a, header = nrrd.read(str(path), index_order="C")
        h = {str(k): str(v) for k, v in header.items()}
        return np.asarray(a), {"loader": "pynrrd", "header": h}
    if low.endswith(".nii") or low.endswith(".nii.gz"):
        import nibabel as nib

        im = nib.load(str(path))
        return np.asarray(im.dataobj), {
            "loader": "nibabel",
            "affine": np.asarray(im.affine).tolist(),
        }
    if low.endswith((".tif", ".tiff")):
        import tifffile

        return np.asarray(tifffile.imread(str(path))), {"loader": "tifffile"}
    if low.endswith(".npy"):
        return np.load(path, mmap_mode="r"), {"loader": "numpy"}
    raise ValueError(f"Unsupported array format: {path}")


def save_preview(a: np.ndarray, path: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    b = np.asarray(a)
    if b.ndim == 3:
        # Use an occupancy maximum projection. For continuous CT data, clip robustly.
        img = np.max(b, axis=0)
    elif b.ndim == 2:
        img = b
    else:
        return
    fig, ax = plt.subplots(figsize=(7, 7))
    if np.issubdtype(img.dtype, np.integer) and np.unique(img[:: max(1, img.shape[0] // 512), :: max(1, img.shape[1] // 512)]).size < 128:
        ax.imshow(img, cmap="nipy_spectral", interpolation="nearest")
    else:
        lo, hi = np.percentile(img[np.isfinite(img)], [1, 99])
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def classify(path: Path) -> str:
    s = str(path).lower()
    if any(x in s for x in ("labelstr", "labels", "mask", "seg")):
        return "label"
    if any(x in s for x in ("imagetr", "images", "volume", "ct")):
        return "image"
    if path.suffix.lower() == ".nml":
        return "skeleton_nml"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset_root")
    ap.add_argument("out_dir")
    ap.add_argument("--max-arrays", type=int, default=24)
    args = ap.parse_args()

    root = Path(args.dataset_root).resolve()
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "previews").mkdir(exist_ok=True)

    files = [p for p in root.rglob("*") if p.is_file()]
    tree_lines: list[str] = []
    by_ext: dict[str, int] = {}
    by_class: dict[str, int] = {}
    for p in sorted(files):
        rel = p.relative_to(root)
        cls = classify(rel)
        by_class[cls] = by_class.get(cls, 0) + 1
        suffix = ".nii.gz" if p.name.lower().endswith(".nii.gz") else p.suffix.lower()
        by_ext[suffix] = by_ext.get(suffix, 0) + 1
        tree_lines.append(f"{p.stat().st_size}\t{cls}\t{rel}")
    (out / "DATASET_TREE.tsv").write_text("bytes\tclass\tpath\n" + "\n".join(tree_lines) + "\n", encoding="utf-8")

    array_exts = (".nrrd", ".nhdr", ".nii", ".nii.gz", ".tif", ".tiff", ".npy")
    candidates = [p for p in files if p.name.lower().endswith(array_exts)]
    # Deterministic: inspect labels first, then images, then other arrays, smaller files first.
    candidates.sort(key=lambda p: ({"label": 0, "image": 1}.get(classify(p.relative_to(root)), 2), p.stat().st_size, str(p)))

    arrays: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    for idx, p in enumerate(candidates[: args.max_arrays]):
        rel = p.relative_to(root)
        rec: dict[str, Any] = {
            "path": str(rel),
            "class": classify(rel),
            "bytes": p.stat().st_size,
            "sha256": sha256(p),
        }
        try:
            a, meta = load_array(p)
            rec.update(meta)
            rec["array"] = summarize_array(a)
            preview = out / "previews" / f"{idx:03d}_{p.name.replace('.', '_')}.png"
            save_preview(a, preview, str(rel))
            if preview.exists():
                rec["preview"] = str(preview.relative_to(out))
        except Exception as exc:  # preserve failure evidence rather than silently skipping
            rec["error"] = repr(exc)
            errors.append({"path": str(rel), "error": repr(exc)})
        arrays.append(rec)

    nmls = [p for p in files if p.suffix.lower() == ".nml"]
    nml_summary: list[dict[str, Any]] = []
    for p in sorted(nmls)[:20]:
        text = p.read_text(encoding="utf-8", errors="ignore")
        nml_summary.append(
            {
                "path": str(p.relative_to(root)),
                "bytes": p.stat().st_size,
                "sha256": sha256(p),
                "node_tags": int(text.count("<node ")),
                "edge_tags": int(text.count("<edge ")),
                "thing_tags": int(text.count("<thing ")),
            }
        )

    result = {
        "dataset_root": str(root),
        "file_count": len(files),
        "total_bytes": int(sum(p.stat().st_size for p in files)),
        "by_extension": dict(sorted(by_ext.items())),
        "by_class": dict(sorted(by_class.items())),
        "array_candidates": len(candidates),
        "arrays_inspected": len(arrays),
        "array_records": arrays,
        "nml_count": len(nmls),
        "nml_records_first20": nml_summary,
        "errors": errors,
        "status": "POSITIVE_CONTROL_DATASET_INVENTORIED" if arrays and not all("error" in r for r in arrays) else "INVENTORY_FAILED",
        "claim_boundary": "This run inventories an archived, manually labeled positive-control dataset. It does not recover new scroll surface, ink, or text and does not by itself constitute a Vesuvius Prize submission.",
    }
    (out / "INVENTORY.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("file_count", "total_bytes", "by_extension", "by_class", "array_candidates", "arrays_inspected", "nml_count", "status")}, indent=2))


if __name__ == "__main__":
    main()
