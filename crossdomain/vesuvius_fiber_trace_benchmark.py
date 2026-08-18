#!/usr/bin/env python3
"""CPU benchmark for conservative papyrus-fiber tracing on official labeled cubes.

The benchmark is deliberately cube-held-out. Five Scroll-1 cubes train a voxel-cost
model, one Scroll-5 cube selects only the polarity of an intensity baseline, and two
other Scroll-5 cubes are sealed test cubes. NML skeletons provide path ground truth.
No result from this script is new scroll text or ink recovery.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import os
import random
import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.ndimage import (
    binary_dilation,
    gaussian_filter,
    gaussian_gradient_magnitude,
    gaussian_laplace,
    uniform_filter,
)
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
import tifffile

SEED = 20260818
RNG = np.random.default_rng(SEED)


@dataclass
class Cube:
    key: str
    scroll: str
    origin_xyz: tuple[int, int, int]
    size: int
    image_path: Path
    label_path: Path
    nml_path: Path


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def discover_cubes(root: Path) -> list[Cube]:
    images = list(root.rglob("imagesTr/*_0000.tif"))
    labels = {p.name: p for p in root.rglob("labelsTr/*.tif")}
    nmls = list(root.rglob("*.nml"))
    out: list[Cube] = []
    nml_by_key: dict[str, Path] = {}
    nml_re = re.compile(r"fibers_(s\w+?)_(\d+)z_(\d+)y_(\d+)x_(\d+)_v\d+\.nml$")
    for p in nmls:
        m = nml_re.search(p.name)
        if not m:
            continue
        scroll, z, y, x, size = m.groups()
        if scroll == "s1a":
            scroll = "s1"
        key = f"{scroll}_{int(z):05d}_{int(y):05d}_{int(x):05d}_{int(size)}"
        nml_by_key[key] = p
    img_re = re.compile(r"(s\d+)_(\d+)_(\d+)_(\d+)_(\d+)_0000\.tif$")
    for im in images:
        m = img_re.search(im.name)
        if not m:
            continue
        scroll, z, y, x, size = m.groups()
        key = f"{scroll}_{int(z):05d}_{int(y):05d}_{int(x):05d}_{int(size)}"
        lp = labels.get(f"{key}.tif")
        npth = nml_by_key.get(key)
        if lp and npth:
            out.append(Cube(key, scroll, (int(x), int(y), int(z)), int(size), im, lp, npth))
    return sorted(out, key=lambda c: c.key)


def parse_nml(cube: Cube) -> list[dict[str, object]]:
    root = ET.parse(cube.nml_path).getroot()
    fibers: list[dict[str, object]] = []
    ox, oy, oz = cube.origin_xyz
    for thing in root.iter():
        if thing.tag.split("}")[-1] != "thing":
            continue
        tid = thing.attrib.get("id", str(len(fibers)))
        nodes: dict[str, np.ndarray] = {}
        edges: list[tuple[str, str]] = []
        for e in thing.iter():
            tag = e.tag.split("}")[-1]
            if tag == "node" and all(k in e.attrib for k in ("id", "x", "y", "z")):
                # Arrays are z,y,x after the official generator transpose.
                x = float(e.attrib["x"]) - ox
                y = float(e.attrib["y"]) - oy
                z = float(e.attrib["z"]) - oz
                nodes[e.attrib["id"]] = np.array([z, y, x], dtype=float)
            elif tag == "edge" and "source" in e.attrib and "target" in e.attrib:
                edges.append((e.attrib["source"], e.attrib["target"]))
        edges = [(a, b) for a, b in edges if a in nodes and b in nodes]
        if nodes and edges:
            fibers.append({"id": tid, "nodes": nodes, "edges": edges})
    return fibers


def line_voxels(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = max(1, int(math.ceil(float(np.linalg.norm(b - a)) * 1.5)))
    p = np.rint(np.linspace(a, b, n + 1)).astype(int)
    return np.unique(p, axis=0)


def rasterize_fibers(cube: Cube, fibers: Sequence[dict[str, object]]) -> np.ndarray:
    mask = np.zeros((cube.size, cube.size, cube.size), dtype=bool)
    for f in fibers:
        nodes = f["nodes"]
        assert isinstance(nodes, dict)
        for a_id, b_id in f["edges"]:  # type: ignore[index]
            pts = line_voxels(nodes[a_id], nodes[b_id])  # type: ignore[index]
            good = np.all((pts >= 0) & (pts < cube.size), axis=1)
            pts = pts[good]
            if len(pts):
                mask[pts[:, 0], pts[:, 1], pts[:, 2]] = True
    return mask


def robust_normalize(image: np.ndarray) -> np.ndarray:
    x = image.astype(np.float32)
    lo, hi = np.percentile(x, [0.5, 99.5])
    if hi <= lo:
        return np.zeros_like(x)
    return np.clip((x - lo) / (hi - lo), 0.0, 1.0)


def feature_volumes(image: np.ndarray) -> list[np.ndarray]:
    x = robust_normalize(image)
    g1 = gaussian_filter(x, 1.0).astype(np.float32)
    g2 = gaussian_filter(x, 2.0).astype(np.float32)
    grad = gaussian_gradient_magnitude(x, 1.0).astype(np.float32)
    lap = np.abs(gaussian_laplace(x, 1.0)).astype(np.float32)
    mean = uniform_filter(x, size=5).astype(np.float32)
    mean2 = uniform_filter(x * x, size=5).astype(np.float32)
    std = np.sqrt(np.maximum(0, mean2 - mean * mean)).astype(np.float32)
    return [x, g1, g1 - g2, grad, lap, std]


def values_at(features: Sequence[np.ndarray], coords: np.ndarray) -> np.ndarray:
    return np.stack([f[coords[:, 0], coords[:, 1], coords[:, 2]] for f in features], axis=1)


def sample_cube(cube: Cube, max_pos: int = 80_000, neg_ratio: int = 4) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    im = tifffile.imread(cube.image_path)
    fibers = parse_nml(cube)
    center = rasterize_fibers(cube, fibers)
    pos = np.argwhere(center)
    if len(pos) > max_pos:
        pos = pos[RNG.choice(len(pos), max_pos, replace=False)]
    near2 = binary_dilation(center, iterations=2)
    near8 = binary_dilation(center, iterations=8)
    hard = np.argwhere(near8 & ~near2)
    far = np.argwhere(~near8)
    want = max(len(pos), 1) * neg_ratio
    nh = min(len(hard), want // 2)
    nf = min(len(far), want - nh)
    hard = hard[RNG.choice(len(hard), nh, replace=False)] if nh else hard[:0]
    far = far[RNG.choice(len(far), nf, replace=False)] if nf else far[:0]
    neg = np.concatenate([hard, far], axis=0)
    feats = feature_volumes(im)
    X = np.concatenate([values_at(feats, pos), values_at(feats, neg)], axis=0)
    y = np.concatenate([np.ones(len(pos), dtype=np.uint8), np.zeros(len(neg), dtype=np.uint8)])
    return X, y, {"positive": len(pos), "hard_negative": nh, "far_negative": nf, "fiber_count": len(fibers)}


def graph_path_segments(cube: Cube, target_len: float = 48.0, stride: float = 32.0, max_segments: int = 20) -> list[dict[str, object]]:
    fibers = parse_nml(cube)
    segs: list[dict[str, object]] = []
    for f in fibers:
        nodes = f["nodes"]
        edges = f["edges"]
        assert isinstance(nodes, dict) and isinstance(edges, list)
        adj: dict[str, list[str]] = {k: [] for k in nodes}
        for a, b in edges:
            adj[a].append(b); adj[b].append(a)
        # Most annotations are single chains. Walk from each endpoint without revisiting edges.
        endpoints = sorted([k for k, v in adj.items() if len(v) == 1])
        if not endpoints:
            endpoints = [sorted(nodes)[0]]
        used_edges: set[tuple[str, str]] = set()
        for start in endpoints:
            chain = [start]
            prev: str | None = None
            cur = start
            while True:
                nxt = [v for v in adj[cur] if v != prev and tuple(sorted((cur, v))) not in used_edges]
                if not nxt:
                    break
                # Deterministic choice at rare branches: longest immediate edge first.
                nxt.sort(key=lambda v: (-float(np.linalg.norm(nodes[v] - nodes[cur])), v))
                v = nxt[0]
                used_edges.add(tuple(sorted((cur, v))))
                chain.append(v); prev, cur = cur, v
            if len(chain) < 3:
                continue
            pts = [nodes[k] for k in chain]
            cum = [0.0]
            for a, b in zip(pts, pts[1:]):
                cum.append(cum[-1] + float(np.linalg.norm(b - a)))
            total = cum[-1]
            s = 0.0
            while s + target_len <= total + 1e-6:
                e = s + target_len
                i0 = max(0, np.searchsorted(cum, s, side="right") - 1)
                i1 = min(len(pts) - 1, np.searchsorted(cum, e, side="left"))
                path_nodes = pts[i0 : i1 + 1]
                if len(path_nodes) >= 2:
                    vox: list[np.ndarray] = []
                    for a, b in zip(path_nodes, path_nodes[1:]):
                        vox.append(line_voxels(a, b))
                    gt = np.unique(np.concatenate(vox, axis=0), axis=0)
                    good = np.all((gt >= 0) & (gt < cube.size), axis=1)
                    gt = gt[good]
                    if len(gt) >= 15:
                        segs.append({"fiber_id": f["id"], "gt": gt, "start": gt[0], "end": gt[-1]})
                s += stride
    # Favor segments with compact bounding boxes and keep deterministic diversity.
    segs.sort(key=lambda r: (int(np.prod(np.ptp(r["gt"], axis=0) + 1)), str(r["fiber_id"])))  # type: ignore[arg-type]
    if len(segs) > max_segments:
        idx = np.linspace(0, len(segs) - 1, max_segments, dtype=int)
        segs = [segs[i] for i in idx]
    return segs


def route_astar(cost: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    # A compact 26-neighbour A* independent of scikit-image internals.
    shape = cost.shape
    start_t = tuple(int(x) for x in start)
    end_t = tuple(int(x) for x in end)
    moves = [(dz, dy, dx) for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dz, dy, dx) != (0, 0, 0)]
    step_len = {m: math.sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]) for m in moves}
    min_cost = max(1e-6, float(np.nanmin(cost)))
    def heur(p: tuple[int, int, int]) -> float:
        return math.dist(p, end_t) * min_cost
    pq: list[tuple[float, float, tuple[int, int, int]]] = [(heur(start_t), 0.0, start_t)]
    came: dict[tuple[int, int, int], tuple[int, int, int]] = {}
    best = {start_t: 0.0}
    while pq:
        _, g, p = heapq.heappop(pq)
        if p == end_t:
            path = [p]
            while p in came:
                p = came[p]; path.append(p)
            return np.asarray(path[::-1], dtype=int)
        if g != best.get(p):
            continue
        for m in moves:
            q = (p[0] + m[0], p[1] + m[1], p[2] + m[2])
            if not (0 <= q[0] < shape[0] and 0 <= q[1] < shape[1] and 0 <= q[2] < shape[2]):
                continue
            ng = g + step_len[m] * 0.5 * (float(cost[p]) + float(cost[q]))
            if ng < best.get(q, float("inf")):
                best[q] = ng; came[q] = p
                heapq.heappush(pq, (ng + heur(q), ng, q))
    raise RuntimeError("No path")


def straight_path(start: np.ndarray, end: np.ndarray) -> np.ndarray:
    return line_voxels(start.astype(float), end.astype(float))


def path_metrics(pred: np.ndarray, gt: np.ndarray, tolerance: float = 2.0) -> dict[str, float]:
    if len(pred) == 0 or len(gt) == 0:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "mean_pred_to_gt": float("inf"), "mean_gt_to_pred": float("inf")}
    dpg, _ = cKDTree(gt).query(pred, k=1)
    dgp, _ = cKDTree(pred).query(gt, k=1)
    precision = float(np.mean(dpg <= tolerance))
    recall = float(np.mean(dgp <= tolerance))
    f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
    return {"precision": precision, "recall": recall, "f1": f1, "mean_pred_to_gt": float(np.mean(dpg)), "mean_gt_to_pred": float(np.mean(dgp)), "p95_pred_to_gt": float(np.percentile(dpg, 95))}


def trace_segment(image: np.ndarray, model: HistGradientBoostingClassifier, seg: dict[str, object], intensity_polarity: str, margin: int = 8) -> dict[str, object] | None:
    gt = np.asarray(seg["gt"], dtype=int)
    lo = np.maximum(0, gt.min(axis=0) - margin)
    hi = np.minimum(np.asarray(image.shape), gt.max(axis=0) + margin + 1)
    shape = hi - lo
    if np.prod(shape) > 1_200_000 or np.any(shape <= 2):
        return None
    roi = image[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]]
    feats = feature_volumes(roi)
    mat = np.stack([f.ravel() for f in feats], axis=1)
    prob = model.predict_proba(mat)[:, 1].reshape(roi.shape).astype(np.float32)
    start = np.asarray(seg["start"], dtype=int) - lo
    end = np.asarray(seg["end"], dtype=int) - lo
    # Ensure endpoints are traversable and attractive.
    cost_model = -np.log(np.clip(prob, 1e-5, 1.0)) + 0.02
    cost_model[tuple(start)] = 0.001; cost_model[tuple(end)] = 0.001
    x = robust_normalize(roi)
    signal = x if intensity_polarity == "bright" else 1.0 - x
    cost_intensity = -np.log(np.clip(signal, 1e-4, 1.0)) + 0.02
    cost_intensity[tuple(start)] = 0.001; cost_intensity[tuple(end)] = 0.001
    try:
        pm = route_astar(cost_model, start, end) + lo
        pi = route_astar(cost_intensity, start, end) + lo
    except RuntimeError:
        return None
    ps = straight_path(np.asarray(seg["start"], int), np.asarray(seg["end"], int))
    probs_on_path = prob[tuple((pm - lo).T)]
    return {
        "fiber_id": str(seg["fiber_id"]),
        "start": np.asarray(seg["start"], int).tolist(),
        "end": np.asarray(seg["end"], int).tolist(),
        "gt": gt.tolist(),
        "model_path": pm.tolist(),
        "intensity_path": pi.tolist(),
        "straight_path": ps.tolist(),
        "model": path_metrics(pm, gt),
        "intensity": path_metrics(pi, gt),
        "straight": path_metrics(ps, gt),
        "confidence_q25": float(np.quantile(probs_on_path, 0.25)),
        "confidence_mean": float(np.mean(probs_on_path)),
        "roi_shape": shape.tolist(),
    }


def aggregate(records: Sequence[dict[str, object]], method: str) -> dict[str, float]:
    if not records:
        return {"n": 0}
    keys = ["precision", "recall", "f1", "mean_pred_to_gt", "mean_gt_to_pred", "p95_pred_to_gt"]
    out: dict[str, float] = {"n": float(len(records))}
    for k in keys:
        vals = [float(r[method][k]) for r in records if k in r[method]]  # type: ignore[index]
        if vals:
            out[f"mean_{k}"] = float(np.mean(vals))
            out[f"median_{k}"] = float(np.median(vals))
    return out


def confidence_curve(records: Sequence[dict[str, object]]) -> list[dict[str, float]]:
    ranked = sorted(records, key=lambda r: float(r["confidence_q25"]), reverse=True)
    out = []
    for frac in (0.25, 0.5, 0.75, 1.0):
        n = max(1, int(math.ceil(len(ranked) * frac)))
        a = aggregate(ranked[:n], "model")
        a["coverage"] = frac
        a["confidence_cutoff"] = float(ranked[n - 1]["confidence_q25"])
        out.append(a)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset_root")
    ap.add_argument("out_dir")
    ap.add_argument("--segments-per-cube", type=int, default=16)
    args = ap.parse_args()
    root = Path(args.dataset_root).resolve()
    out = Path(args.out_dir).resolve(); out.mkdir(parents=True, exist_ok=True)

    cubes = [c for c in discover_cubes(root) if c.size == 256]
    s1 = [c for c in cubes if c.scroll == "s1"]
    s5 = [c for c in cubes if c.scroll == "s5"]
    if len(s1) < 5 or len(s5) < 3:
        raise RuntimeError(f"Expected at least 5 Scroll-1 and 3 Scroll-5 256-cubes; got {len(s1)}, {len(s5)}")
    train = s1[:5]
    val = [s5[0]]
    test = s5[1:3]

    Xs: list[np.ndarray] = []; ys: list[np.ndarray] = []; sample_info = {}
    for c in train:
        X, y, info = sample_cube(c)
        Xs.append(X); ys.append(y); sample_info[c.key] = info
        print("train sample", c.key, info, flush=True)
    X = np.concatenate(Xs); y = np.concatenate(ys)
    model = HistGradientBoostingClassifier(
        learning_rate=0.08,
        max_iter=180,
        max_leaf_nodes=31,
        l2_regularization=1.0,
        random_state=SEED,
    )
    model.fit(X, y)
    train_prob = model.predict_proba(X)[:, 1]
    model_summary = {
        "train_points": int(len(y)),
        "train_positive_fraction": float(np.mean(y)),
        "train_roc_auc": float(roc_auc_score(y, train_prob)),
        "train_average_precision": float(average_precision_score(y, train_prob)),
        "feature_names": ["raw", "gaussian1", "dog1_2", "gradient1", "abs_laplacian1", "local_std5"],
        "sample_info": sample_info,
    }

    # Select only the bright/dark polarity of an intensity-only comparator on the validation cube.
    val_image = tifffile.imread(val[0].image_path)
    val_segments = graph_path_segments(val[0], max_segments=args.segments_per_cube)
    polarity_scores = {}
    for polarity in ("bright", "dark"):
        rr = [trace_segment(val_image, model, s, polarity) for s in val_segments]
        rr = [r for r in rr if r is not None]
        polarity_scores[polarity] = aggregate(rr, "intensity").get("mean_f1", -1.0)
    polarity = max(polarity_scores, key=polarity_scores.get)

    all_records: list[dict[str, object]] = []
    cube_counts = {}
    for c in test:
        image = tifffile.imread(c.image_path)
        segments = graph_path_segments(c, max_segments=args.segments_per_cube)
        recs = []
        for i, seg in enumerate(segments):
            r = trace_segment(image, model, seg, polarity)
            if r is None:
                continue
            r["cube"] = c.key; r["segment_index"] = i
            recs.append(r); all_records.append(r)
        cube_counts[c.key] = len(recs)
        print("test cube", c.key, len(recs), flush=True)

    summary = {
        "protocol": {
            "seed": SEED,
            "train_cubes": [c.key for c in train],
            "validation_cubes": [c.key for c in val],
            "test_cubes": [c.key for c in test],
            "leakage_control": "No Scroll-5 test cube contributes labels, paths, or voxels to model fitting. The single validation cube selects only bright versus dark polarity for the intensity comparator; the proposed model is not retuned on validation or test paths.",
            "path_tolerance_voxels": 2.0,
            "target_segment_arc_length_voxels": 48.0,
        },
        "dataset_inputs": [{"key": c.key, "image_sha256": file_sha256(c.image_path), "label_sha256": file_sha256(c.label_path), "nml_sha256": file_sha256(c.nml_path)} for c in cubes],
        "model": model_summary,
        "intensity_baseline_polarity": polarity,
        "validation_polarity_scores": polarity_scores,
        "test_segment_counts": cube_counts,
        "test": {
            "proposed_model_astar": aggregate(all_records, "model"),
            "intensity_astar": aggregate(all_records, "intensity"),
            "straight_line": aggregate(all_records, "straight"),
            "confidence_coverage": confidence_curve(all_records) if all_records else [],
        },
        "status": "BENCHMARK_COMPLETED" if all_records else "NO_TEST_SEGMENTS",
        "claim_boundary": "This is a cube-held-out benchmark on an archived labeled positive-control dataset. It is not a new unwrapped surface, winding solution, ink detection, readable text, or prize-winning claim. Any contribution claim requires comparison with current official baselines and current scroll data.",
    }
    (out / "TRACE_RECORDS.json").write_text(json.dumps(all_records), encoding="utf-8")
    (out / "RESULT.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (out / "METRICS.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["cube", "segment_index", "confidence_q25", "confidence_mean", "model_precision", "model_recall", "model_f1", "intensity_precision", "intensity_recall", "intensity_f1", "straight_precision", "straight_recall", "straight_f1"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in all_records:
            w.writerow({
                "cube": r["cube"], "segment_index": r["segment_index"], "confidence_q25": r["confidence_q25"], "confidence_mean": r["confidence_mean"],
                "model_precision": r["model"]["precision"], "model_recall": r["model"]["recall"], "model_f1": r["model"]["f1"],
                "intensity_precision": r["intensity"]["precision"], "intensity_recall": r["intensity"]["recall"], "intensity_f1": r["intensity"]["f1"],
                "straight_precision": r["straight"]["precision"], "straight_recall": r["straight"]["recall"], "straight_f1": r["straight"]["f1"],
            })
    print(json.dumps(summary["test"], indent=2))


if __name__ == "__main__":
    main()
