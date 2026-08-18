#!/usr/bin/env python3
"""Corrected entry point for the Vesuvius fiber benchmark.

v1 used numpy.unique on interpolated line voxels. numpy.unique sorts rows, which
silently destroyed path order and could change the intended endpoints. This wrapper
replaces that helper with order-preserving consecutive de-duplication before invoking
the otherwise frozen v1 benchmark. The correction affects ground-truth path ordering
and the straight-line comparator; all split and feature-model rules remain unchanged.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vesuvius_fiber_trace_benchmark as core


def ordered_line_voxels(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = max(1, int(np.ceil(float(np.linalg.norm(b - a)) * 1.5)))
    p = np.rint(np.linspace(a, b, n + 1)).astype(int)
    if len(p) <= 1:
        return p
    keep = np.ones(len(p), dtype=bool)
    keep[1:] = np.any(p[1:] != p[:-1], axis=1)
    return p[keep]


core.line_voxels = ordered_line_voxels

if __name__ == "__main__":
    core.main()
