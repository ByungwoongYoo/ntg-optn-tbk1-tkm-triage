#!/usr/bin/env python3
"""Corrected runner for the v7 protein-held-out clinical stacking benchmark."""
from __future__ import annotations

import numpy as np
import proteingym_clinical_stack_v7 as core


def fixed_weights_per_protein_class(df):
    work = df.reset_index(drop=True)
    w = np.zeros(len(work), dtype=float)
    for idx in work.groupby('protein_file', sort=False).indices.values():
        idx = np.asarray(idx, dtype=int)
        y = work.iloc[idx]['label'].to_numpy(dtype=int)
        for cls in (0, 1):
            ii = idx[y == cls]
            if len(ii):
                w[ii] = 0.5 / len(ii)
    if not np.isfinite(w).all() or w.sum() <= 0:
        raise RuntimeError('Invalid protein/class balancing weights')
    w *= len(w) / w.sum()
    return w


def fixed_oriented(df, models, signs):
    # pandas 3 may return a read-only NumPy view; request a writable copy.
    X = df[models].to_numpy(dtype=float, copy=True)
    for j, model in enumerate(models):
        if signs[model] < 0:
            X[:, j] = 1.0 - X[:, j]
    return X


core.weights_per_protein_class = fixed_weights_per_protein_class
core.oriented = fixed_oriented

if __name__ == '__main__':
    core.main()
