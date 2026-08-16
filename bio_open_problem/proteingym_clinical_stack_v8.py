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


core.weights_per_protein_class = fixed_weights_per_protein_class

if __name__ == '__main__':
    core.main()
