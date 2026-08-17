#!/usr/bin/env python3
"""Final compatibility runner for the ClinVar time-machine analysis.

This keeps the frozen protocol and all estimands unchanged. It adds only a structural
empty-cohort guard so a prespecified cutoff with zero eligible future-resolved VUS is
recorded as n=0 instead of causing scikit-learn to abort on predict_proba([]).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import clinvar_time_machine_v2_runner as compat

# The compatibility runner has already installed its source-format and vectorized
# cluster-bootstrap patches on the shared core module.
tm = compat.tm
_original_apply_temporal_recipe = tm.apply_temporal_recipe


def apply_temporal_recipe_empty_safe(model, test: pd.DataFrame) -> pd.DataFrame:
    if len(test) == 0:
        out = test.copy()
        # Preserve the exact output contract expected by downstream cohort tables.
        for col in [
            "temporal_ensemble_raw",
            "temporal_poet_raw",
            "temporal_ensemble_prob",
            "temporal_poet_prob",
        ]:
            out[col] = pd.Series(index=out.index, dtype=float)
        return out
    return _original_apply_temporal_recipe(model, test)


tm.apply_temporal_recipe = apply_temporal_recipe_empty_safe

if __name__ == "__main__":
    tm.main()
