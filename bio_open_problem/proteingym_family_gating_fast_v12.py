#!/usr/bin/env python3
"""Fast prespecified family-held-out experiment: primary global30 clustering + linear gate.

This wrapper intentionally narrows the candidate set before seeing any result:
- fixed ensembles: mean3, mean5, weighted3_p1, weighted5_p1;
- dynamic gate: ridge reliability gate only;
- direct stack: logistic C=0.1 only.
It reuses the complete-case, cluster-held-out implementation in v12.
"""
from __future__ import annotations

import proteingym_family_gating_v12 as core

core.FIXED_SPECS = [
    core.FixedSpec("mean3", "mean", 3),
    core.FixedSpec("mean5", "mean", 5),
    core.FixedSpec("weighted3_p1", "weighted", 3, 1.0),
    core.FixedSpec("weighted5_p1", "weighted", 5, 1.0),
]
core.GATE_SPECS = [core.GateSpec("ridge_gate_t010", "ridge", 0.10)]
core.DIRECT_SPECS = [core.DirectSpec("direct_logit_c01", 0.1)]

if __name__ == "__main__":
    core.main()
