#!/usr/bin/env python3
"""Compatibility entry point for the canonical B2/Z100 prefix generator.

Adaptive workflows originally referenced the v3 path. They now use the v5
implementation, which retains the same exact prefix semantics and original
B2[2] constraints while adding only the parity and quotient-capacity cuts
proved in `B2_Z100_REDUNDANT_CUTS_PROOF.md`.
"""
from b2_z100_prefix_proof_cnf_v5 import main

if __name__ == '__main__':
    main()
