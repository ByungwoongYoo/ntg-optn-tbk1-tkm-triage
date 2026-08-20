# LANTERN CAMI III official-mapping freeze package

This package supersedes the project’s incorrect sequential Toy pairing assumptions and freezes one unchanged LANTERN-v7 candidate for **active CAMI III blind execution only**.

## Scientific status

- Public Toy micro-tuning is permanently stopped.
- The corrected Toy mapping is stored in `config/official_toy_mapping.*`.
- All corrected Toy individuals were already evaluated by v3; no clean public-Toy holdout remains project-wide.
- V7 samples 18/19 are a cross-individual/coassembly stress test, not a longitudinal or untouched holdout.
- Historical artifacts and failures are retained; current interpretation is governed by `evidence/claim_corrections.*`.
- The candidate is submission-ready in the engineering sense: frozen code, parameters, container recipe, manifest-driven input grouping, FASTA validation, hashes, and recovery procedure are supplied.
- The active challenge has **not** been run or submitted from this package because restricted challenge inputs/portal authorization were not provided.

## Frozen execution

1. Fill `config/challenge_sample_manifest.template.tsv` with the actual challenge files and the official individual IDs supplied with the challenge data.
2. Build the container.
3. Run `run/run_frozen_v7.sh`.
4. Validate all output FASTAs.
5. Review `submission/SUBMISSION_CHECKLIST.md`.
6. Upload only after explicit approval.

No sample-number adjacency or sequential pairing inference exists in the frozen runner.
