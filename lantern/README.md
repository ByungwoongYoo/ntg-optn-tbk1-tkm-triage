# LANTERN — CAMI III execution workspace

LANTERN stands for **Longitudinal Assembly of Novel Taxa with Evidence-based Rescue and Non-leaky Validation**.

This branch is an isolated execution workspace for the CAMI III Longitudinal Human Gut assembly challenge. It is not merged into `main` and must not be described as a successful CAMI method until the pre-specified gates have been completed.

## Verified CAMI III facts (snapshot: 2026-08-20)

- CAMI III Longitudinal Human Gut opened on 2026-07-14.
- Assembly submissions close on 2026-10-14.
- The challenge contains 9 individuals × 4 time points, with short- and long-read data.
- The Toy Longitudinal Human Gut benchmark contains 10 individuals × 2 time points (20 samples), 100 Gbp short reads and 100 Gbp long reads, plus gold standards.
- The active challenge accepts global cross-sample, patient-specific cross-sample, or single-sample assemblies, using short reads, long reads, or both.
- A first-ten-sample submission route is explicitly permitted for methods that cannot process the full dataset.
- Assembly output is FASTA; sequence lines may contain only uppercase `A`, `T`, `C`, `G`, or `N`.
- Active-challenge online gold-standard evaluation is disabled.
- Reproducibility metadata must include software versions, parameters, databases, seeds, and run-combination/selection rules.

Authoritative URLs and the frozen statements are recorded in `config/cami3_official_snapshot_20260820.json`.

## Non-leaky execution order

1. Official-source and file-format smoke test.
2. Toy-data acquisition and exact manifest.
3. Strong single-sample, individual co-assembly, and pooled baselines.
4. Pseudo-novel benchmark construction.
5. LANTERN longitudinal rescue implementation.
6. Ablation and adversarial stress tests.
7. Freeze code/configuration/thresholds.
8. Blind application to CAMI III challenge reads.
9. Submission-format validation and final evidence package.

## Current claim boundary

No baseline result, LANTERN improvement, blind submission, or novel-genome recovery has yet been established on this branch. The first workflow only verifies official data retrieval, archive contents, FASTA integrity, and provenance before computationally expensive assembly work begins.
