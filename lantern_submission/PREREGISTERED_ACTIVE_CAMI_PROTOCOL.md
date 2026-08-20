# Preregistered active CAMI protocol

## Current status

`PAUSE_BLOCKED`

The package may be tested only with synthetic or public Toy inputs until the human participant personally accepts the CAMI data agreement and private compute/storage are confirmed. No restricted data may enter a public runner or artifact.

## Fixed active grouping

| Individual | Samples in timepoint order |
|---|---|
| 1 | 0, 2, 5, 19 |
| 2 | 1, 4, 11, 30 |
| 3 | 3, 8, 20, 24 |
| 4 | 6, 7, 9, 31 |
| 5 | 10, 28, 29, 32 |
| 6 | 12, 14, 18, 34 |
| 7 | 13, 22, 25, 26 |
| 8 | 15, 16, 27, 33 |
| 9 | 17, 21, 23, 35 |

Grouping is read only from explicit metadata. Sequential numbering and similarity inference are forbidden.

## Frozen candidate

- `LANTERN-v7-short-terminal-500-s1`
- algorithm source commit: `d2f44b529d2198784ce7666ba1c98ae44709f981`
- seed: `20260820`
- rule file: `configs/frozen_v7_rule.json`
- no threshold, source, seed, postprocessing or selection change after active input access or output inspection

## Pre-execution gates

All must pass before a restricted run:

1. Explicit mapping loader PASS with 9 individuals, 4 timepoints, 36 exact samples.
2. Populated private input manifest contains all 36 samples once.
3. Every input exists, is nonempty, has valid gzip structure and matches its frozen SHA-256.
4. Available compute is sufficient for the complete prespecified patient-specific run; otherwise record `BLOCKED_RESOURCES` rather than silently switching to a subset.
5. Container digest, source commit, rule SHA-256, mapping SHA-256, manifest SHA-256 and generated plan SHA-256 are frozen.
6. The participant separately authorizes restricted-data staging/execution.

## Execution

Run `run_frozen_v7.sh` unchanged on private local/HPC storage. Preserve:

- command lines and generated `RUN_PLAN.sh`
- stdout/stderr
- start/end time and wall time
- peak memory and hardware/software details where the scheduler provides them
- input, intermediate final-assembly and output hashes
- per-individual FASTA validation
- combined FASTA validation

No result-dependent rerun, seed selection or parameter change is allowed.

## Submission gate

The validated output must remain private until a second, distinct user authorization permits portal submission. Immediately after submission, status is only `SUBMITTED_AWAITING_OFFICIAL_EVALUATION`. A green local or GitHub run is not scientific success.

## Terminal decisions

- `ADVANCE_TO_MANUSCRIPT`: official valid evaluation meets preregistered success and guardrail criteria.
- `STOP_NEGATIVE`: official valid evaluation does not meet criteria; positive assembler claim ends and Toy tuning does not resume.
- `STOP_INVALID`: leakage, post-result change, wrong mapping or invalid submission is found.
- `PAUSE_BLOCKED`: account, agreement, data, resource or submission-policy requirements remain unresolved.
