# Resource and recovery gate

## Status

`BLOCKED_RESOURCES` until the actual private execution host reports CPU, RAM and scratch capacity and the active input manifest reports real byte counts.

The repository does not assume that a public GitHub runner can execute the restricted 36-sample analysis. Active reads and derived outputs are prohibited from GitHub Actions and public artifacts.

## Required preflight calculation

After private staging, calculate:

- total compressed and uncompressed short-read bytes
- total compressed and uncompressed long-read bytes
- available scratch bytes
- available RAM per concurrent individual job
- available CPU threads and scheduler wall-time limits

The full prespecified run covers all 9 explicit individuals and 36 samples. If it cannot run with the available private resources, stop with `BLOCKED_RESOURCES`; do not silently switch to first-ten samples or another subset.

Conservative planning targets, to be replaced by measured private preflight values:

- 32–64 CPU threads per individual job
- 256–512 GB RAM per metaSPAdes hybrid job
- 1–2 TB scratch per individual job
- 8–12 TB aggregate scratch if multiple individuals are retained concurrently

These are engineering estimates, not measured active-CAMI requirements.

## Recovery policy

1. Keep each individual in an isolated work directory.
2. Preserve completed assembler outputs and SHA-256 records.
3. On technical failure, rerun only the failed command with the identical arguments and inputs.
4. Do not alter thresholds, seed, candidate sources or postprocessing.
5. Record original and resumed stdout/stderr, exit codes, timestamps and scheduler resource reports.
6. If an output checksum changes after a nominal resume, invalidate that individual and investigate before proceeding.
7. Combine FASTAs only after all 9 individual validators pass.
8. Never upload restricted logs, manifests, reads or assemblies to a public artifact store.
