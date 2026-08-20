# LANTERN CAMI III frozen submission package

Status: engineering package only. No active CAMI reads have been downloaded or run, no portal submission has occurred, and no positive assembler claim is established.

## Production mapping rule

`lantern_submission/scripts/load_explicit_mapping.py` is the only production grouping path. It requires explicit TSV columns:

- `individual_id`
- `sample_id`
- `timepoint`

It rejects duplicate samples, duplicate individual/timepoint rows, wrong group sizes, missing/extra samples, and mapping/manifest disagreement. It writes deterministic `MAPPING_FREEZE.json` and SHA-256 records. Sample order, sequential numbering, Mash and sequence similarity are never used.

Historical Mash-inference code remains only in Git history/provenance and is not called by this package or its CI.

## Frozen candidate

- Candidate: `LANTERN-v7-short-terminal-500-s1`
- Algorithm source commit: `d2f44b529d2198784ce7666ba1c98ae44709f981`
- Parameters: `configs/frozen_v7_rule.json`
- Result-dependent retuning: prohibited
- Public Toy tuning: stopped

## Private preflight

Populate `configs/input_manifest_TEMPLATE.tsv` with private absolute paths and the actual SHA-256 of every read file. Do not commit the populated file.

```bash
bash lantern_submission/run_frozen_v7.sh \
  --mapping /private/active_mapping.tsv \
  --manifest /private/input_manifest.tsv \
  --out /private/lantern_run \
  --threads 32 \
  --memory-gb 256 \
  --expected-individuals 9 \
  --expected-timepoints 4 \
  --dry-run
```

The dry run verifies explicit mapping membership, exact sample coverage, file existence, gzip integrity and SHA-256; freezes the mapping, inputs and rule; and creates `RUN_PLAN.json`, `RUN_PLAN.sh`, and `MANIFEST.sha256` without invoking assemblers.

## Frozen execution

After separate authorization for restricted-data staging and execution, rerun the same command without `--dry-run`. Outputs remain on private local/HPC storage. The final candidate FASTA is:

```text
/private/lantern_run/submission/LANTERN_CAMI3_ASSEMBLY.fasta
```

It must pass `lantern_cami3/scripts/verify_cami_assembly.py` before any separate portal-submission approval is requested.

## Security boundary

Never upload active reads, populated active manifests, derived active assemblies, checksums that reveal unapproved metadata, or logs containing restricted paths to GitHub Actions artifacts, the public repository, Zenodo or another public service.

A green GitHub Actions run establishes technical package execution only. It is not scientific superiority, CAMI rank, or official evaluation.
