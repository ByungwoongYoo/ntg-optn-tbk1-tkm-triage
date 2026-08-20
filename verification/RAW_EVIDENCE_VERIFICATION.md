# B2[2]/Z100 raw-evidence verification

**Status: PASS**

- GitHub source artifacts: 373
- final v10 t=2,3,4 subbranches: 255
- final v8 t=5,...,45 branches: 41
- v9 g=2,3 branches: 67
- actual g=4,...,7 CNF/DRAT cases with original and fresh VERIFIED logs: 4
- valid 13-set: `[0, 5, 7, 31, 58, 61, 62, 63, 72, 80, 84, 91, 97]`

## Warnings
- historical_self_hash_lines_skipped: `370 artifact-local files hashed themselves while being written; the corrected release manifest supersedes those self-lines.`
- g1_v8_t2_t4_historical_superseded: `{2: {'status': 'EXHAUSTED_NO_WITNESS', 'completed': True, 'timed_out': False}, 3: {'status': 'INCOMPLETE_TIMEOUT', 'completed': False, 'timed_out': True}, 4: {'status': 'EXHAUSTED_NO_WITNESS', 'completed': True, 'timed_out': False}}`
- historical_preupload_manifest_not_authoritative: `{'listed': 4132, 'missing_after_actions_upload': 25, 'missing_examples': ['repository_snapshot/.github/workflows/b2-cube-status-collector-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap-cases-v2.yml', 'repository_snapshot/.github/workflows/b2-gap-dfs-v5.yml', 'repository_snapshot/.github/workflows/b2-gap-dfs-v6.yml', 'repository_snapshot/.github/workflows/b2-gap-sat-v4.yml', 'repository_snapshot/.github/workflows/b2-gap-vector-v3.yml', 'repository_snapshot/.github/workflows/b2-gap1-branches-v7.yml', 'repository_snapshot/.github/workflows/b2-gap1-existing-aggregate-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap1-u-split-v8-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap1-u-split-v8-control-tail-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap1-v-split-v10-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap1-v10-result-collector-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap1-v7-status-collector-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap23-u-split-v9-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap23-v9-result-collector-20260817.yml', 'repository_snapshot/.github/workflows/b2-gap4-drat-proof-v1.yml', 'repository_snapshot/.github/workflows/b2-gap567-drat-20260817.yml', 'repository_snapshot/.github/workflows/b2-hard-gap-reflected-cubes-20260817.yml', 'repository_snapshot/.github/workflows/b2-missing-gap-drat-20260817.yml', 'repository_snapshot/.github/workflows/b2-neighborhood-exact-rerun-20260817.yml', 'repository_snapshot/.github/workflows/b2-z100-fast.yml', 'repository_snapshot/.github/workflows/b2-z100-final-proof-bundle-20260817.yml', 'repository_snapshot/.github/workflows/b2-z100-final-record-check-20260817.yml', 'repository_snapshot/.github/workflows/b2-z100-probe.yml', 'repository_snapshot/.github/workflows/b2-z100-raw-evidence-package-20260818.yml'], 'explanation': 'actions/upload-artifact omitted hidden .github paths. Exact WORKFLOW_AT_RUN snapshots and the corrected release manifest supersede it.'}`

## Failures
- none
