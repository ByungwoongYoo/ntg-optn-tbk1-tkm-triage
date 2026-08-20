# LANTERN v1–v7 corrected audit and CAMI submission freeze

## Canonical mapping

The only accepted CAMI III Toy mapping is:

`0/1, 2/7, 3/5, 4/8, 6/19, 9/12, 10/15, 11/13, 14/16, 17/18`.

Sequential adjacency and similarity-inferred pairing are prohibited. PR #30 used all nine corrected pairs after development on 0/1, so no project-wide untouched public-Toy pair remains. Further v8/v9, threshold-sweep and rescue-intensity tuning on public Toy is stopped.

## Corrected v1–v7 evidence table

| Version/stage | Samples | Classification | Corrected scientific status | Direct evidence |
|---|---|---|---|---|
| v1 core development | 0/1 | development only | Valid same-individual development. GF −0.034942 pp; mean +0.036974 pp; strict success not met. | [PR #21](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/21), [run 32305282908](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32305282908), [artifact 9376226761](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32305282908/artifacts/9376226761) |
| v1 low-depth reapplication | 2/3, 4/5, 6/7, 8/9 | cross-individual stress | Not longitudinal validation. | [Corrected PR #22](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/22) |
| v2 additive development | 0/1 | development only | Valid same-individual development. GF +0.178909 pp; mean +0.197045 pp; chimera relative +1.615%; strict success not met. | [PR #23](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/23), [run 32306563734](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32306563734), [artifact 9377068681](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32306563734/artifacts/9377068681) |
| v2 claimed holdouts | 2/3, 4/5, 6/7, 8/9 | cross-individual stress | Frozen-rule stress only; not independent longitudinal validation. | [Corrected PR #25](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/25) |
| v3 development | 0/1 | development only | Valid same-individual development; C0394 frozen before corrected-pair evaluation. | [PR #27](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/27) |
| v3 corrected multi-pair audit | 2/7, 3/5, 4/8, 6/19, 9/12, 10/15, 11/13, 14/16, 17/18 | valid longitudinal | Only valid project-wide multi-pair longitudinal public-Toy audit. Aggregate strict verdict FAIL; mean GF +0.028590 pp; mean recovery +0.028209 pp; max relative chimera +13.491%. | [Corrected PR #30](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/30), [run 32314301866](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32314301866), [aggregate artifact 9389187289](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32314301866/artifacts/9389187289) |
| v4 hybrid development | 0/1 | development only | Valid same-individual development. GF +0.219212 pp; mean +0.291927 pp; low abundance +0.087327 pp; strict success not met. | [Corrected PR #35](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/35), [run 32332720793](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32332720793), [artifact 9393830269](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32332720793/artifacts/9393830269) |
| v4/v6 14/15 executions | 14/15 | cross-individual stress | 14 belongs with 16; 15 belongs with 10. Final v6 stress verdict NEGATIVE; no longitudinal-ablation claim. | [PR #38](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/38), [PR #43](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/43), [finalizer PR #47](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/47) |
| v5 long-read attempt | 2/3; planned 16/17 and 18/19 | workflow failure | Cross-individual design and workflow failure; no defensible performance result. | [Corrected PR #37](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/37), [run 32333700973](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32333700973) |
| v6 high-depth development | 0/1 | development only | Valid development probe only; no independent Toy holdout remains. | [PR #42](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/42) |
| v7 short-backbone development | 0/1 | development only | Winner `short_terminal_500_s1`; GF +0.175196 pp; mean +0.443214 pp; low +0.285919 pp; strict success not met. | [Corrected PR #48](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/48), [run 32351114616](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32351114616), [artifact 9401163313](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32351114616/artifacts/9401163313) |
| v7 finalizer | 18/19 | cross-individual stress | `PARTIAL`, `strict_success=false`; 18 belongs with 17 and 19 with 6; both were already exposed in v3. | [Re-audited PR #59](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/pull/59), [run 32369519543](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32369519543), [artifact 9408725398](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32369519543/artifacts/9408725398) |

## Independent v7 artifact arithmetic

Artifact 9408725398 has GitHub digest `ffeea8901c0c3b56bdcba82f0a8da0a5a32cb3ac5dcb2d1e245b9f7512c3d703`. The archive and its internal SHA manifest passed integrity checks.

- strongest baseline: `metaspades_hybrid_pair`
- GF gain: `+0.19705271268690616 pp`
- mean recovery gain: `+0.8235269478120657 pp`
- low-abundance gain: `+0.343797130120449 pp`
- relative chimera change: `−0.44187475457456615` (−44.1875%)
- former `longitudinal_ablation` value: `+0.8173546210815257 pp`, corrected label `cross_individual_coassembly_contribution`
- final status: `PARTIAL`; `strict_success=false`

Deficits to the frozen gates:

- GF: `0.500000 − 0.19705271268690616 = 0.30294728731309384 pp`
- mean: `1.000000 − 0.8235269478120657 = 0.1764730521879343 pp`
- low abundance: `2.000000 − 0.343797130120449 = 1.656202869879551 pp`

A green GitHub Actions result indicates technical execution success only.

## Frozen submission candidate

The active-challenge candidate is the exact v7 rule `short_terminal_500_s1`, frozen at source commit `d2f44b529d2198784ce7666ba1c98ae44709f981`. Parameters are in `configs/frozen_v7_rule.json`; grouping is loaded only from explicit metadata by `lantern_v7/infer_longitudinal_pairs.py` (now a metadata loader, not an inference algorithm).

Build and dry-run commands are in `README.md` and `run_frozen_v7.sh`. The package produces per-individual metaSPAdes-short, metaSPAdes-hybrid, MEGAHIT and Flye assemblies; applies the unchanged terminal-extension rule; combines unique prefixed IDs; validates uppercase `ATCGN`; and records hashes and versions.

## Go/no-go

Upload only after the active mapping, input manifest, frozen commit, container image, command plan and final FASTA hashes are approved. If official CAMI evaluation is negative or lacks a robust advantage, terminate the positive assembler claim and convert the study to a negative benchmark/methodology note on recovery–chimera trade-offs, public-benchmark overfitting and sample-mapping/validation errors.
