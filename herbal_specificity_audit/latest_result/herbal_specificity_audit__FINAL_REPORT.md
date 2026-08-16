# Superseded development-pilot report

> **Do not cite this directory as the final analysis.**

The material originally stored in `latest_result/` came from a **500-record development pilot**. That pilot was used to refine deterministic extraction rules and debug the confidence-interval procedure. It is superseded by the held-out full frozen-corpus analysis completed on 2026-08-16.

## Authoritative final result

See:

- `../final_v3/README.md`
- `../final_v3/FINAL_REPORT.md`
- `../final_v3/results_key.json`

Successful GitHub Actions run: `31949335813`  
Run commit: `93f099ecc0112de254f5c312783baf272b07bc7f`  
Artifact ID: `9264344939`  
Artifact SHA-256: `8f73a70289338b6ed53392808124c95a3353c0a3359135a6e33a5c7f5678a769`

## Final prespecified decision

`no_reliable_positive_specificity_signal`

The held-out corpus did not show a reliable positive within-disease similarity advantage. Under the frozen corpus, strict extraction rule, disease labels, and prespecified controls, author-reported hub/core/key target-gene lists did **not** provide a reproducible disease-specific fingerprint.

Primary held-out result:

- 114 articles across 9 disease classes
- within-disease mean Jaccard: 0.1219115689
- between-disease mean Jaccard: 0.1214275645
- delta: 0.0004840045
- one-sided label-permutation p: 0.4335
- stratified article-bootstrap 95% CI: [-0.0124, 0.0140]

The earlier pilot output is retained only as development provenance. Its numerical conclusion must not be presented as the study's final result.
