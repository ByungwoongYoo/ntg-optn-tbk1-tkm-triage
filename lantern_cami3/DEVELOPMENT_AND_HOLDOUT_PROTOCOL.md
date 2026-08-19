# LANTERN development and held-out validation protocol

## Fixed units
- Development individual: Toy samples 0 and 1.
- Untouched held-out individual: Toy samples 2 and 3.
- Optional second held-out individual: Toy samples 4 and 5.
- The challenge dataset is not a development set and remains inaccessible until the human participant accepts CAMI's data agreement.

## Development rules
The v1 method and numerical decision gates were frozen before public Toy gold-standard performance was read. Results from samples 0/1 may be used to diagnose execution defects and make one documented v2 method revision. Every v2 change must cite a failure mode observed on the development individual and must be committed before any sample2/3 gold-standard or read-to-genome mapping is accessed.

## Held-out rules
All sample2/3 assembly, clustering, read-evidence scoring, candidate selection, ablation, scaffold decisions, software versions, parameters, and hashes are frozen before sample2/3 GSA or read-to-genome mapping is downloaded. No v2 threshold may be changed after the held-out result is opened. A negative held-out result remains negative.

## Primary comparison
LANTERN_full is compared with every successfully executed strong baseline. The primary baseline reference is the strongest baseline on genome fraction, with mean per-genome recovery and conservative cross-BINID chimera proxy as tie breakers. Longitudinal contribution is assessed by the prespecified LANTERN_no_longitudinal ablation.

## Secondary comparisons
- short-only, long-only, and hybrid baselines;
- no-long and no-consensus ablations;
- conservative long-read scaffolding as a secondary method, never substituted post hoc for LANTERN_full;
- evaluation-only pseudo-novel target-withholding tiers;
- abundance, domain, GC, fragmentation, and repeat-proxy strata.

## Claim boundary
A held-out Toy improvement demonstrates reproducible performance on a public synthetic benchmark. It is not an official CAMI III rank, not a recovery of an unpublished challenge genome, and not a discovery of a new organism. Those claims require restricted challenge submission and official evaluation.
