# Official CAMI III rules audit

Checked against the CAMI primary pages on 2026-08-21.

## Confirmed

- Active Longitudinal Human Gut data: 9 individuals × 4 timepoints, 36 short-read and 36 long-read samples.
- Explicit active mapping is stored in `configs/active_official_mapping.tsv`.
- Assembly accepts global cross-sample, patient-specific cross-sample, or single-sample modes; short, long, or combined data may be used.
- A first-ten-sample cross-assembly or ten single-sample assemblies is permitted for methods unable to process the full dataset.
- Assembly deadline: 2026-10-14.
- Active-challenge online evaluation is disabled.
- Assembly output is FASTA; headers use the published character set and sequences must match `[ATCGN]+` in uppercase.
- Reproducibility metadata must describe exact versions, parameters, databases, seed behavior and any run combination/selection.
- Active data are restricted, may not be redistributed, published, deposited, or used outside the challenge before restrictions are lifted.
- Dataset download requires human login and acceptance of the CAMI Data Download Agreement; the public page currently instructs users to log in with a university email address.

Primary pages:

- https://cami-challenge.org/datasets/human-gut/
- https://cami-challenge.org/cami-iii-challenges/
- https://cami-challenge.org/submit/
- https://cami-challenge.org/file-formats/
- https://cami-challenge.org/faq/
- https://cami-challenge.org/schedule/

## Not confirmed by the public pages

- Whether patient-specific assemblies are uploaded as nine FASTAs, one archive, one combined FASTA, or nine submissions.
- Whether candidate and frozen baseline can both be submitted and the exact submission/replacement limit.
- Whether the first-ten route can validly assess a longitudinal patient-specific method despite incomplete individuals.
- Whether an independent researcher can use a personal email rather than a university/institutional email.
- The precise embargo treatment of participant-derived assemblies and detailed logs.
- The official result-notification date and final metric set beyond the portal’s current MetaQUAST 5.2.0 statement.

These questions are captured in `CAMI_SUPPORT_INQUIRY_DRAFT.md`. Until answered, submission structure remains `PAUSE_BLOCKED`.
