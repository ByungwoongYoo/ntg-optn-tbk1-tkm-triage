# Prospective human missense VUS validation protocol

**Protocol freeze date:** 2026-08-17  
**Status:** protocol frozen; prediction cohort not yet generated  
**Purpose:** test future, not retrospective, clinical-variant generalization.

## Question

For variants classified as *uncertain significance* at the freeze date, can a fully frozen missense predictor assign calibrated pathogenicity probabilities that outperform frozen component predictors when later independent evidence resolves the classification?

## Important data constraint discovered before cohort construction

The ProteinGym v1.3 clinical score matrix used in the retrospective benchmark contains only 32,000 Pathogenic and 30,727 Benign rows. It contains no VUS rows. Therefore the ProteinGym matrix alone cannot be used to create a prospective VUS panel. A new score-generation or all-variant score source is required before predictions can be frozen.

## Cohort inclusion

At the eventual prediction-freeze date, include variants that satisfy all of the following:

1. ClinVar molecular consequence is missense SNV.
2. Reference assembly, chromosome, position, reference/alternate allele, transcript and protein consequence can be normalized without ambiguity.
3. Aggregate ClinVar classification is VUS.
4. No Pathogenic/Likely pathogenic or Benign/Likely benign submission is present at baseline.
5. The frozen predictor and every prespecified comparator produce a score.
6. Variant is not in the retrospective ProteinGym clinical training/evaluation set.

Exclude conflicting classifications, somatic-only records, variants lacking an unambiguous protein consequence, and records whose current classification relies on the same computational model as the primary evidence source in a way that would create direct circularity.

## Frozen prediction objects

For every included variant, publish before outcomes are known:

- normalized genomic and protein identifier;
- scores from every constituent predictor;
- fixed-ensemble score;
- dynamic-gating score, when the family-held-out experiment supports it;
- calibrated probability;
- abstain/non-abstain decision;
- reason for abstention;
- software version, source release, complete code and SHA-256 manifest.

## Future outcome definition

Primary outcome is a post-freeze reclassification to Pathogenic/Likely pathogenic or Benign/Likely benign supported by one of:

1. ClinGen expert panel or practice guideline;
2. ClinVar review status of at least two stars with no conflict and documented assertion criteria;
3. an independently published, well-validated multiplexed functional assay, analysed separately from the clinical endpoint.

Variants that remain VUS are not counted as benign and do not enter the primary discrimination analysis.

## Primary metrics

- Brier score and log loss;
- calibration intercept and slope;
- AUROC and AUPRC;
- likelihood ratios in prespecified score bins;
- coverage-risk curve with abstention;
- performance stratified by gene, inheritance mechanism, review status and evidence source.

## Success criterion

A prospective advance requires all of the following:

1. at least 200 independently resolved variants;
2. lower Brier score and log loss than the frozen best component;
3. AUROC gain with 95% CI above zero or equivalent superiority in a prespecified proper scoring rule;
4. calibration slope between 0.8 and 1.2 after no post-outcome refitting;
5. at 80% coverage, high-confidence classification error below 10%;
6. the direction is retained in expert-panel-only and functional-assay-only sensitivities.

## Failure criterion

The prospective claim fails if the primary proper scoring rule does not improve, calibration deteriorates materially, superiority disappears in expert-panel-only records, or the result depends on computational evidence that helped create the future ClinVar classification.

## Current boundary

This protocol does not claim that the future validation has occurred. It fixes the rules needed to prevent later outcome-driven cohort selection, threshold tuning or model replacement.
