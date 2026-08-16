# Does the disease label matter in herbal network pharmacology?

## Authoritative final v3 held-out corpus report

**Analysis date:** 2026-08-16  
**Data source:** Europe PMC Open Access full-text XML  
**Primary unit:** one article-level, author-reported hub/core/key target-gene set  
**Author:** Byungwoong Yoo, Independent Researcher  
**Prespecified decision:** `no_reliable_positive_specificity_signal`

## Final conclusion

The held-out corpus did not show a reliable positive within-disease similarity advantage. Under the frozen corpus, strict extraction rule, disease labels, and prespecified controls, the reported hub/core/key target lists did **not** provide a reproducible disease-specific fingerprint.

This conclusion concerns the disease specificity of a published computational evidence object. It does not establish clinical efficacy or inefficacy, causal mechanism, biochemical binding, product identity, dose, bioavailability, or safety.

## Frozen question

For article `i`, let `G_i` be its explicit author-reported target-gene set and `Y_i` its title-mapped disease label. The primary estimand was:

```text
Delta_J = mean[J(G_i,G_j) | Y_i = Y_j]
        - mean[J(G_i,G_j) | Y_i != Y_j]

J(A,B) = |A intersect B| / |A union B|
```

A reproducible disease fingerprint requires a positive effect that exceeds the label-permutation null, has a bootstrap interval above zero, and survives the prespecified robustness controls.

## Development and held-out separation

The deterministic first 500 metadata-eligible records were used only to refine the extraction lexicon, remove false-positive algorithm abbreviations, and repair the confidence-interval procedure. They were excluded from every primary numerical result below.

## Corpus flow

| Stage | N |
|---|---:|
| Europe PMC query hits | 4,427 |
| Selected for full-text XML | 2,881 |
| XML downloaded | 2,879 |
| Explicit target lists extracted | 2,152 |
| Strict exact-count lists in the complete run | 650 |
| Held-out extracted target lists | 1,776 |
| Held-out strict exact-count lists | 519 |
| Primary specific-label analysis | 114 |
| Disease classes | 9 |

The frozen query was:

```text
OPEN_ACCESS:Y AND FIRST_PDATE:[2015-01-01 TO 2026-08-16]
AND TITLE:"network pharmacology"
AND (decoction OR formula OR herb OR herbal OR phytochemical
OR "traditional Chinese medicine" OR "traditional medicine" OR natural OR plant)
```

Reviews, meta-analyses, bibliometric studies, editorials, protocols, and records without an explicit finite target list were excluded by coded rules.

## Extraction checks

| Check | Result |
|---|---:|
| Held-out extracted articles | 1,776 |
| Held-out strict lists | 519 |
| Strict share | 29.22% |
| Author-stated count matched extracted list | 100% by strict inclusion rule |
| Independent official-symbol-only exact agreement | 505/519 (97.30%) |

The 100% author-count figure is an eligibility condition, not an independent accuracy estimate. The official-symbol-only extractor is the independent computational agreement check. Its 14 discordant articles were removed in a separate sensitivity analysis.

## Primary result

| Quantity | Result |
|---|---:|
| Within-disease mean Jaccard | 0.1219115689 |
| Between-disease mean Jaccard | 0.1214275645 |
| `Delta_J` | 0.0004840045 |
| Relative within-vs-between increase | 0.3986% |
| Within-disease pairs | 722 |
| Between-disease pairs | 5,719 |
| Label permutations | 5,000 |
| One-sided permutation p | 0.4335 |
| Two-sided permutation p | 0.9186 |
| Stratified article-bootstrap 95% CI | [-0.0124, 0.0140] |

The arithmetic is:

```text
0.1219115689 - 0.1214275645 = 0.0004840045
0.0004840045 / 0.1214275645 = 0.003986 = 0.3986%
```

The observed advantage was therefore effectively zero and statistically indistinguishable from random reassignment of disease labels.

## Prespecified controls

| Analysis | Delta | One-sided p |
|---|---:|---:|
| List-size/publication-era-stratified permutation | 0.0005 | 0.4531 |
| Frequency-preserving article-gene incidence null | 0.0005 | 0.9003 |
| Dice similarity | 0.0001 | 0.4661 |
| Cosine similarity | -0.0009 | 0.5225 |
| Independent official-symbol agreement subset | -0.0007 | 0.5243 |
| Exact gene-set deduplication | 0.0005 | 0.4337 |
| Remove 5 most recurrent genes | -0.0026 | 0.6892 |
| Remove 10 most recurrent genes | 0.0019 | 0.3038 |
| Remove 20 most recurrent genes | 0.0013 | 0.3403 |

The frequency-preserving null retained every article's list length and every gene's article frequency. Its null mean was 0.0041, compared with the observed 0.0005. Thus the observed disease grouping explained no more overlap than the generic hub-frequency structure.

## Disease recoverability from genes alone

A class-weighted linear SVM received only the binary article-by-gene matrix, with TF-IDF fitted inside repeated stratified cross-validation. It received no title, herb, formula, journal, year, abstract, or pathway information.

| Metric | All genes | Top 10 recurrent genes removed |
|---|---:|---:|
| Articles | 114 | 114 |
| Classes | 9 | 9 |
| Accuracy | 0.1517 | 0.1305 |
| Balanced accuracy | 0.1561 | 0.1361 |
| Macro-F1 | 0.1358 | 0.1177 |
| Balanced-accuracy chance | 0.1111 | 0.1111 |
| Label-permutation p | 0.0796 | 0.2886 |

The gene lists did not support statistically significant recovery of the original disease label. The small apparent signal weakened after the most recurrent hubs were removed.

## Generic hub concentration

Across the 519 held-out strict lists:

- article-gene edges: 5,016;
- unique genes: 477;
- top-10 edge share: 37.7%;
- top-20 edge share: 57.5%;
- gene-frequency Gini: 0.804;
- genes appearing in at least 10% of articles: 23.

The most recurrent genes were AKT1 (265/519, 51.1%), TNF (231/519, 44.5%), IL6 (226/519, 43.5%), TP53 (198/519, 38.2%), and EGFR (196/519, 37.8%). Each was distributed across many unrelated disease labels.

## What this analysis resolves

The following bounded empirical question is closed for this frozen public corpus and extraction contract:

> Do explicit author-reported hub/core/key lists carry a reproducible disease-specific signature after disease-label permutation, list-size/year stratification, recurrent-hub removal, an independent symbol-extraction subset, exact-set deduplication, and an article-gene frequency-preserving null?

**Answer: no reliable positive disease-specific signature was detected.**

## What it does not resolve

- Whether a particular herb or formula is clinically effective.
- Whether a particular listed target is causal.
- Whether experimentally measured, tissue-specific, patient-derived, or multi-omics networks retain disease specificity.
- Whether network pharmacology outside the frozen open-access title query behaves identically.
- Whether an individual paper committed misconduct.

## Reproducibility record

- Successful GitHub Actions run: `31949335813`
- Job: `95170261504`
- Run commit: `93f099ecc0112de254f5c312783baf272b07bc7f`
- Artifact: `9264344939`
- Artifact SHA-256: `8f73a70289338b6ed53392808124c95a3353c0a3359135a6e33a5c7f5678a769`
- Base random seed: `20260816`

The complete artifact contains executable code, frozen protocol, derived article and article-gene tables, exact source spans, exclusions, results, figures, manifests, and file hashes. Raw Europe PMC XML is identified by PMCID and checksum but is not redistributed.

## Press-safe statement

> In a held-out audit beginning with 4,427 open-access records, reported hub/core/key gene lists from 114 strictly eligible studies across nine diseases were no more similar within the same disease than across different diseases, and an AI classifier could not reliably recover the disease from those genes alone.

Do not convert this result into claims that herbal medicines are ineffective, that all network pharmacology is invalid, or that every screened paper was wrong.
