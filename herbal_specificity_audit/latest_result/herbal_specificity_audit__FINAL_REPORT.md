# Does the disease label matter in herbal network pharmacology?

## Final computational audit report

**Analysis date:** 2026-08-16  
**Mode:** deterministic pilot  
**Data source:** Europe PMC Open Access full-text XML  
**Primary unit:** one article-level, author-reported hub/core/key target-gene set  
**Author:** Byungwoong Yoo, Independent Researcher

## Executive finding

Author-reported target sets retained a statistically detectable disease-specific signal under both label permutation and the frequency-preserving incidence null. The magnitude and recoverability results determine whether that signal is practically useful.

The conclusion is deliberately bounded. This study evaluates whether the **reported computational target lists** carry recoverable information about the named disease. It does not determine whether an herb is clinically effective, whether a listed target is causal, or whether any individual formula is safe.

## Frozen question and estimand

For article `i`, let `G_i` be the explicit hub/core/key gene set and `Y_i` its disease label. The primary statistic was:

```text
Delta_J = mean[J(G_i,G_j) | Y_i = Y_j] - mean[J(G_i,G_j) | Y_i != Y_j]
J(A,B) = |A intersect B| / |A union B|
```

A positive, non-random `Delta_J` is the minimum expected signature of disease specificity.

## Corpus flow

| Stage | Count |
|---|---:|
| Europe PMC query hits | 4,427 |
| Metadata after deduplication | 4,427 |
| Selected for XML | 500 |
| XML downloaded | 499 |
| Explicit target lists extracted | 375 |
| High-confidence lists | 207 |
| Specific disease mapped | 220 |
| Primary specific-label analysis | 43 |
| Specific disease classes in primary analysis | 6 |

The frozen query was:

```text
OPEN_ACCESS:Y AND FIRST_PDATE:[2015-01-01 TO 2026-08-16] AND TITLE:"network pharmacology" AND (decoction OR formula OR herb OR herbal OR phytochemical OR "traditional Chinese medicine" OR "traditional medicine" OR natural OR plant)
```

Reviews, meta-analyses, bibliometric studies, editorials, protocols, and articles without an explicit finite target list were excluded by predeclared rules.

## Primary result: specific disease labels

| Quantity | Result |
|---|---:|
| Articles | 43 |
| Disease classes | 6 |
| Within-disease mean Jaccard | 0.1183 |
| Between-disease mean Jaccard | 0.0899 |
| `Delta_J` | 0.0284 |
| Stratified-bootstrap 95% CI | [0.0995, 0.1987] |
| Label permutations | 1000 |
| One-sided permutation p | 0.0030 |
| Two-sided permutation p | 0.0030 |

### List-size and publication-era stratified permutation

| Quantity | Result |
|---|---:|
| Stratified `Delta_J` | 0.0284 |
| One-sided p | 0.0020 |
| Null 95% interval | [-0.010668629977701582, 0.018434086172299986] |

### Frequency-preserving incidence null

This null randomizes article–gene edges while preserving every article's list length and every gene's total article frequency. It therefore asks whether disease grouping explains more overlap than the generic hub-frequency structure alone.

| Quantity | Result |
|---|---:|
| Observed `Delta_J` | 0.0284 |
| Randomized incidence matrices | 50 |
| Null mean | 0.0019 |
| Null 95% interval | [-0.008292502166206521, 0.01416925165720156] |
| One-sided p | 0.0196 |

## Can an AI recover the disease from the target genes alone?

A linear support-vector classifier received only the article-by-gene binary matrix. TF–IDF weighting was fitted inside the cross-validation pipeline. No title, herb name, journal, year, pathway name, or abstract text was supplied.

| Metric | All genes | After removing top 10 recurrent genes |
|---|---:|---:|
| Articles | 38 | 38 |
| Classes | 5 | 5 |
| Accuracy | 0.4286 | 0.2988 |
| Balanced accuracy | 0.4333 | 0.2867 |
| Macro-F1 | 0.3587 | 0.2305 |
| Balanced-accuracy chance | 0.2000 | 0.2000 |
| Label-permutation p | 0.0645 | 0.0968 |

Classification above chance would indicate some disease information, but not necessarily biological specificity: journals, databases, formula traditions, or repeated pipelines can create a learnable signature. Failure to exceed chance is stronger evidence that the reported lists are label-invariant.

## Recurrent generic targets

| Rank | Gene | Articles | Article share | Disease labels | Normalized disease entropy |
|---:|---|---:|---:|---:|---:|
| 1 | AKT1 | 100 | 48.3% | 34 | 0.858 |
| 2 | TNF | 77 | 37.2% | 30 | 0.837 |
| 3 | IL6 | 68 | 32.9% | 28 | 0.825 |
| 4 | EGFR | 65 | 31.4% | 24 | 0.773 |
| 5 | TP53 | 61 | 29.5% | 26 | 0.804 |
| 6 | ESR1 | 54 | 26.1% | 23 | 0.765 |
| 7 | PTGS2 | 53 | 25.6% | 20 | 0.731 |
| 8 | STAT3 | 47 | 22.7% | 19 | 0.699 |
| 9 | MAPK1 | 46 | 22.2% | 22 | 0.759 |
| 10 | SRC | 43 | 20.8% | 14 | 0.642 |
| 11 | CASP3 | 42 | 20.3% | 18 | 0.725 |
| 12 | MMP9 | 40 | 19.3% | 19 | 0.738 |
| 13 | JUN | 40 | 19.3% | 20 | 0.732 |
| 14 | IL1B | 35 | 16.9% | 17 | 0.719 |
| 15 | VEGFA | 34 | 16.4% | 17 | 0.719 |
| 16 | HSP90AA1 | 33 | 15.9% | 14 | 0.646 |
| 17 | MAPK3 | 27 | 13.0% | 14 | 0.657 |
| 18 | PPARG | 24 | 11.6% | 10 | 0.579 |
| 19 | MAPK8 | 23 | 11.1% | 9 | 0.538 |
| 20 | HIF1A | 21 | 10.1% | 9 | 0.539 |

Summary concentration measures:

- Total article–gene edges: **1722**
- Unique genes: **327**
- Top-10 edge share: **35.7%**
- Top-20 edge share: **54.2%**
- Gene-frequency Gini: **0.703**
- Genes reported in at least 10% of articles: **20**

High frequency plus high disease entropy identifies a gene that is repeatedly reported across many unrelated disease labels—the pattern expected for a generic network hub rather than a disease fingerprint.

## Disease representation

| Disease label | Articles |
|---|---:|
| hepatocellular_carcinoma | 10 |
| depression | 9 |
| lung_cancer | 7 |
| breast_cancer | 6 |
| covid_19 | 6 |
| colorectal_cancer | 5 |

## What was solved, and what was not

### Solved by this audit

The analysis provides a reproducible answer to the following bounded question:

> Do explicit author-reported hub/core/key gene lists from open-access herbal network-pharmacology papers contain a statistically recoverable disease-specific signature under label permutation, list-size/year stratification, recurrent-hub removal, and an article–gene frequency-preserving null?

Every included article, target symbol, source span, exclusion reason, code path, random seed, and file checksum is in the accompanying package.

### Not established

- clinical efficacy or inefficacy of any herb or formula;
- causal involvement of a target;
- biochemical binding;
- product identity, dose, bioavailability, or safety;
- validity of network pharmacology outside the frozen corpus and extraction contract.

## Limitations

1. Target lists were machine-extracted from XML. The primary set used strict, predeclared high-confidence rules, but machine extraction can still miss lists embedded only in images or supplementary files.
2. Disease labels came from a frozen title lexicon. Unmapped or multiply framed diseases were omitted from the specific-label analysis.
3. The study evaluates what articles explicitly reported, not every target generated internally by their pipelines.
4. Open-access indexing and title wording create coverage selection.
5. Repeated use of the same public databases can be a property of the field rather than misconduct by individual authors.
6. A non-specific target list does not imply that the intervention has no biological or clinical effect.

## Reproducibility files

- `tables/articles.csv`: one row per extracted article.
- `tables/article_gene_edges.csv`: one row per article–gene edge.
- `tables/exclusions.csv`: all recorded exclusions.
- `tables/extraction_candidates.csv`: candidate lists and exact source spans.
- `results/results.json`: complete machine-readable results.
- `results/flow.json`: corpus flow.
- `figures/`: publication-ready plots.
- `MANIFEST.json`: SHA-256 hashes and sizes.
- `audit.py`: complete executable pipeline.

## Final bounded conclusion

Author-reported target sets retained a statistically detectable disease-specific signal under both label permutation and the frequency-preserving incidence null. The magnitude and recoverability results determine whether that signal is practically useful.
