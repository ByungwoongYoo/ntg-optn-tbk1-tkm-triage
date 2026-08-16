# Does the Disease Label Matter? A Corpus-Scale Specificity Audit of Herbal Network Pharmacology

**Byungwoong Yoo**  
Independent Researcher, Republic of Korea

## Abstract

**Background:** Herbal network-pharmacology studies commonly interpret hub or core targets as disease mechanisms, but the disease specificity of those target sets has not been tested at corpus scale.  
**Methods:** We froze a Europe PMC Open Access query through 2026-08-16, extracted explicit author-reported hub/core/key target lists from article XML using an NCBI human-gene dictionary, and mapped diseases from titles with a predeclared lexicon. The primary statistic compared mean within-disease and between-disease Jaccard similarity. Significance was assessed by disease-label permutation, a list-size/publication-era-stratified permutation, and a bipartite incidence null preserving article list lengths and gene frequencies. A cross-validated linear classifier tested whether disease labels could be recovered from gene sets alone.  
**Results:** The primary analysis included 43 articles across 6 repeatedly represented disease labels. Mean Jaccard similarity was 0.1183 within diseases and 0.0899 between diseases, for a difference of 0.0284 (one-sided permutation p=0.0030). The top 10 genes accounted for 35.7% of all article–gene edges. Full null-control and classification results are reported in the accompanying tables.  
**Conclusions:** The result is restricted to the specificity of reported computational target lists and does not address clinical efficacy or causal mechanism. The corpus-level tests quantify whether disease labels leave a recoverable signature after controlling for generic hub reuse.

## Introduction

Network pharmacology is widely used to propose multi-component, multi-target mechanisms for herbal medicines. A recurring analytic sequence intersects predicted compound targets with disease-associated genes, constructs a protein–protein interaction network, and selects high-degree or central nodes as hub targets. The resulting targets are often discussed as disease-specific mechanistic explanations. However, highly annotated and highly connected genes can recur across many diseases, and the same databases and centrality algorithms may reproduce a generic set of inflammatory, apoptotic, and proliferative hubs.

This study asks a deliberately falsifiable question: if disease labels are removed, do the reported target lists retain enough structure to identify the disease? A disease-specific target set should be more similar to other studies of the same disease than to studies of different diseases, should outperform label-permuted controls, and should retain signal after article list length and global gene frequency are controlled.

## Methods

### Protocol and corpus

The query, date window, inclusion rules, disease lexicon, target-list anchors, random seed, and primary estimand were fixed in code before full-corpus outcome inspection. Europe PMC Open Access full-text XML was searched with the query reproduced in the final report. Reviews, meta-analyses, bibliometric studies, protocols, and papers without a finite explicit target list were excluded.

### Target extraction

Paragraphs, table cells, and captions containing predeclared phrases such as “hub genes,” “core targets,” and “key genes” were searched. Candidate gene symbols were normalized against the NCBI Homo sapiens gene information file. The primary analysis retained high-confidence candidates defined by an explicit count match or a dense gene list in an abstract, Results, Conclusion, or table context. The highest-scoring high-confidence set was selected per article. Exact source spans were retained for audit.

### Disease labels

Disease labels were assigned from article titles by a frozen, longest-match lexicon. Labels with fewer than the predeclared minimum number of articles were excluded from disease-level analysis but remained in corpus-wide hub-frequency summaries.

### Statistical analysis

Jaccard similarity was calculated for every article pair. The primary effect was the mean within-disease similarity minus the mean between-disease similarity. Label permutation preserved target lists and disease-group sizes. A stratified permutation operated within list-size quartiles and publication-era bins. A bipartite double-edge-swap null preserved each article's target-list length and each gene's total frequency. Sensitivity analyses removed the 5, 10, and 20 most frequently reported genes. Disease recoverability was evaluated with repeated stratified cross-validation using TF–IDF-weighted binary gene features and a class-weighted linear support-vector classifier. All tests used seed 20260816.

## Results

See `FINAL_REPORT.md`, `results/results.json`, and the accompanying figures and tables. The numerical outputs in the Abstract were generated directly from the frozen analysis object, not entered manually.

## Discussion

The central interpretation depends on the direction and robustness of the corpus-level signal. A null result indicates that the reported hub/core/key lists do not constitute a reproducible disease fingerprint under the tested controls. A positive result indicates some disease-associated structure, but it remains necessary to distinguish biological specificity from shared databases, repeated pipelines, and publication subcultures. Generic hub concentration is therefore interpreted jointly with label-permutation, frequency-preserving randomization, hub-removal sensitivity, and out-of-sample disease recoverability.

The study does not evaluate whether any herbal intervention is effective. It evaluates whether a common computational evidence object—the reported target list—supports the disease-specific interpretation often attached to it.

## Data and code availability

All derived article-level data, source locators, exclusions, code, environment, random seeds, and checksums are included in the reproducibility package. Europe PMC full texts are not redistributed; the package contains identifiers and short extraction spans sufficient to locate each claim in the source article.
