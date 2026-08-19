# K. pneumoniae–colistin discovery: fixed decision protocol

Date frozen: 2026-08-19  
Study lead: Byungwoong Yoo  
Scope: public-data genomic association study; no experimental causal claim

## Claim boundary

The strongest claim available from this workflow, before laboratory validation, is:

> A previously unrecognized genomic marker associated with colistin resistance in *Klebsiella pneumoniae*, independently replicated after population-structure control and not immediately explained by the audited known mechanisms.

The workflow does **not** establish a new resistance gene, a causal mechanism, a clinical diagnostic, or treatment guidance. No positive result may be manufactured by relaxing gates after seeing outcomes.

## Fixed promotion gates

A candidate is promoted only if all applicable gates are passed:

1. **Phenotype integrity**
   - no unresolved R/S conflict for the same assembly;
   - binary label provenance retained;
   - MIC/method/standard metadata audited where available;
   - sensitivity analysis excluding ambiguous reinterpretations.

2. **Discovery significance**
   - discovered without access to the source-held-out validation cohort;
   - population structure controlled with Mash/pyseer sensitivity models;
   - positive direction in all prespecified structure dimensions;
   - multiple-testing-adjusted significance or a prespecified discovery-screen threshold followed by independent validation.

3. **Untouched validation**
   - same marker queried in a BioProject/source-held-out cohort;
   - effect direction concordant;
   - odds ratio greater than 1 with 95% confidence interval excluding 1;
   - validation multiplicity controlled across discovery-selected markers.

4. **Transportability**
   - positive random-effects summary across at least three source groups;
   - positive random-effects summary across at least three countries where metadata permit;
   - leave-one-source and leave-one-country sensitivity checks;
   - no single study, country, sequence type, or outbreak clone explains the complete signal.

5. **Lineage robustness**
   - within-lineage or matched-lineage analysis where sample size permits;
   - association remains after sequence-type/Mash-cluster adjustment;
   - genomic relatedness and sample duplication audited.

6. **Known-mechanism exclusion**
   - independent sequence-level audit for acquired `mcr` genes;
   - `mgrB` truncation, frameshift, deletion, promoter/intergenic disruption, and insertion-sequence events;
   - coding and relevant regulatory changes in `pmrAB`, `phoPQ`, `crrAB`, `arn/eptA` and other curated colistin pathways;
   - Kleborate, AMRFinderPlus, CARD/RGI, ResFinder/PointFinder or equivalent results retained with versions;
   - promoted association re-estimated after removing genomes with a plausible known mechanism.

7. **Assembly and context integrity**
   - assembly quality, species assignment, contamination and anomalous length audited;
   - candidate supported by coherent genomic context, not contig-edge or low-complexity artifact;
   - sequence context recurs independently and maps consistently;
   - unitigs representing the same locus collapsed before counting independent discoveries.

8. **Novelty audit**
   - exact and relaxed sequence searches against current NCBI nucleotide/protein resources;
   - CARD, AMRFinderPlus reference catalogues and relevant curated databases searched;
   - locus/gene/mutation and mechanistic literature searched;
   - prior reports classified as identical, adjacent, mechanistically related, or genuinely unreported;
   - novelty wording revised downward whenever an earlier report is found.

9. **Independent-database replication**
   - candidate evaluated in a second database or independently assembled cohort after accession-level deduplication;
   - discovery records are not silently recycled as validation records;
   - failure to obtain an independent cohort is disclosed and prevents the strongest claim.

## Prespecified pivots if the primary analysis is negative

Pivots are scientific sensitivity analyses, not routes to force a positive result. Each pivot must preserve a fresh holdout and multiplicity control.

### P1 — phenotype harmonization
Restrict to isolates with explicit MIC and accepted broth-dilution terminology; harmonize breakpoint-year/standard where metadata allow. This may reduce power but improves label specificity.

### P2 — lineage-stratified analysis
Run within the largest adequately powered sequence types/Mash clusters, followed by a meta-analysis across lineages. This targets convergent signals hidden by between-lineage structure.

### P3 — gene-level rare-variant burden
Aggregate disruptive and rare nonsynonymous variants by gene or regulatory region rather than testing each rare allele separately. Burden definitions must be fixed before validation.

### P4 — structural and intergenic variation
Audit insertion sequences, promoter disruptions, gene interruptions, contig breaks and local copy-number changes that SNP-only or coding-only analyses may miss, especially around `mgrB` and lipid-A regulatory loci.

### P5 — pangenome/intergenic unitigs
Use discovery-only unitigs and context mapping to find gene-presence, promoter, mobile-element or unannotated intergenic signals. Highly correlated unitigs are collapsed to loci before validation.

### P6 — epistasis and conditional effects
Test a small, biologically prespecified set of candidate-by-known-mechanism interactions only after the main-effect analysis. No exhaustive interaction fishing.

### P7 — independent pathogen–drug pair
If *K. pneumoniae*–colistin yields no surviving unexplained marker, return to the frozen global target-ranking table and advance the next eligible pair. The new pair receives a new discovery/validation split and a separate protocol; evidence is not pooled post hoc with this study.

## Stop rules

The study is closed as a rigorous negative result if:

- no candidate survives untouched validation;
- apparent signals disappear after lineage/source control;
- all surviving signals map to known mechanisms;
- independent-database replication fails;
- phenotype or assembly quality is insufficient to support the claim.

A negative result must not be relabelled as a discovery. It may still support a methodological report on how apparently unexplained AMR collapses under stringent audit.

## Public-release gate

A media/preprint package is released as a positive discovery only after gates 1–9 are documented in a candidate dossier containing:

- accession-level cohort manifest and exclusions;
- discovery and untouched-validation statistics;
- source/country/lineage sensitivity analyses;
- known-mechanism audit;
- sequence/context files and annotations;
- database/literature novelty table;
- reproducible code, environment, SHA-256 manifest and immutable release DOI;
- explicit statement that functional causality remains unproven without laboratory validation.

Until then, public wording is limited to cohort construction, methodological progress, or a negative audit result.
