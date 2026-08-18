# Protocol freeze v1 — lineage-aware discovery of residual AMR markers

**Frozen:** 2026-08-19  
**Primary public release:** EMBL–EBI AMR Portal 2026-07  
**Independent replication universe:** NCBI Pathogen Detection / AST and MicroBIGG-E snapshots available at validation time  
**Status:** prospective computational protocol; no marker result existed when this document was committed

## 1. Scientific question

Among pathogen isolates with a resistant antimicrobial-susceptibility phenotype, which whole-genome sequence features are reproducibly associated with resistance after excluding isolates whose phenotype is already matched to a catalogued AMRFinderPlus determinant for the same antibiotic ontology?

The study is designed to discover **genetic marker candidates**, not to establish a causal resistance gene or a clinically deployable diagnostic.

## 2. Frozen claim boundary

A successful computational result may support only the following form of claim:

> A previously uncatalogued sequence marker was reproducibly associated with residual resistance in multiple lineage-aware, non-random validation sets.

The study will not claim that a marker causes resistance, confers resistance when transferred, is sufficient for resistance, predicts treatment outcome, or should alter patient care without independent laboratory validation.

A phenotype lacking an exact ontology-matched AMRFinderPlus row is termed a **candidate residual phenotype**. It is not automatically biologically unexplained because ontology mismatch, incomplete catalog coverage, AST heterogeneity, epistasis, gene expression, copy number, assembly error, and non-genetic mechanisms may remain.

## 3. Stage 1: organism–antibiotic pair selection

All eligible organism–antibiotic ontology pairs in the pinned release are screened before viewing whole-genome feature associations.

### 3.1 Phenotype mapping

Primary phenotype mapping:

- Resistant: `resistant` or `R`
- Susceptible: `susceptible`, `sensitive`, or `S`
- Intermediate and non-susceptible calls are excluded from the primary analysis

Sensitivity mapping additionally treats non-susceptible calls as resistant.

Duplicate rows are collapsed to one BioSample–assembly–species–antibiotic-ontology call. Groups containing discordant R/S calls are excluded from the primary analysis.

### 3.2 Calibration and eligibility gates

The conservative Stage-1 calibration subset requires an assembly to appear in the genotype release, demonstrating that at least one AMRFinderPlus result was generated for that assembly. An exact matching `antibiotic_ontology` row is treated as a catalogued determinant for that drug.

A pair is eligible only if it has:

- at least 200 unambiguous labeled assemblies;
- at least 50 resistant assemblies;
- at least 50 susceptible assemblies;
- at least 25 resistant assemblies without an exact ontology-matched determinant;
- at least 10 resistant assemblies with an exact ontology-matched determinant; and
- a Haldane–Anscombe corrected odds ratio of at least 2.0 for the catalogued determinant–resistance association.

Pairs are ranked by a frozen score combining residual-resistant count, susceptible count, country coverage, year coverage, and class balance. The top three eligible pairs, rather than only the top pair, advance when computationally feasible. This reduces dependence on a single idiosyncratic organism or drug.

### 3.3 Stage-1 no-go rule

The project does not proceed to a positive marker claim if no pair passes the calibration gates. Thresholds may be relaxed only in a separately versioned exploratory protocol; such results cannot be presented as confirmatory findings under this freeze.

## 4. Cohort construction for each advancing pair

### 4.1 Primary discovery contrast

- Cases: resistant assemblies without an exact ontology-matched AMRFinderPlus determinant
- Controls: susceptible assemblies without an exact ontology-matched AMRFinderPlus determinant

Susceptible assemblies with a known determinant are excluded from the primary analysis because they may represent silent genes, breakpoint differences, annotation error, or incomplete expression.

### 4.2 Mandatory metadata and QC

Each included assembly must have a retrievable genome annotation or genome sequence and a retrievable AMRFinderPlus report or an auditable no-hit report. Assemblies are excluded for failed retrieval, duplicate sequence, gross contamination, incompatible taxonomy, or prespecified assembly-quality failure.

Where available, the audit records assembly length, contig count, N50, ambiguous-base fraction, GC fraction, collection year, country, host, source, AST standard, AST method, and associated publication/project.

## 5. Non-random data partitioning

Random train/test splits are not accepted as the primary validation because closely related bacterial isolates can leak across folds.

For each pair, the following partitions are constructed before association testing:

1. **Discovery set:** older years and projects not assigned to the holdouts.
2. **Temporal holdout:** the latest adequately sized collection years.
3. **Geographic holdout:** one or more countries absent from discovery, chosen by deterministic sample-size rules.
4. **Project/publication holdout:** entire associated publication or sequencing project withheld.
5. **Lineage holdout:** one or more genetic clusters withheld when cluster sizes permit.

A holdout must contain at least 20 cases and 20 controls to serve as a formal replication set. Smaller partitions are reported descriptively only.

## 6. Genomic feature discovery

### 6.1 Primary feature representation

The primary genome-wide representation is compacted de Bruijn graph unitigs generated with a pinned `unitig-caller` container. Unitigs are chosen because they can represent genes, alleles, intergenic sequence, indels, and local structural variation without relying on a single reference genome.

### 6.2 Primary association model

The primary association is a fixed-effects microbial GWAS in `pyseer` with:

- binary residual-resistance phenotype;
- Mash-distance-derived population-structure covariates selected by a frozen scree/variance rule;
- country, collection period, AST standard, and major source category covariates when estimable without separation;
- minor allele frequency at least 1%;
- present in at least 10 assemblies;
- not present in more than 99% of assemblies; and
- significance threshold based on Bonferroni correction by the number of unique presence/absence patterns.

A mixed-model sensitivity analysis and lineage-effect analysis are mandatory. Features showing severe model non-convergence or complete separation without stable penalized estimates are not promoted.

### 6.3 Clustered locus reporting

Highly correlated significant unitigs are clustered into loci before interpretation. The primary reporting unit is a locus or haplotype block, not the number of individual overlapping unitigs.

## 7. Frozen replication endpoint

A locus qualifies as a **replicated residual-AMR marker candidate** only if all conditions below hold:

1. Discovery association passes the genome-wide corrected threshold.
2. Discovery odds ratio is at least 2.0 in the same direction.
3. The locus has at least 80% sequence identity over at least 80% of its query length when mapped across datasets, unless the primary feature is an exact short variant represented separately.
4. At least two independent formal holdouts show the same effect direction.
5. At least one holdout has nominal two-sided p < 0.05 after testing only the frozen discovery loci.
6. The pooled holdout odds ratio is at least 2.0 and its 95% confidence interval excludes 1.0.
7. Association remains directionally consistent in at least two major genetic lineages or passes a within-lineage meta-analysis without dominant single-lineage dependence.
8. Leave-one-country-out and leave-one-project-out analyses do not reduce the association to a single source.
9. The locus is not explained by a frozen negative-control artifact analysis.
10. Independent NCBI Pathogen Detection replication is directionally concordant when an adequately labeled, non-overlapping cohort exists.

A marker failing any condition may be reported as a discovery-only signal but not as the primary positive finding.

## 8. Artifact and falsification tests

Mandatory tests include:

- exact and class-level known-determinant remapping;
- strict-R versus broad-R phenotype definitions;
- CLSI-only and EUCAST-only subsets when sample sizes permit;
- MIC/zone-measurement reclassification where breakpoints and units are auditable;
- exclusion of low-quality and unusually fragmented assemblies;
- leave-one-country, leave-one-year, leave-one-project, and leave-one-lineage-out analyses;
- negative-control antibiotics with no plausible biological relationship;
- permutation within major lineage/project strata;
- feature missingness and assembly-contiguity association checks;
- proximity to contig ends, transposases, integrons, plasmid replicons, and known mobile elements;
- assessment of co-selection with known resistance loci for other drug classes; and
- duplicate/near-duplicate genome removal sensitivity.

## 9. Novelty audit

A replicated locus is called **previously uncatalogued** only after audit against the frozen or current versions recorded in the evidence manifest of:

- NCBI AMRFinderPlus database and organism-specific mutation catalogues;
- CARD/RGI, subject to redistributable-access constraints;
- ResFinder and PointFinder;
- NCBI nucleotide and protein databases;
- UniProt, Pfam, InterPro, and conserved-domain annotations;
- published organism–drug GWAS and resistance-mechanism literature; and
- the source publications associated with the contributing isolates.

Similarity to an uncharacterized protein is not evidence of a new resistance mechanism. A locus overlapping a known gene but containing a previously unreported allele is described as a **candidate resistance-associated allele**, not a new gene.

## 10. Independent replication and overlap control

The external replication cohort will be drawn from NCBI Pathogen Detection AST-linked isolates. BioSample IDs, assembly accessions, SRA runs, publication IDs, and near-identical genomes are used to remove overlap with the EBI discovery cohort.

The external dataset is not used to alter discovery thresholds, redefine the marker, select the favorable allele, or change the direction of effect.

## 11. Negative and null outcomes

The following are valid final outcomes and will be reported rather than hidden:

- no calibrated organism–drug pair;
- no genome-wide significant locus;
- discovery signal that fails non-random holdouts;
- lineage-confounded signal;
- association explained by a known determinant after deeper remapping;
- association explained by AST or assembly artifacts;
- replicated marker that is not novel after database/literature audit; or
- robust residual resistance with no single reproducible genetic marker.

## 12. Reproducibility and evidence preservation

The final evidence package must contain:

- exact public-data release labels and upstream checksums;
- all inclusion/exclusion manifests;
- frozen scripts and containers;
- software versions and Git commit SHAs;
- raw association outputs and model diagnostics;
- holdout definitions fixed before feature testing;
- locus sequences and coordinate mappings;
- database-search outputs and dates;
- an independent verifier for headline counts and effect estimates;
- a machine-readable claim-boundary file; and
- SHA-256 hashes for all released artifacts.

No press release, preprint, Zenodo record, or formal novelty claim is issued until the confirmatory gates in Sections 7–10 have been independently re-run from the frozen package.
