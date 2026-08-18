# Public-data discovery of unexplained antimicrobial resistance markers

This isolated branch starts a reproducible search for genomic markers associated with phenotypic antimicrobial resistance that is **not matched to a known AMRFinderPlus determinant for the same antibiotic** in public AMR resources.

## Evidence boundary

The first-stage term `portal-unexplained resistance` means:

1. a strict resistant phenotype is present for an assembly–antibiotic pair;
2. no genotype record is joined to that same assembly and antibiotic ontology in the EMBL-EBI AMR Portal release used for analysis.

It does **not** yet mean that all known mechanisms have been biologically excluded. Ontology gaps, incomplete assemblies, expression mechanisms, promoter changes, loss-of-function mutations, porin/efflux effects, and phenotype noise remain possible. Any candidate marker must survive population-structure control, lineage-held-out validation, geographic/temporal replication, known-mechanism re-audit, and literature/database novelty checking before public release.

## Staged workflow

1. **Global census and target selection** using the latest stable EMBL-EBI AMR Portal phenotype/genotype Parquet release and NCBI Pathogen Detection scale/provenance records.
2. **Strict phenotype deduplication** and anti-join against known same-antibiotic AMRFinderPlus calls.
3. **Genome retrieval** for a pre-specified discovery/validation subset.
4. **Population-structure-aware bacterial GWAS** using unitigs, gene presence/absence, loss-of-function and nonsynonymous variants.
5. **Independent replication** across held-out lineages, countries, years, and source studies.
6. **Mechanistic plausibility and novelty audit** against AMRFinderPlus, CARD, ResFinder, PointFinder, NCBI literature and current primary studies.
7. **Conservative result package** with raw evidence, code, manifests, a technical report and a press-safe claim boundary.

No causal, clinical, diagnostic, treatment, or surveillance claim is made at project initialization.