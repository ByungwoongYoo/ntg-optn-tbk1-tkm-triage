# AMR residual marker discovery

Goal: identify reproducible whole-genome markers associated with phenotypic antimicrobial resistance among isolates whose phenotype is not matched by a catalogued AMRFinderPlus determinant.

Stages:

1. Residual-pair screen across a pinned EMBL-EBI AMR Portal release.
2. Whole-genome feature discovery only after a pair passes volume and known-determinant calibration gates.
3. Lineage-aware validation across held-out countries, years, projects, and genetic clusters.
4. Novelty audit against AMRFinderPlus, CARD, ResFinder/PointFinder, literature, and protein/domain databases.
5. Independent sequence and phenotype replication before causal language.

A resistant isolate with no exact ontology-matched AMRFinderPlus row is only a *candidate residual*. It is not automatically biologically unexplained.
