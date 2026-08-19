# AMR unexplained-resistance project — research status

Last updated: 2026-08-19  
Repository branch: `amr-kp-colistin-portal-unitig-gwas-20260819`

## Objective

Identify, or rigorously fail to identify, a previously unrecognized genomic marker associated with colistin resistance in *Klebsiella pneumoniae* using public phenotype-linked genomes. A positive result must survive known-mechanism exclusion, population-structure correction, untouched source-held-out validation, cross-source/country/lineage checks, independent-database replication, and a current novelty audit.

## Current gate tracker

| Gate | State | Evidence / next action |
|---|---|---|
| Global pathogen–drug target ranking | Complete | *K. pneumoniae*–colistin selected from the frozen comparison set. |
| Phenotype and duplicate/conflict audit | Complete with limitation | No within-assembly R/S contradiction in the strict raw labels; explicit reference-BMD vocabulary is sparse, so method metadata remain a sensitivity-analysis limitation. |
| Portal-residual cohort construction | Complete | Resistant genomes not immediately explained by the prespecified AMR Portal marker mapping were isolated as an operational discovery pool; this is not a novelty claim. |
| Targeted colistin-pathway variant analysis | Running / audit | Corrected exact assembly-version handling and source-held-out evaluation workflow are in place. Results require independent review before interpretation. |
| Independent sequence-level known-mechanism audit | Running / repair | Acquired `mcr`, `mgrB` disruption, curated chromosomal changes and relevant regulatory loci must be removed or conditioned on. Kleborate/Kaptive compatibility and complete output integrity are explicit execution gates. |
| Whole-genome discovery-only unitig GWAS | Triggered | Discovery is restricted to the source-disjoint training subset; validation unitigs are queried only after selection. Mash/pyseer sensitivity models and cross-source/country replication are prespecified. |
| Within-lineage / leave-one-source sensitivity | Pending upstream candidates | Performed only on discovery-selected candidates; no post-hoc threshold relaxation. |
| Genomic-context and assembly-artifact audit | Pending upstream candidates | Candidate context, contig edge, low complexity, insertion sequence, annotation and recurrence will be checked. |
| Independent database replication | Pending | NCBI Pathogen Detection/BV-BRC or another non-overlapping phenotype-linked cohort after accession-level deduplication. |
| CARD/AMRFinderPlus/NCBI/literature novelty audit | Pending | Conducted only for candidates surviving statistical and mechanism gates. |
| Positive public claim | Not authorized | No candidate currently has permission to be described as a new marker, gene, or mechanism. |

## Decision rule

The analysis is optimized for validity, not for obtaining a positive p-value. Execution errors are repaired and scientifically justified sensitivity analyses are pursued, but a candidate is never promoted by changing the threshold, validation split, phenotype definition, or covariate set after seeing the desired result.

If no marker survives, the study closes as a rigorous negative audit or advances the next pathogen–drug pair under a separately frozen protocol. If a marker survives all computational gates, it is described only as an independently replicated resistance-associated marker until functional laboratory validation establishes causality.
