# Panax publication-readiness decision

## Decision

The analysis is worth continuing, but it is not submission-ready. The defensible article type is a reproducible computational virus-discovery Original Article or Short Report based on public transcriptome reanalysis. The present evidence does not support a report of three novel viruses or species.

Working title:

> Partial Picornavirales-like RdRP sequence clusters recovered from public *Panax notoginseng* root transcriptomes: a reproducible computational audit

Working claim allowed only if the enhanced fragment, mapping-seed-free de novo,
and current-database gates retain all three clusters:

> Three divergent partial RdRP-containing sequence clusters with Picornavirales-like homology were recovered and raw-read-supported in public *Panax notoginseng*-associated root RNA-seq libraries.

Current conservative statement: three predefined clusters received support in
the legacy candidate-panel mapping audit; the enhanced retained count is
pending and may be lower.

The three labels denote computational sequence clusters. A1 and A2 are not counted as two virus species: their selected 408-aa RdRP segments are 82.353% identical, and the available partial region is not a taxonomic delimiter.

## Claims excluded

- three novel viruses or species
- viruses infecting *Panax notoginseng*
- healthy-versus-root-rot association or causation
- active replication, pathogenicity, or transmission
- complete genomes
- formal genus, family, or species assignment

## P0 submission blockers

1. Complete and freeze the current nr, nt, TSA, environmental, viral, and nonviral search gate with validated controls, search dates, parameters, and database provenance.
2. Perform mapping-seed-free per-library assembly and recover non-source assemblies. A2 requires an independent second assembly; if that fails, downgrade or remove it from the primary discovery count.
3. Add proper-pair, insert/orientation, duplicate-aware, candidate-diagnostic read, endpoint, soft-clip, chimera, and coverage-discontinuity checks.
4. Compete reads against *Panax* nuclear and organellar sequences, known *Panax* viruses, fungal/oomycete/root-community and other plausible eukaryotic decoys, plus UniVec/adapters.
5. Treat the `Zh/Bh` phenotype crosswalk as unresolved unless an authoritative mapping is obtained. If it remains unresolved, omit every healthy/root-rot comparison. Retain the HiSeq 4000 versus NovaSeq 6000 platform discrepancy.
6. Extend sequence architecture where possible to Pro-Pol, capsid proteins, additional ORFs, and putative segments. Rebuild the homolog panel with current top hits and disclose B's unstable tree context.
7. Obtain a durable data-release plan: permanent DOI archive for evidence/code and an appropriate INSDC/TSA/TPA accession route for public-read-derived sequences. Expiring GitHub Actions artifacts are not a Data Availability solution.

## P1 analyses

- A1/A2 haplotype, assembly-bubble, and recombination sensitivity
- broad dereplicated current-homolog Pro-Pol/RdRP tree
- trimming, model, taxon-sampling, and B-placement sensitivity
- capsid/structural-segment search and segment co-occurrence
- support-threshold sensitivity table rather than binary support alone
- pool-level effect sizes only if an authoritative phenotype crosswalk is obtained; no confirmatory root-rot association test is supportable with the current pooled 3-versus-3 design

## Candidate disposition

| Candidate | Current disposition | Why | Publication condition |
|---|---|---|---|
| A1 | Legacy retain; enhanced pending | Strong PF00680 evidence and multi-run mapping support; an older artifact contains a related RdRP, but that assembly has not yet passed the mapping-seed-free per-library gate | Resolve proper-pair/duplicate sensitivity, obtain non-source mapping-seed-free assembly context, and explain cross-run sequence heterogeneity. |
| A2 | Provisional | Broadest mapping support and computationally separable from A1 | Must obtain a credible independent non-source assembly; otherwise downgrade or exclude from the primary count. |
| B | Legacy strongest mapping candidate; enhanced pending | Longest partial ORF and near-full high-identity non-source **masked-consensus mapping** reproduction; this is not yet a de novo recovery claim | Complete de novo, decoy, and proper-pair checks and broaden architecture/homolog analyses because phylogenetic context is unstable. |

## Tables required

1. Run, sample, archive token, authoritative phenotype status, pooling, platform, and FASTQ provenance.
2. Candidate lengths, ORF/RdRP coordinates, open boundaries, adapter handling, sequence hashes, and independent assemblies.
3. Search database, date/version, top hit, identity, coverage, E-value, controls, and technical completeness.
4. Per-run proper unique fragments, breadth at 1x/5x/10x, median depth, diagnostic reads, and duplicate sensitivity.
5. Domain, motif, architecture, completeness, Pro-Pol/capsid evidence, and taxonomic claim boundary.
6. Decoy, threshold, alignment-trimming, model, and tree-sensitivity results.

## Figures required

1. Discovery, validation, downgrade, and exclusion flow.
2. Candidate architecture and per-run coverage tracks.
3. Current-homolog maximum-likelihood tree with support values and a B sensitivity panel.
4. Run-token support heatmap. A phenotype plot is allowed only if the crosswalk is authoritatively resolved.

## APC-free journal strategy

All zero-APC recommendations below refer to the conventional subscription route, not optional open access.

1. **Virus Genes** — preferred after the blockers are closed. NGS, transcriptomic/metagenomic, and phylogenetic work is in scope, but the study must be complete and sequence accession/deposition is expected. <https://link.springer.com/journal/11262/aims-and-scope> and <https://link.springer.com/journal/11262/how-to-publish-with-us>
2. **VirusDisease** — practical fallback for plant-virus characterization under the same claim boundary and deposition requirements. <https://link.springer.com/journal/13337/aims-and-scope> and <https://link.springer.com/journal/13337/how-to-publish-with-us>
3. **Journal of Plant Pathology** — conditional fallback only if the archive phenotype crosswalk and plant-pathology relevance are resolved. <https://link.springer.com/journal/42161/aims-and-scope> and <https://link.springer.com/journal/42161/how-to-publish-with-us>
4. **Archives of Virology** — not recommended for the present partial candidates; its sequence-only route expects complete or fundamentally distinct genome organization with biological significance. <https://link.springer.com/journal/705/aims-and-scope>
5. **bioRxiv** — zero-fee preprint option after the P0/P1 work is complete and the candidate wording is locked. <https://www.biorxiv.org/about/FAQ>

Taxonomic architecture decisions should follow the current ICTV family criteria rather than a single RdRP-core tree: <https://ictv.global/report/chapter/secoviridae/secoviridae>.

Potential public-read-derived sequence submission routes should be confirmed with INSDC before submission: <https://www.ncbi.nlm.nih.gov/genbank/submit_types/>.
