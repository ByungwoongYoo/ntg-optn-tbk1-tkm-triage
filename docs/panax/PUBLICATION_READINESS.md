# Panax publication-readiness decision

## Decision

The six-library reproducibility audit is complete and successful, but the study
is not yet submission-ready. All three predefined candidates (A1, A2, and B)
passed both the enhanced fragment-support gate and the mapping-seed-free de novo
recovery gate. The defensible article type remains a reproducible computational
sequence-discovery Original Article or Short Report based on public
transcriptome reanalysis. The evidence does not support a report of three novel
viruses or species.

Working title:

> Partial Picornavirales-like RdRP sequence clusters recovered from public *Panax notoginseng* root transcriptomes: a reproducible computational audit

Working claim allowed only after the current-database, competitive-decoy, and
deposition blockers below are closed:

> Three divergent partial RdRP-containing sequence clusters with Picornavirales-like homology were recovered and raw-read-supported in public *Panax notoginseng*-associated root RNA-seq libraries.

Current conservative statement:

> Three predefined partial RdRP-containing sequence clusters passed independent raw-read fragment-support and mapping-seed-free de novo recovery gates across six public *Panax*-associated root RNA-seq libraries.

This is a sequence-cluster reproducibility result. It is not evidence that the
clusters are three virus species, that they infect *Panax notoginseng*, or that
they are associated with disease.

The three labels denote computational sequence clusters. A1 and A2 are not counted as two virus species: their selected 408-aa RdRP segments are 82.353% identical, and the available partial region is not a taxonomic delimiter.

## Verified gate status

The enhanced six-library workflow run [`32482881928`](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32482881928)
completed successfully. Its aggregate artifact,
`PANAX_SIX_LIBRARY_ENHANCED_AUDIT_20260821` (`9457903511`), had independently
verified ZIP SHA-256
`875c013ad25b3486400c7e65ed359b6774b4cf42564c890543e663431a29afa3`.
The evidence-integrity, fragment-audit, fragment-candidate, de novo-audit, and
de novo-candidate components passed. The aggregate publication gate correctly
remained `technical_incomplete` because the competitive decoy, separate
current-database, and durable-deposition components are still open.

| Candidate | Independent fragment runs, validated / strong | De novo source recovery | De novo non-source recoveries | Six-library disposition |
|---|---:|---|---|---|
| A1 | 3 / 2 | DRR853912 | DRR853910, DRR853911 | Pass |
| A2 | 6 / 4 | DRR853912 | DRR853907, DRR853910, DRR853911 | Pass |
| B | 3 / 2 | DRR853910 | DRR853912 | Pass |

The separate current-sequence workflow run `32482881859` is still technically
incomplete. Local database/contamination checks, PF00680 domain analysis, and
UniVec screening passed; phylogenetic QC completed technically with no recorded
failure. Among the A1/A2/B candidate queries, no hit in the technically
complete remote modes met the workflow-defined near-identical criterion of
query coverage >=80% with protein identity >=90% or nucleotide identity
>=95%; the positive controls succeeded where included. However,
`protein_nonviral` attempts 1-6 produced structurally
invalid remote archives, attempt 7 ended with a transient remote-transport
failure, and the search budget was then exhausted. Remaining nucleotide modes
have not produced a complete validated gate. Those remote technical failures
are neither biological negative results nor evidence of novelty. No
current-database absence or novelty claim is allowed until the entire gate
completes with valid controls and provenance.

## Claims excluded

- three novel viruses or species
- viruses infecting *Panax notoginseng*
- healthy-versus-root-rot association or causation
- active replication, pathogenicity, or transmission
- complete genomes
- formal genus, family, or species assignment

## P0 submission blockers

1. Complete and freeze the current nr, nt, TSA, environmental, viral, and nonviral search gate with validated controls, search dates, parameters, and database provenance. The present `protein_nonviral` failure is a remote structural/technical failure and must not be reported as a no-hit result.
2. Complete competitive read assignment against *Panax* nuclear and organellar sequences, known *Panax* viruses, fungal/oomycete/root-community and other plausible eukaryotic decoys. The passed UniVec/local-contamination checks do not replace this broader decoy analysis.
3. Obtain a durable data-release plan: a permanent DOI archive for evidence/code and an appropriate INSDC/TSA/TPA accession route for public-read-derived sequences. GitHub Actions artifact `9457903511` is verified analysis evidence, but an expiring workflow artifact is not a Data Availability solution.
4. Freeze phenotype-neutral metadata and wording. Treat the `Zh/Bh` phenotype crosswalk as unresolved unless an authoritative mapping is obtained; otherwise omit every healthy/root-rot comparison. Retain the HiSeq 4000 versus NovaSeq 6000 platform discrepancy.
5. Finalize the current homolog/architecture context, including Pro-Pol, capsid/additional-ORF or putative-segment evidence where recoverable, and disclose the uncertainty of B's phylogenetic placement.

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
| A1 | Retain as a partial sequence cluster; publication conditional | Fragment gate passed in 3 independent runs (2 strong); mapping-seed-free de novo recovery passed in source DRR853912 and non-source DRR853910/DRR853911 | Close the study-wide current-database, decoy, architecture/context, and deposition blockers; explain cross-run sequence heterogeneity without assigning a species. |
| A2 | Retain as a partial sequence cluster; publication conditional | Fragment gate passed in all 6 runs (4 strong); de novo recovery passed in source DRR853912 and non-source DRR853907/DRR853910/DRR853911 | Close the study-wide blockers and keep A1/A2 as computationally distinct clusters rather than asserting two species. |
| B | Retain as a partial sequence cluster; publication conditional | Fragment gate passed in 3 independent runs (2 strong); de novo recovery passed in source DRR853910 and non-source DRR853912 | Close the study-wide blockers and broaden architecture/homolog sensitivity because phylogenetic placement remains uncertain. |

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
