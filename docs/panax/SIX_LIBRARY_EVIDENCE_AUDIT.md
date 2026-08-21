# Six-library Panax candidate evidence audit

## Audit verdict

The existing six-library artifact correctly calculates read support for three competitive reference sequences. All archived metrics, depth summaries, masked consensuses, presence calls, pairwise identities, and checksum manifests were independently recalculated without discrepancy.

This verifies only that `A1-like`, `A2-like`, and `B-like` partial sequences receive quality-filtered read support in public *Panax notoginseng*-associated root RNA-seq. It does not verify three virus species, three complete genomes, six independent plant infections, a true *Panax* host, active replication, root-rot association, pathogenicity, or transmission.

Artifact ZIP SHA-256:

`39af4a9a1351b22e67d13e08ebe3314883a37e97b58f5e3858f3d4fe29e89abf`

Artifact locator: GitHub Actions run
[`32454265572`](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32454265572),
artifact `PANAX_SIX_LIBRARY_CLEAN_AUDIT_20260819`. GitHub Actions retention is temporary; the checksum
identifies the audited ZIP but does not replace permanent deposition.

Source commit:

`f6678e538237627e478438bcbac61156d1f27a5b`

## Fixed workflow definitions

- Quality-filtered depth: mapping quality >=20, base quality >=20, overlapping mates counted once.
- Strict support: at least 20 MAPQ>=20 read records and breadth at 1x >=0.80.
- Strong support: at least 50 MAPQ>=20 read records, breadth at 1x >=0.90, and breadth at 5x >=0.50.
- Masked consensus: haploid SNP-only consensus with positions below quality-filtered 5x depth changed to `N`.
- `mapped_mapq20` is an alignment/read-record count, not a read-pair count.

These are fixed technical rules, not validated clinical, biological, prevalence, or taxonomic thresholds.

## Recalculated support summary

| Candidate | Zh-token strict / strong | Bh-token strict / strong | Evidence interpretation |
|---|---:|---:|---|
| A1 | 0/3 / 0/3 | 3/3 / 3/3 | Strong A1-like support in all three Bh-token runs; lower support in all Zh-token runs remains below the fixed strict threshold. |
| A2 | 3/3 / 1/3 | 3/3 / 3/3 | Strict support in all six runs, but ZhB and ZhC are weak positives and low-level carryover is not excluded. |
| B | 0/3 / 0/3 | 3/3 / 2/3 | Strong in BhA and BhC; BhB is weak and threshold-sensitive. |

`Zh` and `Bh` are archive tokens only. No condition-level statistical comparison is performed.

## Run-level evidence

| Candidate | Run (token) | MAPQ>=20 records | Breadth 1x / 5x | Mean depth | Strict / strong |
|---|---|---:|---:|---:|---|
| A1 | DRR853907 (ZhA) | 54 | 67.7% / 38.5% | 3.67x | false / false |
| A1 | DRR853908 (ZhB) | 111 | 74.1% / 46.6% | 6.98x | false / false |
| A1 | DRR853909 (ZhC) | 35 | 73.8% / 10.4% | 2.34x | false / false |
| A1 | DRR853910 (BhA) | 384 | 93.6% / 82.1% | 24.62x | true / true |
| A1 | DRR853911 (BhB) | 103,714 | 97.1% / 92.9% | 6,794.07x | true / true |
| A1 | DRR853912 (BhC) | 490 | 100% / 100% | 31.27x | true / true |
| A2 | DRR853907 (ZhA) | 236 | 100% / 81.4% | 11.30x | true / true |
| A2 | DRR853908 (ZhB) | 43 | 86.8% / 8.8% | 2.16x | true / false |
| A2 | DRR853909 (ZhC) | 57 | 94.2% / 21.2% | 2.97x | true / false |
| A2 | DRR853910 (BhA) | 30,997 | 100% / 100% | 1,678.53x | true / true |
| A2 | DRR853911 (BhB) | 26,895 | 100% / 100% | 1,460.44x | true / true |
| A2 | DRR853912 (BhC) | 7,138 | 100% / 100% | 376.48x | true / true |
| B | DRR853907 (ZhA) | 50 | 74.1% / 6.9% | 1.80x | false / false |
| B | DRR853908 (ZhB) | 16 | 19.0% / 3.1% | 0.43x | false / false |
| B | DRR853909 (ZhC) | 53 | 73.9% / 10.3% | 1.98x | false / false |
| B | DRR853910 (BhA) | 454 | 100% / 99.7% | 17.06x | true / true |
| B | DRR853911 (BhB) | 48 | 81.2% / 9.5% | 1.83x | true / false |
| B | DRR853912 (BhC) | 401 | 100% / 96.4% | 15.26x | true / true |

## High-coverage consensus checks

- B has the cleanest non-source **masked-consensus mapping** reproduction: BhA and BhC overlap at 3,309 of 3,436 positions after 5x masking and are 99.365% identical. This is not de novo reassembly evidence.
- A2 has broad support, but high-coverage consensuses include sequence populations separated by approximately 2-6% across runs.
- A1 has strong multi-run A1-like support, but BhC versus BhA/BhB masked consensuses are approximately 93% identical. The evidence therefore supports an A1-like sequence family rather than the identical selected A1 sequence in every run.
- A1 and A2 RdRP segments are 82.353% identical across 408 aligned amino acids and share no exact 31-mers in the selected nucleotide references. This supports computational separability, not separate virus species.

Masked-consensus identities remain descriptive. The workflow does not apply an additional SNP QUAL or within-sample allele-fraction filter after the 5x mask, so they are not used for strain phylogeny, transmission inference, or evolutionary-rate analysis.

## Alternative explanations not yet excluded

- Low-level A2 support in Zh-token runs is compatible with genuine low abundance, index hopping, or cross-sample contamination; lane and index metadata are unavailable.
- A1 BhA has a low proper-pair fraction in the existing candidate-only competitive mapping and requires proper-pair-only sensitivity analysis.
- Duplicate flags are zero because no duplicate-marking stage was run. The very high A1 depth in BhB cannot be interpreted as unique-molecule abundance.
- The competitive reference contains only A1, A2, and B. It is not a comprehensive *Panax* nuclear/organelle, fungal, oomycete, soil, vector, and known-virus decoy.
- Source-run mapping is circular support for an assembled source contig. Non-source reconstruction and contig-spanning evidence are required for stronger corroboration.
- The compact artifact omits raw FASTQ and BAM files. It records their hashes but cannot by itself reproduce mapping from the reads or inspect read names and fragment endpoints.

## Legacy mapping-only candidate-level decision

The dispositions below describe the audited legacy mapping artifact only. They
do not pre-empt the pending enhanced proper-fragment, mapping-seed-free de novo,
and current-database gates.

| Candidate | Retain? | Defensible current interpretation | Required upgrade before manuscript submission |
|---|---|---|---|
| A1 | Legacy mapping-only retain | Divergent A1-like partial Picornavirales/RdRP sequence family with strong multi-run read support | Proper-pair and duplicate-aware remapping; non-source reassembly; decoy competition; explain approximately 93% cross-run consensus identity. |
| A2 | Legacy mapping-only provisional | Broadest read support and distinct from A1 at the selected sequence level | Obtain an independent non-source assembly; otherwise downgrade or exclude it from the primary discovery count; audit low-level carryover. |
| B | Legacy mapping-only retain | Partial Picornavirales-like sequence with near-full non-source BhC masked-consensus mapping reproduction; not de novo reassembly | Duplicate-aware/proper-pair and decoy remapping; broaden homolog/architecture analysis because tree placement is unstable. |

No candidate is called a formal virus species or a complete genome.

## Publication gate

The existing package is suitable for internal sharing with the above caveats, but it is not yet a submission-grade discovery package. The publication gate remains open until the added read-level controls, complete current-database sequence gate, source-metadata correction, permanent evidence deposition plan, and candidate architecture analyses are complete.
