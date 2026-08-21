# Six-library Panax candidate evidence audit

## Enhanced audit verdict

The enhanced six-library workflow completed successfully for all six complete
paired FASTQ run objects. Its evidence-integrity, duplicate-aware
proper-fragment, mapping-seed-free de novo, candidate-fragment, and candidate
de novo core components all passed technically.

All three predefined references are therefore **conditionally retained** as
partial Picornavirales-like RdRP-containing sequence clusters pending the
separate current-database gate. This is not evidence for three formal virus
species, three complete genomes, six independent plant infections, a true
*Panax* host, active replication, root-rot association, pathogenicity, or
transmission.

## Immutable artifact provenance

- Workflow run:
  [`32482881928`](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32482881928)
- Source commit: `605e5f18788eab15c9c2f9887c0e750f6058d552`
- Final artifact: `PANAX_SIX_LIBRARY_ENHANCED_AUDIT_20260821`
  (`9457903511`)
- ZIP size: `171,198,483` bytes
- ZIP SHA-256:
  `875c013ad25b3486400c7e65ed359b6774b4cf42564c890543e663431a29afa3`
- ZIP validation: valid archive; 388 top-level checksum entries passed and
  none failed. The nested fragment and de novo manifests passed 8/8 and 4/4,
  respectively.
- Evidence inventory: six validated run manifests, six reference manifests,
  18 masked consensuses, and all six full de novo assemblies retained and
  hash-linked.

GitHub Actions retention is temporary. The digest identifies the audited ZIP
but does not replace durable public deposition.

## Predeclared technical rules

These are workflow-defined sequence-support rules, not biological detection
limits, prevalence thresholds, or taxonomic criteria.

### Proper-fragment selection

- Both mates must be primary, QC-passing alignments to the same candidate,
  form a proper pair, and each have MAPQ >=30.
- Results are reported before duplicate marking and after `samtools markdup`
  exclusion.
- A validated run/reference requires at least 10 preduplicate proper fragments
  with preduplicate breadth at 1x >=0.80, plus at least five nonduplicate
  proper fragments with nonduplicate breadth at 1x >=0.50.
- A strong run/reference requires at least 25 preduplicate proper fragments,
  preduplicate breadth at 1x >=0.90 and at 5x >=0.50, plus at least 10
  nonduplicate proper fragments with nonduplicate breadth at 1x >=0.80.
- A candidate passes with at least two eligible independent validated runs and
  at least one eligible independent strong run.
- Exact sequence/endpoint fingerprints are compared across runs. A smaller run
  in a predeclared suspect or unevaluable shadow comparison cannot supply
  independent support. Shadow QC is not proof of index hopping or
  contamination.

### Mapping-seed-free de novo recovery

- Every assembly uses the complete paired FASTQ library without candidate
  baiting, mapping seeds, reference guidance, or target-read selection.
- Source-run recovery requires one collinear contig with query coverage >=0.90
  and coordinate-weighted nucleotide identity >=0.98.
- Non-source recovery requires one collinear contig with query coverage >=0.70
  and coordinate-weighted nucleotide identity >=0.90.
- Maximum permitted internal query and subject gaps are 150 nt and 500 nt.
- A candidate passes only with at least one source-run and one non-source-run
  recovery.

## Candidate-level gate result

| Sequence cluster | Eligible validated runs | Eligible strong runs | Fragment gate | Source de novo recovery | Non-source de novo recovery | De novo gate | Disposition |
|---|---|---|---|---|---|---|---|
| `PNX_Picorna_A1_ref` | `DRR853910`, `DRR853911`, `DRR853912` | `DRR853911`, `DRR853912` | pass | `DRR853912` | `DRR853910`, `DRR853911` | pass | conditionally retain pending current-DB gate |
| `PNX_Picorna_A2_ref` | all six runs | `DRR853907`, `DRR853910`, `DRR853911`, `DRR853912` | pass | `DRR853912` | `DRR853907`, `DRR853910`, `DRR853911` | pass | conditionally retain pending current-DB gate |
| `PNX_Picorna_B_ref` | `DRR853910`, `DRR853911`, `DRR853912` | `DRR853910`, `DRR853912` | pass | `DRR853910` | `DRR853912` | pass | conditionally retain pending current-DB gate |

Independent validated/strong counts are A1 3/2, A2 6/4, and B 3/2. The audit
contained 18 run/reference metric rows, 37,939 exact fragment fingerprints,
and 45 cross-run comparisons. Thirty-five comparisons had no predeclared
shadow signal; ten were unevaluable because the smaller run had fewer than ten
unique fingerprints. No comparison was warning or suspect. Unevaluable
run/reference states were excluded where required rather than treated as clean
evidence.

## De novo recovery detail

Coverage and identity below are the single-contig query coverage and
coordinate-weighted nucleotide identity used by the predeclared gate.

| Sequence cluster | Run | Role | Coverage | Identity | Gate |
|---|---|---|---:|---:|---|
| A1 | `DRR853912` | source | 1.000000 | 1.000000 | pass |
| A1 | `DRR853910` | non-source | 1.000000 | 0.919014 | pass |
| A1 | `DRR853911` | non-source | 0.786217 | 0.920665 | pass |
| A2 | `DRR853912` | source | 1.000000 | 1.000000 | pass |
| A2 | `DRR853907` | non-source | 1.000000 | 0.970016 | pass |
| A2 | `DRR853910` | non-source | 1.000000 | 0.938412 | pass |
| A2 | `DRR853911` | non-source | 1.000000 | 0.969611 | pass |
| B | `DRR853910` | source | 1.000000 | 0.999127 | pass |
| B | `DRR853912` | non-source | 1.000000 | 0.991851 | pass |

Per-run candidate recovery counts were 1/3, 0/3, 0/3, 3/3, 2/3, and 3/3 for
`DRR853907` through `DRR853912`, respectively. A failed candidate recovery in a
technically complete run is a sequence-level non-recovery under these
thresholds, not a workflow failure.

## Candidate structure and cross-run context

- A1 is a 1,988-nt single-source partial contig from `DRR853912`.
- A2 is a 2,468-nt single-source partial contig from `DRR853912`.
- B is a 3,436-nt reference derived from `DRR853910_21434` after removal of a
  validated 22-nt 5-prime terminal adapter prefix; its N-terminal ORF boundary
  remains open.
- The source-versus-non-source recoveries above demonstrate reproduction from
  complete independently archived run objects. Separate accessions do not
  establish independent plants or infections; the source article reports five
  plants pooled per replicate.
- Proper-nonduplicate masked-consensus identities are descriptive QC only.
  They are not used for transmission, strain phylogeny, or evolutionary-rate
  inference.

## Historical legacy audit

The earlier candidate-only mapping package remains useful as a historical
screen, but it is no longer the candidate decision source. Its ZIP SHA-256 was
`39af4a9a1351b22e67d13e08ebe3314883a37e97b58f5e3858f3d4fe29e89abf`
from run
[`32454265572`](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32454265572).
That package did not mark duplicates, require both-mate MAPQ30 proper
fragments, inspect exact fingerprints, or perform complete-read de novo
assembly. The enhanced artifact supersedes its legacy retain/provisional
labels.

## Remaining alternative explanations and blockers

- The current competitive mapping panel contains only A1, A2, and B. It lacks
  a version-pinned *Panax notoginseng* nuclear/organelle panel, rRNA, known
  *Panax*-associated viruses, plausible fungal/oomycete/root-community loci,
  PhiX, and a complete adapter/vector panel. MAPQ from this three-target panel
  therefore does not establish broader biological specificity.
- No warning/suspect exact-fingerprint shadow pattern was found, but exact
  fingerprints cannot prove or exclude index hopping, especially for
  comparisons with insufficient smaller-run evidence.
- `Zh` and `Bh` remain archive tokens. No condition-level statistical
  comparison is performed because their mapping to the article's `CK`/`ROT`
  labels is unresolved.
- Candidate sequences are partial. Architecture, additional segments, true
  host, replication, and formal taxonomic status remain unresolved.
- The current-database sequence gate is a separate fail-closed workflow and
  must complete before a manuscript-level retention/novelty statement is
  frozen.
- Durable public sequence and evidence deposition has not been completed.

## Publication gate

The six-library evidence package is technically complete and passes its core
candidate-support gates, but the overall publication package remains
`technical_incomplete` and `submission_ready=false` because the competitive
decoy audit, current-database sequence gate, and durable public deposition are
not complete.

The six-library result supports this bounded wording:

> Three predefined divergent partial RdRP-containing sequence clusters with
> Picornavirales-like homology were recovered and raw-read-supported in public
> *Panax notoginseng*-associated root RNA-seq libraries.

Use in a submission remains conditional on closing the current-database,
competitive-decoy, and durable-deposition blockers and must not be expanded to
a formal species, true-host, infection, replication, root-rot, causality,
pathogenicity, transmission, or agricultural/medical claim.
