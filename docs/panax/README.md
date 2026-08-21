# Panax public-sequence candidate audit

This directory separates the Panax public-transcriptome analysis from the
repository's original OPTN-TBK1 project. It documents a computational sequence
audit, not a claim of virus discovery, host infection, or root-rot association.

## Current decision

The enhanced six-library audit is complete, and its core evidence components
passed. All three
predefined partial Picornavirales-like RdRP-containing **sequence clusters**
(`A1-like`, `A2-like`, and `B-like`) passed both the independent proper-fragment
gate and the mapping-seed-free de novo source/non-source recovery gate. They
are conditionally retained pending the separate current-database sequence
gate. They are not three established viruses or species.

- `A1-like`: fragment support in three independent eligible runs, including
  two strong runs; source recovery in `DRR853912` and non-source recovery in
  `DRR853910` and `DRR853911`.
- `A2-like`: fragment support in all six runs, including four strong runs;
  source recovery in `DRR853912` and non-source recovery in `DRR853907`,
  `DRR853910`, and `DRR853911`.
- `B-like`: fragment support in three independent eligible runs, including two
  strong runs; source recovery in `DRR853910` and non-source recovery in
  `DRR853912`.

The package remains **not submission-ready** because the current-database gate
has not completed, the competitive host/organelle/microbial/laboratory decoy
panel is incomplete, and durable public sequence/evidence deposition has not
been completed.

## Audit map

- [`SOURCE_CONTEXT_AUDIT.md`](SOURCE_CONTEXT_AUDIT.md): source study, archive
  crosswalk, pooling, platform discrepancy, and corrected phenotype boundary.
- [`RUN_CROSSWALK.tsv`](RUN_CROSSWALK.tsv): exact archive identifiers.
- [`SIX_LIBRARY_EVIDENCE_AUDIT.md`](SIX_LIBRARY_EVIDENCE_AUDIT.md): independent
  enhanced proper-fragment and mapping-seed-free de novo audit, immutable
  artifact provenance, candidate gates, and remaining limitations.
- [`PUBLICATION_READINESS.md`](PUBLICATION_READINESS.md): candidate disposition,
  mandatory blockers, manuscript scope, and journal strategy.

## Current working claim boundary

The six-library artifact supports the following wording at the read-support
and de novo-recovery level:

> Three predefined divergent partial RdRP-containing sequence clusters with
> Picornavirales-like homology were recovered and raw-read-supported in public
> *Panax notoginseng*-associated root RNA-seq libraries.

Use of that sentence in a submission remains conditional on closing the
current-database, competitive-decoy, and durable-deposition blockers and must
be accompanied by the unresolved-host and phenotype limitations. It does not
establish the true host, active
replication, pathogenicity, transmission, condition association, genome
completeness, or a formal taxon. `Zh` and `Bh` remain archive tokens because no
authoritative mapping to the source article's `CK` and `ROT` condition labels
was recovered.

## Enhanced artifact provenance

- Workflow run:
  [`32482881928`](https://github.com/ByungwoongYoo/ntg-optn-tbk1-tkm-triage/actions/runs/32482881928)
- Artifact: `PANAX_SIX_LIBRARY_ENHANCED_AUDIT_20260821` (`9457903511`)
- Source commit: `605e5f18788eab15c9c2f9887c0e750f6058d552`
- ZIP size: `171,198,483` bytes
- ZIP SHA-256:
  `875c013ad25b3486400c7e65ed359b6774b4cf42564c890543e663431a29afa3`
- Internal checksum validation: 388 passed, 0 failed; all six full assemblies
  were retained and hash-linked.
