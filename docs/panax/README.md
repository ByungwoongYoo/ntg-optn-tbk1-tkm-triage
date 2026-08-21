# Panax public-sequence candidate audit

This directory separates the Panax public-transcriptome analysis from the
repository's original OPTN-TBK1 project. It documents a computational sequence
audit, not a claim of virus discovery, host infection, or root-rot association.

## Current decision

The analysis should continue, but it is not submission-ready. The legacy
mapping artifact contains three predefined partial Picornavirales-like
RdRP-containing **sequence clusters** (`A1-like`, `A2-like`, and `B-like`). The
retained count remains pending until the enhanced six-run fragment, de novo,
and current-database gates finish. They are not three established viruses or
species.

- `A1-like`: legacy-retain; enhanced disposition pending.
- `A2-like`: provisional and enhanced disposition pending; a non-source de novo assembly is required;
  otherwise downgrade or exclude it from the primary candidate count.
- `B-like`: legacy strongest mapping candidate; enhanced disposition pending, while its broad
  phylogenetic context remains unstable.

## Audit map

- [`SOURCE_CONTEXT_AUDIT.md`](SOURCE_CONTEXT_AUDIT.md): source study, archive
  crosswalk, pooling, platform discrepancy, and corrected phenotype boundary.
- [`RUN_CROSSWALK.tsv`](RUN_CROSSWALK.tsv): exact archive identifiers.
- [`SIX_LIBRARY_EVIDENCE_AUDIT.md`](SIX_LIBRARY_EVIDENCE_AUDIT.md): independent
  recalculation of the original six-run mapping artifact and its limitations.
- [`PUBLICATION_READINESS.md`](PUBLICATION_READINESS.md): candidate disposition,
  mandatory blockers, manuscript scope, and journal strategy.

## Conditional working claim

> If all three enhanced candidate gates pass: three divergent partial
> RdRP-containing sequence clusters with
> Picornavirales-like homology were recovered and raw-read-supported in public
> *Panax notoginseng*-associated root RNA-seq libraries.

Until then, the safe statement is that three predefined clusters received
legacy candidate-panel mapping support and the final retained count is pending.
Neither wording establishes the true host, active replication,
pathogenicity, transmission, condition association, genome completeness, or a
formal taxon. `Zh` and `Bh` remain archive tokens because no authoritative
mapping to the source article's `CK` and `ROT` condition labels was recovered.
