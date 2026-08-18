# Public-sequencing virus-discovery audit (isolated branch)

This directory begins an auditable search for RNA-virus-like sequences hidden in public sequencing data. It is intentionally isolated on branch `virus-discovery-20260818`; the repository's original NTG project and `main` branch are not modified.

## Scientific boundary

A metadata hit, an RdRP-like sequence, or even a complete virus-like genome is **not** evidence that the nominal host was infected or that the sequence causes disease. Host attribution, replication, pathogenicity, and biological effect require additional evidence. Computational outputs are graded conservatively.

## Phase 1

`scan_recent_runs.py` searches ENA for RNA-seq runs from medicinal plants and medicinal fungi first released on or after 2021-01-01. This date deliberately starts after the original Serratus 2020-scale sweep. Projects are ranked only for tractability and suitability for downstream screening; no viral classification is made.

Planned evidence gates after metadata triage:

1. sequence-level RdRP evidence with conserved motifs;
2. current database similarity and novelty checks;
3. raw-read reassembly and read-back mapping;
4. independent-run or independent-project replication where available;
5. contamination and index-hopping checks;
6. phylogenetic placement with explicit uncertainty;
7. SHA-256 manifest and exact provenance for every claim.
