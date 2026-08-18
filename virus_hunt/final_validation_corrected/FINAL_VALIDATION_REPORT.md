# Final Panax-associated RdRP validation report

Generated (UTC): 2026-08-18T14:37:38.408662+00:00

## Overall result

**`NO_CANDIDATE_MET_THE_FULL_SEQUENCE_LEVEL_THRESHOLD`**

- Strict divergent RdRP lineage candidates: **0**
- Supported but archived-sample independence limited: **0**
- Current near-identity detections: **0**
- Nonviral/ambiguous preferred matches: **0**
- ENA independence grade: **`distinct_archived_samples_and_experiments`**
- Required current-database searches complete: **False**

The strongest permitted interpretation is a sequence-level divergent RNA-dependent RNA polymerase (RdRP) lineage candidate in public Panax-associated transcriptomic data. This report does **not** establish a new virus species, active viral replication, infection of Panax, pathogenicity, transmissibility, or any medical effect.

## Per-lineage evidence

| Lineage | Decision | Assembly runs | Raw-supported runs | PALMdb identity range | Current viral hit | Protein id/qcov | NT id/qcov | Tree |
|---|---|---:|---:|---:|---|---:|---:|---|
| `PNX_Duplo_A` | `database_audit_incomplete` | 3 | 3 | 67.7-67.7% | N/A | 99.725/100 | 100.000/100 | true |
| `PNX_Duplo_B` | `database_audit_incomplete` | 3 | 3 | 71.3-71.3% | N/A | 100.000/100 | 100.000/98 | true |
| `PNX_Picorna_A` | `database_audit_incomplete` | 2 | 3 | 60.2-60.5% | N/A | 69.438/100 | NA/NA | true |
| `PNX_Picorna_B` | `database_audit_incomplete` | 2 | 2 | 55.0-55.0% | N/A | 41.350/75 | NA/NA | true |

## Strict evidence gates

A strict candidate required: complete A/B/C palm motifs; concordant PSSM/HMM and PALMdb support; strong cross-run assembly recurrence; raw-read coverage of at least 90% at mean depth at least 5× with at least 20 MAPQ≥20 reads in at least two runs; successful current viral, nonviral, unrestricted-protein, and nucleotide database searches; no protein hit at ≥90% identity over ≥80% of the query; no nucleotide hit at ≥95% identity over ≥80% of the query; and no substantially stronger nonviral explanation.

## Independence boundary

The ENA audit grade is `distinct_archived_samples_and_experiments`. Distinct archive accessions do not by themselves prove independent plants, independent infections, or the true host.

## What remains necessary for a formal virus discovery

Formal taxonomic or host claims would require expert viral-taxonomy review, complete or substantially extended genome architecture, stronger contamination/index-hopping exclusion, and ideally independent biological sampling with strand-aware or targeted validation.
