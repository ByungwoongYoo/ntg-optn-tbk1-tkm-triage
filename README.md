# NTG OPTN–TBK1 TKM in silico triage

## Separate formal-logic result: Ulrich u4

The `u4_resolution/` directory contains an explicit infinite syntactic
countermodel showing that Ulrich's `u4` formula is not a single axiom for
positive implicational logic. The formal proof, exact problem correspondence,
and independent-verification design are documented in
[`u4_resolution/README.md`](u4_resolution/README.md) and
[`u4_resolution/U4_FORMAL_SCOPE_20260820.md`](u4_resolution/U4_FORMAL_SCOPE_20260820.md).

This result is logically and methodologically separate from the biomedical
docking study described below. Historical verification records remain in the
repository; strengthened releases must be made only from an all-green commit of
the dedicated Lean 4.33 workflow.

Corrected v2 reserved DOI:
[10.5281/zenodo.22031656](https://doi.org/10.5281/zenodo.22031656).
All versions belong under [concept DOI
10.5281/zenodo.21987698](https://doi.org/10.5281/zenodo.21987698); the historical
v1 is [10.5281/zenodo.21987699](https://doi.org/10.5281/zenodo.21987699).
The invalid identifier `10.5281/zenodo.21987529` must not be cited.

This repository contains the **processed data, analysis outputs, and figure source files** for a
hypothesis-generating *in silico* triage study of traditional Korean medicine (TKM)-derived
compounds at the OPTN–TBK1 interface region, in the context of normal-tension glaucoma (NTG).

> **Scope / disclaimer.** This work does **not** claim therapeutic efficacy, protein–protein
> interaction (PPI) disruption, or E50K-selective inhibition. It prioritizes *interface-adjacent
> candidate compounds* and reports methodological findings, including a blind-center grid-placement
> artifact, a glycoside size-bias case, and the limits of rigid docking for mutation-state
> selectivity. All results are computational prioritizations and require experimental validation.

## What is in this repository

```
data/
  final_tables/          Final decision table, paper Table 2 (processed), decision README
  docking_scores/        Per-pose scores with ligand efficiency
  pose_contact_metrics/  Pose/contact summary, RMSD summary, raw pose-contact table, verdicts
figures/
  Figure1..5 (PNG previews) + source_svg/ (editable vector sources)
docs/
  software_versions.txt  Exact tool versions used
  methods_notes.md       Short methods notes (grids, thresholds, definitions)
```

This repository intentionally provides **processed data and analysis files**, not every raw
docking output. Source protein structures are the experimental complexes **PDB 5EOF** (wild-type
OPTN(26–103)/TBK1 CTD) and **PDB 5EOA** (E50K), publicly available from the RCSB PDB.

## Key methodological points

- Only the residue-50 side-chain-averaged grid (identical wild-type/E50K centers) constitutes a
  controlled mutation-state comparison; independent blind-center grids were up to ~11.9 Å apart.
- Pose reproducibility = fraction of receptor systems with a low-RMSD pose cluster
  (median pairwise heavy-atom RMSD ≤ 2.0 Å).
- Ligand efficiency = docking score divided by heavy-atom count.

## Software

Python 3.12.13, Uni-Dock v1.1.3, Open Babel 3.1.0, UCSF ChimeraX 1.12. Exact versions are in
`docs/software_versions.txt`.

## Citation

A Zenodo DOI will be added here after repository archiving. Associated manuscript: Yoo B.
"Hypothesis-generating computational triage of traditional Korean medicine-derived compounds as
candidate OPTN–TBK1 interface-region compounds in normal-tension glaucoma" (submitted).

## License

Data and documentation are released under CC BY 4.0 (see `LICENSE`).
