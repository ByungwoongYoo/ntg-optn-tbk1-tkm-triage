# Methods notes (concise)

**Structures.** PDB 5EOF (wild-type OPTN(26–103)/TBK1 CTD) and 5EOA (E50K). Waters/heteroatoms
removed; hydrogens added; protonation at pH 7.4 (Open Babel 3.1.0). Two chain pairs (AC, BD)
treated as independent receptor systems; E50K and BD superposed onto the wild-type AC frame.

**Ligands.** Built from PubChem canonical SMILES with RDKit (ETKDG + MMFF94). Protoberberines
(coptisine, berberine) prepared as +1 cations. Controls: quercetin (PAINS-suspect flavonoid),
amlexanox (TBK1 ATP-site reference), acacetin (linarin aglycone), berberine (protoberberine ref).

**Docking.** Uni-Dock v1.1.3, Vina scoring, exhaustiveness 32. Two grids: (1) blind-center
(28 Å box) from independent cavity detection per structure — centers NOT identical between states
(AC 11.91 Å apart, BD 1.12 Å); (2) residue-50 side-chain-averaged (26 Å box) — identical centers
(0.00 Å), the only controlled WT/E50K comparison. Core panel 6 compounds × 10 seeds; focused panel
5 compounds × 3 seeds, across four receptor systems.

**Metrics.** Contact cutoff 4.0 Å. Dual-chain contact = ≥1 OPTN + ≥1 TBK1 residue. Strong cleft =
≥2 OPTN + ≥2 TBK1 and ≥5 total unique contact residues. Residue-50 proximity flag ≤ 8 Å.
Ligand efficiency = score / heavy-atom count. Pose reproducibility = fraction of systems with
median pairwise heavy-atom RMSD ≤ 2.0 Å.

**Not done.** MM-GBSA / MM-PBSA free-energy refinement was deliberately omitted (interface-adjacent,
limitedly reproducible poses → false precision). No experimental validation.
