# Stage 16 Pose/Contact Analysis

Purpose: post-docking qualitative/structural triage for OPTN-TBK1 interface-region candidate binders.

Inputs used:
- Stage14A Uni-Dock core6 blind seed10 results
- Stage14A Uni-Dock core6 sidechain seed10 score summaries (pose files unavailable in uploaded snapshot)
- Stage15B focused selected-7 blind seed3 results
- Stage15B focused selected-7 sidechain seed3 results
- receptor PDBs from Stage14A/Stage15B packages

Important limitation:
- This is an automated geometric contact screen using the top model in each PDBQT file.
- It is not a wet-lab validation and not MM-GBSA.
- Some control compounds have pose analysis from available blind-grid outputs only.
- Manual PyMOL/ChimeraX visualization remains required before MM-GBSA.

Definitions:
- Dual-chain contact: ligand has at least one contact residue on OPTN chain and one on TBK1 chain within 4 Å.
- Strong cleft contact: at least 2 OPTN residues + at least 2 TBK1 residues + at least 5 total unique receptor contact residues.
- Res50 proximity: minimum heavy-atom distance to OPTN residue 50 <= 8 Å.
- Pose reproducibility is an approximate same-atom-order RMSD across seed top poses.

Current interpretation:
- Coptisine remains the main candidate for MM-GBSA consideration if visual cleft burial is confirmed.
- Cryptotanshinone and Tanshinone IIA are new scaffold candidates but need scaffold-bias/flat aromatic stacking visual review.
- Linarin is a method/glycoside size-bias caution case rather than a main candidate.
