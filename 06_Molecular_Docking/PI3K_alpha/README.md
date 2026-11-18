# PI3Kα Molecular Docking — Tryptanthrin Analogues

**Responsible member:** Ayana  
**Target:** Human PI3Kα / PIK3CA  
**Crystal structure:** 4JPS (catalytic chain A and regulatory chain B retained in initial preparation)  
**Reference ligand:** 1LT (alpelisib)  
**Software:** RDKit, PDBFixer, Meeko, AutoDock Vina 1.2.7  

## Aim
Dock the 10 PI3Kα-specific candidates selected in the ML workstream, examine poses and residue contacts, and rank computational hypotheses with transparent validation limitations.

## Key validation and settings

- Atom-name-mapped crystallographic 1LT SDF prepared using the PDB chemical-component definition and deposited 4JPS coordinates.
- 1LT redocking: **top-ranked symmetry-aware heavy-atom RMSD 0.6039 Å** relative to the crystal coordinates; no pose superposition for scoring.
- AutoDock Vina: **1.2.7**, exhaustiveness **32**, modes **9**, seed **42**, default scoring function.
- Grid centre (Å): **(-1.319, -9.513, 16.948)**; box size (Å): **(18.946, 19.659, 25.114)**.
- Receptor preprocessing reconstructed missing atoms; Meeko required **distance-restricted exclusion of ten problematic residues** more than 20 Å outside the provisional docking box. Earlier checks showed all **85 original pocket heavy atoms within 5 Å** present in prepared receptor, with **0.0 Å displacement**. This does not validate the entire receptor chemistry/protonation.
- These experimental-pose-recovery results support the chosen redocking protocol; they do not demonstrate experimental inhibition or accurate affinity estimation.

## Candidate docking

All **10** analogues in `ligands/PI3Kalpha_selected_top10.csv` were prepared and docked. The three most favourable Vina scores were:

| Rank | Analog_ID | Vina score (kcal/mol) |
|---:|---|---:|
| 1 | TRY_6db5e7a7802259ab | -10.290 |
| 2 | TRY_c88335488290afb4 | -10.261 |
| 3 | TRY_c6b27d68ce6ec718 | -10.175 |

The top-three pose-contact assessment recorded **15, 13, and 14 residues within 4 Å**, respectively. Proximity alone does not establish hydrogen bonding or inhibitory activity. ML applicability-domain similarity was low (~0.40–0.41) for these candidates; see the complete results report and ranking for flags and caveats. Small score differences should not be interpreted as meaningful without further validation.

## Folder guide

- `structures/`: original 4JPS, repaired heavy-atom receptor, prepared receptor PDBQT, atom-mapped crystal/reference ligand.
- `ligands/`: 10 prepared SDF and PDBQT files, original IDs, preparation status.
- `validation/`: reference redocking poses, affinities, RMSD and redocking report.
- `results/`: all 10 docked ligand PDBQTs, final ranking CSV, score figures, scientific report.
- `analysis/`: three top-pose PDBs and HTML visualizations, atom/residue contact tables and quality metrics.
- `logs/`: conversion and Vina logs for auditability.
- `notebooks/`: executed notebooks from redocking through integrated review; `development/` preserves the ligand mapping and receptor repair working records. Failed preparation attempts in the working records are historical; the final prepared receptor and validation files above are the results used for docking.

Start with `results/PI3Kalpha_Final_Results_Report.md` and `results/PI3Kalpha_final_candidate_ranking.csv`.

## Reproducibility and limitations

This is an **evidence and code handover**, not a verified, single-click end-to-end pipeline. The notebooks reflect the actual staged Colab workflow and may request earlier ZIPs or require paths to be adapted to this reorganized folder. The completed end-to-end pipeline has **not** been rerun from this GitHub folder layout. Retain raw preparation records in the development notebooks; do not represent excluded residues as repaired.

No experimental activity, binding free energy, selectivity or clinical benefit is established. Receptor mutations, chain/cofactor choices, docking-box definition, ligand protonation, low ML similarity and excluded residues should be considered when interpreting results. Further orthogonal modelling and experimental testing would be needed to validate these candidates.
