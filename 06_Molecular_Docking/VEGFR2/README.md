# VEGFR2 docking report (PDB 3VHE) - Abhi

## Methods
- Receptor: PDB 3VHE (Crystal structure of human VEGFR2 kinase domain with a novel pyrrolopyrimidine inhibitor.), resolution 1.55 A, chain A; crystal ligand 42Q.
  First alternate conformer kept; waters, ligands and unrelated hetero atoms removed; hydrogens added at pH 7.4 (PDBFixer); missing loops not rebuilt; PDBQT by Meeko.
- Docking box (from crystal ligand centroid): centre [-24.968, -1.141, -10.515], size [22.0, 22.0, 22.0] A.
- Ligands: 10 ML-shortlisted analogs for VEGFR2; RDKit ETKDGv3 + MMFF optimisation; Meeko PDBQT.
- AutoDock Vina 1.2.7 (Vina scoring), exhaustiveness 32, 9 modes, seed 42.

## Redocking validation
{'ligand': '42Q', 'rmsd_model1': 0.55, 'rmsd_best_of_modes': 0.55, 'best_mode': 1, 'vina_model1': -12.722, 'cutoff': 2.0, 'passed': True}

## Docking and ranking (lower composite score = better)
| Composite_rank | Analog_ID | Vina_best_kcal_mol | ML_pIC50 | QED | SA_score | PAINS_flag | Brenk_flag |
|---|---|---|---|---|---|---|---|
| 1 | TRY_7b1bcb6a14c88658 | -11.29 | 6.56 | 0.58 | 2.26 | False | False |
| 2 | TRY_c8cdc786b96d72b0 | -12.03 | 6.26 | 0.56 | 2.46 | False | False |
| 3 | TRY_54049cabf57d8af1 | -9.62 | 6.39 | 0.56 | 2.56 | False | False |
| 4 | TRY_c4f6d2b11da089fb | -9.11 | 6.18 | 0.6 | 2.34 | False | False |
| 5 | TRY_50a05a127072ea56 | -9.66 | 6.31 | 0.56 | 2.53 | False | False |
| 6 | TRY_83e3aaeb142feba2 | -10.63 | 6.33 | 0.54 | 2.54 | False | False |
| 7 | TRY_2327d14018ee1ca2 | -9.22 | 6.34 | 0.52 | 2.37 | False | False |
| 8 | TRY_aff7ade78c2db7c6 | -9.0 | 6.17 | 0.58 | 2.4 | False | False |
| 9 | TRY_aec2ddccc96a3bad | -8.7 | 6.39 | 0.52 | 2.64 | False | False |
| 10 | TRY_4f2df5013a804bb2 | -8.29 | 6.31 | 0.52 | 2.63 | False | False |

Selected for further work: TRY_7b1bcb6a14c88658, TRY_c8cdc786b96d72b0, TRY_54049cabf57d8af1

## Interaction analysis (top 3, Model 1)
H-bond = donor...acceptor <= 3.5 A and D-H...A angle >= 120 deg. Hydrophobic/aromatic entries are distance-based flags only.
           Analog_ID n_hbond     hbond_residues n_hydrophobic_res n_contact_res
TRY_7b1bcb6a14c88658       2 CYS919:A, LEU840:A                 7            15
TRY_c8cdc786b96d72b0       1           CYS919:A                 7            19
TRY_54049cabf57d8af1       1          ILE1044:A                 3            17

## Limitations
Vina scores are approximate and are not comparable between different targets. Ligand protonation/tautomer states were not enumerated.
All results are computational predictions and are not experimentally confirmed activity.
