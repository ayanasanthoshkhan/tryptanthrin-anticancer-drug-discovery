
# TDO2 Molecular Docking Analysis

## 1. Objective
To evaluate ten selected tryptanthrin-derived compounds against
Tryptophan 2,3-dioxygenase (TDO2) using molecular docking,
machine-learning predictions and drug-likeness properties.

## 2. Protein Structure Preparation
- Target: TDO2
- PDB ID: 8VTQ
- Chain: D
- Reference ligand: A1ADX
- Heme cofactor and iron retained.
- Other crystallographic ligands and waters removed.

## 3. Receptor Validation
Reference ligand A1ADX was redocked using AutoDock Vina.

- Exhaustiveness: 32
- Number of poses: 9
- Random seed: 42
- Best docking score: -9.465 kcal/mol
- Heavy-atom RMSD: 1.488 Angstrom
- Acceptance threshold: 2.0 Angstrom
- Validation: PASSED

## 4. Ligand Preparation
Ten compounds were selected from the TDO2 shortlist.

- SMILES validated using RDKit.
- 3D structures generated using ETKDG.
- Geometry optimized using MMFF or UFF.
- PDBQT preparation performed using Meeko.
- All ten compounds successfully prepared.

## 5. Molecular Docking

AutoDock Vina settings:
- Exhaustiveness: 32
- Number of poses: 9
- Random seed: 42

### Docking Results

| Analog ID | Docking score (kcal/mol) |
|---|---:|
| TRY_91b948175b25d8f6 | -9.837 |
| TRY_ebbb95c0d5acdd8d | -9.592 |
| TRY_3ba689cad029a158 | -9.437 |
| TRY_c557a68e91178c2d | -9.305 |
| TRY_de48e599145298e2 | -9.289 |
| TRY_b8ef54ecfea0d3cf | -9.243 |
| TRY_6ebce10324125fed | -9.070 |
| TRY_a38fbe6b04a5d8cf | -8.771 |
| TRY_44c5e8a1eb8f55fe | -8.758 |
| TRY_b38affb4ddf93e6b | -8.542 |

## 6. ML and Drug-Likeness Analysis

### Top Three Candidates

| Analog ID | Consensus pIC50 | Docking score | QED |
|---|---:|---:|---:|
| TRY_91b948175b25d8f6 | 6.54089 | -9.837 | 0.491391 |
| TRY_ebbb95c0d5acdd8d | 6.53251 | -9.592 | 0.495720 |
| TRY_3ba689cad029a158 | 6.50649 | -9.437 | 0.540405 |

All three candidates had negative PAINS and Brenk flags.

## 7. Protein-Ligand Interaction Analysis

Important nearby residues:
- HIS76
- PHE72
- LEU147
- ARG144
- ALA150
- SER151
- GLY152

The three compounds showed potential polar contacts involving
HIS76 and GLY152, hydrophobic contacts and proximity to HEM.

These are preliminary distance-based observations and do not
confirm hydrogen bonding or direct iron coordination.

## 8. Discussion

TRY_91b948175b25d8f6 obtained the most favourable predicted
docking score among the ten screened compounds.

TRY_ebbb95c0d5acdd8d ranked second, while
TRY_3ba689cad029a158 ranked third and showed the highest
QED among the shortlisted candidates.

The three compounds showed similar binding-site contact patterns.

## 9. Limitations

- Docking scores are predictions, not experimental binding affinities.
- Vina has limitations in modelling heme iron coordination.
- Distance-based contacts do not confirm interaction energetics.
- ML predictions depend on the training dataset and applicability domain.
- Experimental validation is required.

## 10. Conclusion

Ten tryptanthrin-derived compounds were successfully docked
against TDO2 using a receptor containing the heme cofactor.

Reference ligand redocking achieved an RMSD of 1.488 Angstrom,
satisfying the selected validation threshold.

Three compounds were prioritized using docking scores,
ML predictions, drug-likeness and preliminary interaction analysis.

TRY_91b948175b25d8f6 was the leading computational candidate.

Experimental studies are required to confirm its inhibitory activity.
