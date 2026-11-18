# EGFR Molecular Docking

**Responsible member:** Nakshathra

## 1. Objective

To perform molecular docking of selected tryptanthrin analogues against the Epidermal Growth Factor Receptor (EGFR) and prioritise potential candidates for further investigation.

**Target:** EGFR  
**PDB structure:** 1M17  
**Chain:** A  
**Docking software:** AutoDock Vina 1.2.7  
**Ligand preparation:** RDKit and Meeko

## 2. Workflow

1. Retrieved the EGFR crystal structure (PDB: 1M17).
2. Prepared and cleaned the receptor structure.
3. Identified the reference ligand binding-site coordinates.
4. Generated 3D conformers for the selected analogues using RDKit.
5. Prepared ligand and receptor PDBQT files using Meeko.
6. Performed molecular docking using AutoDock Vina.
7. Ranked candidates based on docking scores.
8. Examined selected docked poses and residue-level contacts.

## 3. Ligand Shortlist

Ten compounds were selected from the ML-prioritised shortlist for EGFR docking.

All ten ligands were successfully prepared and docked.

## 4. Docking Results

| Compound ID | Vina score (kcal/mol) |
|---|---:|
| TRY_d75d8cfc4947cbd8 | -9.278 |
| TRY_0268c2efe76f90a9 | -9.098 |
| TRY_97b89eb811badd69 | -9.075 |
| TRY_4abc3bf0799ed766 | -9.052 |
| TRY_2327d14018ee1ca2 | -9.014 |
| TRY_fed8d7044fa845ad | -8.920 |
| TRY_51906f2bc26ac156 | -8.802 |
| TRY_01fe40c9b5aea27c | -8.669 |
| TRY_2632ae4667cd5093 | -8.634 |
| TRY_0328953f8b3ac4a0 | -8.256 |

The most favourable docking score was obtained by **TRY_d75d8cfc4947cbd8 (-9.278 kcal/mol)**.

## 5. Preliminary Pose Analysis

The top three candidates were examined for residue-level contacts within 4 Å.

- TRY_d75d8cfc4947cbd8: contacts with six unique residues.
- TRY_97b89eb811badd69: contacts with four unique residues.
- TRY_0268c2efe76f90a9: contacts with one unique residue.

These are distance-based contacts and should not be interpreted as confirmed hydrogen bonds or experimentally validated interactions.

## 6. Preliminary Candidate

TRY_d75d8cfc4947cbd8 is the leading computational candidate based on its docking score and observed residue-contact profile.

This is a computational prioritisation result, not evidence of confirmed EGFR inhibition or anticancer activity.

## 7. Important Validation Limitation

The initial reference-ligand redocking exercise produced a reported RMSD of 1.716 Å. However, the extracted ligand's inferred chemical connectivity did not match the expected identity of erlotinib.

Therefore, the redocking result is considered **provisional and not valid confirmation of the docking protocol**. Reference-ligand identity and connectivity must be corrected before claiming successful protocol validation.

## 8. Folder Structure

- `notebooks/` – EGFR docking workflow notebook.
- `structures/` – Receptor and reference-ligand structures.
- `ligands/` – Prepared analogue structures in SDF and PDBQT formats.
- `results/` – Docking scores, rankings, docked poses and contact-analysis files.
- `analysis/` – Additional analysis and validation outputs.

## 9. Next Steps

- Resolve reference-ligand identity and redocking validation.
- Review binding interactions using appropriate molecular visualisation tools.
- Compare EGFR results with the other selected anticancer targets.
- Integrate docking results with ADMET and ML prioritisation.
- Select candidates for further investigation.

---

*All docking scores are computational predictions and require further validation.*
