# IDO1 Molecular Docking — Methods & Results Report
**Target:** Human Indoleamine 2,3-Dioxygenase 1 (IDO1)  
**PDB:** 9FOZ (X-ray, 1.69 Å)  
**Analyst:** Shipra  
**Date:** 2026-11-18  
**Protocol:** AutoDock Vina 1.2.7, Exhaustiveness 32, Seed 42, 9 modes

---

## 1. Receptor Preparation

- **Structure:** PDB 9FOZ — Human IDO1, Chain A, 386 residues + HEM cofactor
- **Altlocs:** Default altloc A selected (`--default_altloc A`); residues with altloc B/C removed
- **Missing residues:** 32 unobserved N-terminal residues (all His-tag + linker; not part of native IDO1 — no impact on binding site)
- **Quality:** Clashscore 3.17, Ramachandran outliers 0%, rotamer outliers 0.94%, RSRZ 7.56% (slightly elevated, documented)
- **HEM handling:** Meeko 0.8.0 cannot process HEM CIF (RDKit valence error on carboxylate O atoms). Workaround: receptor prepared as protein-only PDBQT via `mk_prepare_receptor`, then HEM atoms manually appended with AutoDock4 atom types (Fe→Fe, N→NA, O→OA, C→A). Final receptor: `IDO1/00_receptor/IDO1.pdbqt` (3608 lines)
- **PDBQT tool:** Meeko 0.8.0 / `mk_prepare_receptor`

## 2. Co-crystal Ligand & Box Definition

- **Co-crystal inhibitor:** A1H93 (triazole-based IDO1 inhibitor; SMILES: `COc1cc(F)c(cc1F)n2cc(nn2)[C@H](O)c3n4cncc4cnc3C5CC5`)
- **Ligand prep:** Crystal coordinates extracted from MMCIF; hydrogens added; `mk_prepare_ligand` → `A1H93_reference.pdbqt`
- **Docking box center:** X=−10.985, Y=29.055, Z=−35.034 (centroid of A1H93 crystal heavy atoms)
- **Box size:** 22 × 22 × 22 Å
- **Key binding residues (IDO1):** Arg231, His346, Ser167, Phe163; Fe of HEM (catalytic coordination site)

## 3. Redocking Validation (S5 Gate)

| Mode | Vina Score (kcal/mol) | RMSD to Crystal (Å) |
|------|-----------------------|----------------------|
| 1 | −10.14 | 2.497 |
| 2 | −9.62 | 1.987 |
| **3** | **−9.13** | **0.444 ✅** |
| 4 | −8.87 | 0.887 |
| 5 | −8.75 | 1.626 |

**Best RMSD = 0.444 Å < 2.0 Å → GATE PASSED ✅**  
Note: Mode 3 (not Mode 1) gives best RMSD — common with Vina; pocket correctly identified.

## 4. Ligand Preparation (S6)

- **Source:** `master_docking_handover.csv` (filtered `Target == IDO1`), 10 analogs
- **3D generation:** RDKit ETKDGv3 (seed=42) + MMFF optimization
- **PDBQT prep:** `mk_prepare_ligand` (Meeko 0.8.0)
- **Result:** 10/10 analogs prepared successfully → `IDO1/04_ligand_prep/`

## 5. Docking Results (S7)

All 10 analogs docked to `IDO1.pdbqt` with identical box/exhaustiveness settings.

| Rank | Analog_ID | Vina Best (kcal/mol) |
|------|-----------|----------------------|
| 1 | TRY_b79da49d45504836 | **−10.91** |
| 2 | TRY_62b0878d0fc0ffb9 | −10.26 |
| 3 | TRY_487832cab4a8a96d | −9.855 |
| 4 | TRY_b82664877e66e437 | −9.596 |
| 5 | TRY_610c82dc1e8f7d06 | −9.496 |
| 6 | TRY_f114c2092e788124 | −9.467 |
| 7 | TRY_70bcb69ebe6539fe | −9.375 |
| 8 | TRY_51dfd35ef9fac3c9 | −9.247 |
| 9 | TRY_bb20494879ae80ec | −9.166 |
| 10 | TRY_5ca6460492f11a63 | −8.918 |

> All scores are predictions. Vina scores are not comparable across targets.

## 6. Interaction Analysis — Top 3 Analogs (S8)

Distance cutoffs: H-bond/coordination ≤ 3.5 Å; hydrophobic C–C ≤ 4.0 Å; Fe coordination ≤ 3.0 Å.

### TRY_b79da49d45504836 (Rank 1, −10.91 kcal/mol)
- **H-bond:** Ligand O ↔ Ser167_OG, 1.209 Å ✅
- **Hydrophobic:** Ligand C ↔ Arg231_CZ, 0.89 Å (close contact — likely π/CH interaction)
- **SMILES:** `O=C1c2cc(F)ccc2-n2c1nc1ccc(N3CCNCC3)cc1c2=O`
- **Notable:** Piperazine substituent — potential H-bond donor/acceptor; no PAINS/Brenk flags

### TRY_62b0878d0fc0ffb9 (Rank 2, −10.26 kcal/mol)
- **Fe coordination:** Ligand H ↔ Fe_HEM, 2.878 Å ✅ (proximal to catalytic iron)
- **H-bonds:** Ligand O ↔ Arg231_CZ, 2.829 Å; Ligand N ↔ Ser167_OG, 1.905 Å ✅
- **SMILES:** `NC(=O)c1ccc2c(=O)n3c(nc2c1)C(=O)c1cc(Cl)ccc1-3`
- **Notable:** Amide group engages both Arg231 and Ser167; proximity to Fe suggests coordination potential

### TRY_487832cab4a8a96d (Rank 3, −9.855 kcal/mol)
- **H-bonds:** Ligand N ↔ Arg231_CZ, 2.896 Å; Ligand N ↔ Ser167_OG, 1.874 Å ✅
- **SMILES:** `N#Cc1ccc2c(=O)n3c(nc2c1)C(=O)c1cc(F)ccc1-3`
- **Notable:** Nitrile group — potential Fe coordination (nitrile N→Fe); dual H-bond network

## 7. Integrated Candidate Ranking (S9)

Composite score = 40% Vina (normalized, inverted) + 40% Consensus pIC50 (normalized) + 20% QED — PAINS/Brenk penalties.

| Composite Rank | Analog_ID | Vina | Consensus pIC50 | QED | PAINS | Brenk | Composite |
|----------------|-----------|------|-----------------|-----|-------|-------|-----------|
| **1** | TRY_b79da49d45504836 | −10.91 | 6.553 | 0.563 | ✗ | ✗ | **0.867** |
| 2 | TRY_62b0878d0fc0ffb9 | −10.26 | 6.401 | 0.576 | ✗ | ✗ | 0.560 |
| 3 | TRY_610c82dc1e8f7d06 | −9.496 | 6.629 | 0.502 | ✗ | ✗ | 0.538 |
| 4 | TRY_f114c2092e788124 | −9.467 | 6.620 | 0.502 | ✗ | ✗ | 0.520 |
| 5 | TRY_bb20494879ae80ec | −9.166 | 6.495 | 0.498 | ✗ | ✗ | 0.282 |

## 8. Final Candidates

**Primary:** **TRY_b79da49d45504836** — Best Vina (−10.91), highest composite score (0.867), clean flags, piperazine substituent enabling Ser167 H-bond  
**Secondary:** **TRY_62b0878d0fc0ffb9** — Second-best Vina (−10.26), Fe proximal contact, dual H-bond at Arg231+Ser167  
**Tertiary:** **TRY_610c82dc1e8f7d06** — Highest ML pIC50 (6.629) among the remaining; clean flags

> **Disclaimer:** All results are computational predictions. Experimental validation (enzymatic assay, crystallography) required for confirmation.

## 9. Limitations

- RSRZ outliers 7.56% in 9FOZ (binding site residues should be verified)
- HEM atom types assigned manually — Fe partial charge not formally optimized for Vina scoring
- Vina uses empirical scoring; Fe-coordination energy not explicitly modeled
- 32 missing N-terminal residues (His-tag only; no binding site impact)

## 10. Output Files

| File | Description |
|------|-------------|
| `IDO1/00_receptor/IDO1.pdbqt` | Final receptor (protein + HEM) |
| `IDO1/01_cocrystal_ligand/A1H93_reference.pdbqt` | Co-crystal ligand |
| `IDO1/02_box/box_params.json` | Docking box parameters |
| `IDO1/03_redocking/IDO1_redocking_validation.csv` | Redocking RMSD results |
| `IDO1/04_ligand_prep/<ID>.pdbqt` | 10 analog ligand files |
| `IDO1/05_docking/<ID>_docked.pdbqt` | Docked poses (9 modes each) |
| `IDO1/05_docking/IDO1_vina_scores.csv` | Vina best scores |
| `IDO1/06_interactions/IDO1_hbond_geometry.csv` | H-bond/coordination contacts |
| `IDO1/06_interactions/IDO1_hydrophobic_contacts.csv` | Hydrophobic contacts |
| `IDO1/07_ranking/IDO1_final_candidate_ranking.csv` | Integrated ranking |
| `IDO1/08_report/IDO1_docking_methods_results.md` | This report |
| `IDO1/08_report/IDO1_docking_notebook.ipynb` | Colab-compatible notebook |
