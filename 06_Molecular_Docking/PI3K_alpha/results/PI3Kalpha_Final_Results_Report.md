# PI3Kα (4JPS) — Integrated docking results

**Target:** PIK3CA · **Structure:** 4JPS · **Reference ligand:** alpelisib (1LT)

## Reference redocking and protocol

- Symmetry-aware crystallographic-coordinate heavy-atom RMSD: **0.6039 Å**
- AutoDock Vina 1.2.7, exhaustiveness 32, num_modes 9, seed 42
- Grid center (Å): [-1.319, -9.513, 16.948]; grid dimensions (Å): [18.946, 19.659, 25.114]
- Receptor preparation limitation: 10 incomplete/excess-bond residues excluded >20 Å outside the provisional box; original 85 pocket heavy atoms retained at 0.0 Å displacement (prior verification)

## Ten-analogue results

| Rank | Analog_ID | Vina score (kcal/mol) | ML consensus pIC50 | Max train Tanimoto |
|---:|---|---:|---:|---:|
| 1 | `TRY_6db5e7a7802259ab` | -10.290 | 6.29 | 0.410 |
| 2 | `TRY_c88335488290afb4` | -10.261 | 6.25 | 0.403 |
| 3 | `TRY_c6b27d68ce6ec718` | -10.175 | 6.37 | 0.413 |
| 4 | `TRY_6bf367653fed7969` | -9.473 | 6.27 | 0.403 |
| 5 | `TRY_c2155e898021da1c` | -9.383 | 6.44 | 0.400 |
| 6 | `TRY_bafc3535a3d34b71` | -9.288 | 6.22 | 0.418 |
| 7 | `TRY_40b0c0e6de99494a` | -9.280 | 6.34 | 0.433 |
| 8 | `TRY_f751442f427b1d4d` | -9.170 | 6.30 | 0.400 |
| 9 | `TRY_e514ff1b3f600530` | -8.910 | 6.22 | 0.452 |
| 10 | `TRY_783966fbcbaa4f5c` | -8.716 | 6.22 | 0.403 |

## Top-three pose analysis

- **TRY_6db5e7a7802259ab**: 15 residues and 62 heavy-atom pairs within 4 Å; minimum distance 3.226 Å; pairs below 2 Å: 0.
- **TRY_c88335488290afb4**: 13 residues and 60 heavy-atom pairs within 4 Å; minimum distance 2.985 Å; pairs below 2 Å: 0.
- **TRY_c6b27d68ce6ec718**: 14 residues and 59 heavy-atom pairs within 4 Å; minimum distance 2.745 Å; pairs below 2 Å: 0.

## Scientific interpretation and limitations

- Docking scores are computational rankings for this receptor only, not measured affinity or confirmed cancer inhibition.
- Contact counts and polar atom distances do not establish hydrogen bonds; donor/acceptor chemistry and angles need separate validation.
- The top-three compounds have limited similarity (~0.40–0.41) to training structures; ML extrapolation remains uncertain.
- QED, SA, PAINS and Brenk information is inherited from the upstream screening table and does not confirm drug safety or synthesizability.
- Ten distant problematic receptor residues were excluded during preparation; document this and the engineered receptor mutations.
- Small differences between Vina scores should not be considered significant without further uncertainty analysis.

## Output files

- `PI3Kalpha_final_candidate_ranking.csv` — all ten docked analogues with ML/QED/SA and caution flags.
- `PI3Kalpha_top3_integrated_summary.csv` — three top poses, contact counts and minima.
- `PI3Kalpha_final_docking_ranking.png` — docking score figure.
- `source_visualizations/` — Step 9 pose files, HTML views and figures.
- Earlier Step 6–9 artifacts are needed for more complete provenance and end-to-end reproduction.
