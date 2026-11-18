# PI3Kα — Step 11 Interaction Geometry and Provisional Lead Review

Reference: 4JPS; the receptor and analogue docked poses are unchanged from the original handover.
This analysis uses explicit polar hydrogens retained in PDBQT and no coordinate alignment.

## Criteria
Hydrogen-bond geometry: D···A <=3.5 Å, H···A <=2.7 Å, D–H···A >=120°. Proximity is not experimental confirmation.
Hydrophobic contacts: ligand/protein carbon-pair distance <=4.5 Å; preliminary proximity only.
Aromatic proximity: ligand aromatic-typed atom within 5.0 Å of aromatic-residue side-chain carbon; π-stacking orientation not established.

## Provisional shortlist
- TRY_6db5e7a7802259ab: Vina -10.290 kcal/mol; consensus predicted pIC50 6.294; QED 0.541; SA 2.49; maximum training Tanimoto 0.410; geometrically plausible D–H···A contacts 0.
- TRY_c88335488290afb4: Vina -10.261 kcal/mol; consensus predicted pIC50 6.254; QED 0.556; SA 2.44; maximum training Tanimoto 0.403; geometrically plausible D–H···A contacts 0.
- TRY_c6b27d68ce6ec718: Vina -10.175 kcal/mol; consensus predicted pIC50 6.366; QED 0.560; SA 2.51; maximum training Tanimoto 0.413; geometrically plausible D–H···A contacts 2.

## Limitations
- Predictions are computational and do not establish PI3Kα inhibition or anticancer activity.
- Vina score differences are small and not equivalent to experimental affinity differences.
- Protein/ligand protonation, hydrogen assignment and tautomeric states may alter hydrogen-bond geometry.
- Carbon proximity does not establish favourable hydrophobic energy; aromatic geometry needs centroid/plane-level review.
- Low similarity to ML training structures limits confidence in predicted pIC50.
- Further structural checks, chemical synthesis feasibility and experimental testing are needed.
- Existing original workflow validation RMSD was 0.6039 Å; no receptor or docking parameters were altered in this step.
