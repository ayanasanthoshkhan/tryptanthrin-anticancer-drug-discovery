# Step 13 — Repeated scaffold-split robustness (Random Forest)

Run UTC: 2026-11-15T15:43:02.680321+00:00
RDKit: 2026.09.2; scikit-learn: 1.9.1
Targets: EGFR, VEGFR2, PI3Kalpha, IDO1, TDO2
New seeds: [13, 23, 37]; n_estimators: 200
Same RF settings as baseline: min_samples_leaf=2, max_features=sqrt, random_state=42; CPU n_jobs=4.
Scaffold-disjoint split: 70/15/15 approximate, two-stage GroupShuffleSplit by Murcko scaffold; exact sample proportions can vary.
Each model fits ONLY its own new train partition; validation/test labels are not used in fitting.
Original seed-42 saved models, test predictions, candidate rankings, and docking handover were untouched.

## Findings

| Target | Repeats | Mean test MAE ± SD | Mean test RMSE ± SD | Mean test R² ± SD | Test R² range |
|---|---:|---:|---:|---:|---:|
| EGFR | 3 | 0.6714 ± 0.0187 | 0.8847 ± 0.0266 | 0.5538 ± 0.0113 | 0.5450–0.5666 |
| VEGFR2 | 3 | 0.6030 ± 0.0221 | 0.7936 ± 0.0297 | 0.5070 ± 0.0391 | 0.4719–0.5491 |
| PI3Kalpha | 3 | 0.5530 ± 0.0197 | 0.7581 ± 0.0248 | 0.6069 ± 0.0408 | 0.5598–0.6308 |
| IDO1 | 3 | 0.6347 ± 0.0139 | 0.8076 ± 0.0275 | 0.5720 ± 0.0253 | 0.5435–0.5919 |
| TDO2 | 3 | 0.5255 ± 0.0239 | 0.6986 ± 0.0114 | 0.3389 ± 0.0546 | 0.2980–0.4009 |

## Interpretation and limitations

- These are robustness refits on alternative scaffold splits, not external or prospective experimental validation.
- Only Random Forest was repeated; do not report XGBoost as robustly revalidated from this script.
- Different runs have different test structures, sizes and activity distributions, so variability includes split composition.
- Identical fingerprint vectors may occur across distinct structures even in scaffold-disjoint splits; overlap counts are reported.
- The model is trained exclusively with target-specific ChEMBL labels; assay heterogeneity and salts/tautomers are not resolved.
- Three repeats illustrate variation but do not yield a calibrated uncertainty interval.
- Retain original seed-42 metrics as baseline; avoid model tuning or ranking based on any new held-out test scores.
- Structural applicability to generated triptanthrin molecules remains a separate, major concern.
