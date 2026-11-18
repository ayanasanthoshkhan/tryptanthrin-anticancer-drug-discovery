# Step 12 — Fingerprint Collision Sensitivity

This is a read-only, post-hoc diagnostic using existing fixed-test predictions.
No models were retrained, no predictions were recalculated, and no source data were modified.

## Results

| Target | Model | Affected test molecules | Original test R2 | Without exact train-FP overlap R2 | Original MAE | Without overlap MAE |
|---|---|---:|---:|---:|---:|---:|
| EGFR | RF | 38/1282 | 0.6287 | 0.6289 | 0.6192 | 0.6216 |
| EGFR | XGB | 38/1282 | 0.6339 | 0.6354 | 0.6083 | 0.6101 |
| VEGFR2 | RF | 17/1254 | 0.5590 | 0.5576 | 0.5964 | 0.5981 |
| VEGFR2 | XGB | 17/1254 | 0.5562 | 0.5545 | 0.5890 | 0.5907 |
| PI3Kalpha | RF | 24/970 | 0.6191 | 0.6152 | 0.5491 | 0.5507 |
| PI3Kalpha | XGB | 24/970 | 0.6610 | 0.6569 | 0.5135 | 0.5162 |
| IDO1 | RF | 17/435 | 0.6244 | 0.6283 | 0.6341 | 0.6379 |
| IDO1 | XGB | 17/435 | 0.6326 | 0.6379 | 0.6464 | 0.6494 |
| TDO2 | RF | 2/164 | 0.5051 | 0.5054 | 0.4576 | 0.4556 |
| TDO2 | XGB | 2/164 | 0.5565 | 0.5554 | 0.4237 | 0.4226 |

## Interpretation

- Counts of affected **test molecules** can exceed counts of distinct overlapping fingerprint patterns.
- Equal 2048-bit ECFP4 input vectors can occur for distinct structures; this is not automatically a canonical-structure duplicate.
- The collision-excluded test subset has a different composition and may have different activity variance. A changed R2 does **not** by itself quantify bias caused by overlap.
- The collision-only subset may be very small. Treat its scores as descriptive, not robust.
- The original fixed test metrics remain the primary reported baseline, accompanied by this sensitivity check and scaffold-disjointness limitations.
- This analysis does not assess stereochemistry-aware fingerprints, tautomer/salt standardization, assay heterogeneity, calibration, or performance on externally measured triptanthrin analogs.
- Avoid selecting or tuning models against these post-hoc test metrics.
