# Step 11 — ML Validation Audit

This is a **non-destructive** audit of pre-existing model inputs, scaffold assignments, and test predictions.
No models were tuned or selected based on the original fixed test results.

## Audit checks

- Verified raw dataset, fingerprint matrix, metadata, activity labels, and split assignments align by row and molecule ID.
- Recomputed Bemis–Murcko scaffolds; checked overlap between every pair of train, validation and test partitions.
- Checked exact RDKit canonical structure overlap and identical 2048-bit Morgan fingerprint overlap across partitions. **Different molecules can have identical folded fingerprints**: this is a representation-level collision, not necessarily exact chemical-structure leakage.
- Checked pIC50 distribution, test activity extrapolation, extreme-value flags, and performance by activity subgroup.
- No salts or tautomers were standardized: alternative representations may remain a source of unmeasured chemical equivalence.

## Summary

| Target | Audit | Train/test scaffold overlap | Structure overlap | FP overlap | Test extremes | Test outside train label range |
|---|---|---:|---:|---:|---:|---:|
| EGFR | REVIEW_FINGERPRINT_COLLISIONS | 0 | 0 | 38 | 0 | 0 |
| VEGFR2 | REVIEW_FINGERPRINT_COLLISIONS | 0 | 0 | 16 | 0 | 0 |
| PI3Kalpha | REVIEW_FINGERPRINT_COLLISIONS | 0 | 0 | 22 | 0 | 0 |
| IDO1 | REVIEW_FINGERPRINT_COLLISIONS | 0 | 0 | 16 | 2 | 1 |
| TDO2 | REVIEW_FINGERPRINT_COLLISIONS | 0 | 0 | 1 | 0 | 0 |

## Interpretation

- Continuous pIC50 prediction is regression; **classification class imbalance** is not the appropriate metric. The activity-bin table quantifies rare-label representation instead.
- Cross-target structural overlap is not in itself leakage: models have different endpoint targets. Contamination is checked within each target.
- **REVIEW_FINGERPRINT_COLLISIONS** indicates identical model inputs for chemically distinct molecules across splits. Review the collision-examples CSV; a chemically scaffold-disjoint split may still contain representation-level duplication.
- A clean split does not prove transferability to generated triptanthrin analogs; use the Step 6–9 structural applicability warnings.
- The pIC50 value is a curated median of per-measurement transformed activity. Raw IC50 median may differ under this nonlinear transform.
- Model agreement does not quantify calibrated uncertainty. No external validation is claimed.
- Test-subgroup slices are descriptive diagnostics, not independent evaluation sets.

## Next checks

- Confirm fifth vs sixth project target scope with teammates.
- Review source-unit consistency, molecule standardization, assay heterogeneity, and extreme activity records.
- Where feasible, use a separately sourced external experimental dataset for prospective validity.
- Document any selection or model hyperparameter tuning that used the same fixed test, and avoid further test-driven tuning.
