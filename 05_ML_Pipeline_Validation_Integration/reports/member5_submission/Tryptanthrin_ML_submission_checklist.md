# ML Screening and Validation — Submission checklist

Run Step 15 from project root; this is a file-integrity and completeness check, not a scientific certification.

- Step 12 maximum absolute ΔR² (collision-excluded minus full): **0.0054**
- Notebook: **Valid Jupyter structure: 11 cells, 5 code cells; EXECUTION NOT VERIFIED**
- Applicability plot regenerated with external legend: **YES**

## File inventory

- [x] ML scripts: `scripts`
- [x] saved trained models: `models`
- [x] original target metrics: `reports\all_targets_model_comparison.csv`
- [x] applicability summary: `reports\applicability_domain\similarity_summary.csv`
- [x] full analog predictions: `reports\analog_predictions\all_target_predictions.csv`
- [x] docking shortlist: `reports\diversity_shortlist\master_docking_handover.csv`
- [x] Step 11 leakage audit: `reports\final_validation\leakage_and_dataset_audit.csv`
- [x] Step 12 collision sensitivity: `reports\final_validation\fingerprint_collision_sensitivity_metrics.csv`
- [x] Step 13 robustness: `reports\final_validation\repeated_scaffold_robustness\rf_robustness_summary.csv`
- [x] Step 14 notebook: `notebooks\Member5_results_review.ipynb`

## Human checks before handover

- [ ] Open and run each code cell of `notebooks/Member5_results_review.ipynb` in the correct Python environment.
- [ ] Verify scientific terminology, figure captions, citations and assay metadata with Member 4.
- [ ] Confirm five versus six targets with the team and guide.
- [ ] Confirm the published references and parent triptanthrin activity comparison; do not claim experimental hits.
- [ ] Share report, figures, saved source code, full predictions and audited 50-entry handover.
- [ ] Do not include virtual environment folders or large raw data archives in GitHub by default.
