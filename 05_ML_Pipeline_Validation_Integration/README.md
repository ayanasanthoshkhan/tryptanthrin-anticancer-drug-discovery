# Tryptanthrin — ML Screening, Validation and Integration

## Team handover (results-review archive)

This archive contains the completed computational ML workstream of the tryptanthrin anticancer research project: training and validation scripts, an executed results-review notebook, evaluation figures and reports, an integrated prediction table, and a preliminary candidate shortlist.

### Modelled targets
EGFR, VEGFR2, PI3Kα, IDO1 and TDO2. The original research work plan mentioned six targets; confirm scope with the team.

### Start here
1. Read `reports/member5_submission/Tryptanthrin_ML_Final_Report_REVIEWED.md`.
2. Review `reports/member5_submission/final_metrics_summary.csv` and the three figures in that folder.
3. Inspect `reports/analog_predictions/all_target_predictions.csv` (56,775 analog–target predictions).
4. Inspect `reports/diversity_shortlist/master_docking_handover.csv` (50 target-specific selections; 49 distinct Analog IDs).
5. Open `notebooks/Tryptanthrin_ML_results_review.ipynb`. It reads saved CSV/PNG results and does **not** retrain the models. Run from the handover root or `notebooks/` folder with Jupyter, IPython and matplotlib installed.

### Reproducibility limits
- Source ChEMBL datasets, the filtered input analog library, training fingerprints and trained model binaries are **not included** in this lightweight WhatsApp ZIP. Training and prediction scripts therefore will NOT run end-to-end from this archive alone. Obtain the source datasets and model artifacts from the team / original working project before attempting reproduction.
- `requirements.txt` is the original starter requirements file and is **not a complete environment lockfile**. Notebook review may need `matplotlib`, `ipython`, `ipykernel`; complete training requires additional libraries used by the scripts (e.g., scikit-learn, xgboost, pandas, numpy and RDKit).
- `reports/all_targets_model_comparison.csv` in this archive is a **review convenience table reconstructed from** the saved `reports/member5_submission/final_metrics_summary.csv` (test metrics only). It is not a regenerated model evaluation.
- Historical internal `member5_*` report folder names and some source script names are retained to preserve generated workflow references. Presentation-facing report titles and filenames use project terminology.

### Scientific limitations
Predictions are in-silico hypotheses, not experimentally validated anticancer activity. Structural applicability is especially limited for EGFR, VEGFR2 and PI3Kα; TDO2 shows weaker robustness under alternative scaffold splits. Review structural alerts and verify published references and parent-compound comparisons before selecting candidates for further studies.
