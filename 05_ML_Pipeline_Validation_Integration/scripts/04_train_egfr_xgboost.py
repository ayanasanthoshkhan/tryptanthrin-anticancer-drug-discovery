"""Step 4 — EGFR XGBoost regression, using Step 3's EXACT scaffold split.

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts/04_train_egfr_xgboost.py

Read-only inputs: Step 2 NPZ and metadata; Step 3 split assignments and metrics.
The original curated datasets and Random Forest model are never changed.
"""
from pathlib import Path
import csv
import json

import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import xgboost as xgb
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / 'data' / 'processed' / 'EGFR_ECFP4_2048.npz'
META = ROOT / 'data' / 'processed' / 'EGFR_fingerprint_metadata.csv'
RF_DIR = ROOT / 'reports' / 'egfr_random_forest'
SPLITS = RF_DIR / 'split_assignments.csv'
RF_METRICS = RF_DIR / 'metrics.json'
OUT = ROOT / 'reports' / 'egfr_xgboost'
MODEL = ROOT / 'models' / 'EGFR_XGB_scaffold_baseline.joblib'
SEED = 42


def load_split(metadata, y):
    with SPLITS.open(newline='', encoding='utf-8-sig') as file:
        assignments = list(csv.DictReader(file))
    if len(assignments) != len(y):
        raise ValueError('Step 3 split assignments do not cover every fingerprint row.')

    groups = {'train': [], 'validation': [], 'test': []}
    scaffolds = {name: set() for name in groups}
    visited = set()
    for record in assignments:
        i = int(record['fingerprint_row'])
        group = record['split']
        if i < 0 or i >= len(y) or i in visited or group not in groups:
            raise ValueError('Invalid or duplicate split assignment.')
        visited.add(i)
        if record['molecule_chembl_id'] != metadata[i]['molecule_chembl_id']:
            raise ValueError(f'Molecule-ID mismatch at fingerprint row {i}.')
        if not np.isclose(float(record['pIC50']), float(y[i]), atol=1e-6):
            raise ValueError(f'Activity label mismatch at row {i}.')
        groups[group].append(i)
        scaffolds[group].add(record['scaffold_smiles'])
    if visited != set(range(len(y))) or any(len(v) == 0 for v in groups.values()):
        raise ValueError('Incomplete or empty train/validation/test split.')
    if any(scaffolds[a] & scaffolds[b] for a, b in
           [('train', 'validation'), ('train', 'test'), ('validation', 'test')]):
        raise ValueError('Detected overlapping scaffolds between splits.')
    return {key: np.asarray(indices, dtype=int) for key, indices in groups.items()}


def scores(actual, predicted):
    return {
        'n': int(len(actual)),
        'MAE': float(mean_absolute_error(actual, predicted)),
        'RMSE': float(np.sqrt(mean_squared_error(actual, predicted))),
        'R2': float(r2_score(actual, predicted)),
    }


def main():
    for path in [FEATURES, META, SPLITS, RF_METRICS]:
        if not path.is_file():
            raise FileNotFoundError(f'Missing: {path}. Complete Steps 2 and 3 first.')
    with np.load(FEATURES) as contents:
        X = contents['X']
        y = contents['y']
    with META.open(newline='', encoding='utf-8-sig') as file:
        metadata = list(csv.DictReader(file))
    if not (X.ndim == 2 and X.shape[1] == 2048 and X.shape[0] == len(y) == len(metadata)):
        raise ValueError('Unexpected feature/label/metadata dimensions.')
    if not np.isfinite(y).all():
        raise ValueError('Non-finite pIC50 labels.')
    if any(int(row['fingerprint_row']) != i for i, row in enumerate(metadata)):
        raise ValueError('Fingerprint row indices are not aligned.')
    split = load_split(metadata, y)
    train, val, test = (split[name] for name in ('train', 'validation', 'test'))
    rf_report = json.loads(RF_METRICS.read_text(encoding='utf-8'))
    if any(rf_report['counts'][name] != len(split[name]) for name in split):
        raise ValueError('Random Forest metrics and scaffold splits disagree.')

    print('EGFR XGBoost — same exact scaffold split as Step 3')
    print(f'Molecules: {len(y):,} | train={len(train):,} validation={len(val):,} test={len(test):,}')
    print('Training XGBoost with validation early stopping... ', flush=True)
    model = xgb.XGBRegressor(
        objective='reg:squarederror',
        n_estimators=1200, learning_rate=0.04, max_depth=5,
        min_child_weight=3, subsample=0.85, colsample_bytree=0.85,
        reg_lambda=5.0, reg_alpha=0.1,
        tree_method='hist', n_jobs=4, random_state=SEED,
        eval_metric='rmse', early_stopping_rounds=40,
    )
    model.fit(X[train], y[train], eval_set=[(X[val], y[val])], verbose=False)
    print(f'Best boosting iteration (zero-based): {model.best_iteration}')

    results = {}
    predictions = {}
    for name, indices in split.items():
        pred = model.predict(X[indices])
        predictions[name] = pred
        results[name] = scores(y[indices], pred)
        r = results[name]
        print(f'{name:10s} n={r["n"]:5d} MAE={r["MAE"]:.4f} RMSE={r["RMSE"]:.4f} R2={r["R2"]:.4f}')

    rf_test = rf_report['scores']['test']
    xgb_test = results['test']
    print('\nHeld-out TEST comparison (same molecule indices):')
    print(f'Random Forest: MAE={rf_test["MAE"]:.4f} RMSE={rf_test["RMSE"]:.4f} R2={rf_test["R2"]:.4f}')
    print(f'XGBoost:       MAE={xgb_test["MAE"]:.4f} RMSE={xgb_test["RMSE"]:.4f} R2={xgb_test["R2"]:.4f}')

    OUT.mkdir(parents=True, exist_ok=True)
    MODEL.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, MODEL)
    report = {
        'target': 'EGFR',
        'labels': 'stored molecule-level median pIC50 (median of per-record pIC50)',
        'fingerprints': 'Morgan radius 2 / 2048 bits, original unstandardized SMILES',
        'split': 'exact reuse of Step 3 split_assignments.csv',
        'seed': SEED, 'xgboost_version': xgb.__version__,
        'best_iteration_zero_based': int(model.best_iteration),
        'metrics': results,
        'rf_test_metrics': rf_test,
        'notes': [
            'Validation used for early stopping; test never used in fitting.',
            'Scaffold holdout does not establish reliability on triptanthrin analogs.',
            'Outlier flags were not excluded.',
            'Any later hyperparameter tuning based on test results would invalidate the test as untouched holdout.',
        ],
    }
    (OUT / 'metrics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (OUT / 'test_predictions.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['fingerprint_row', 'molecule_chembl_id', 'actual_pIC50',
                         'predicted_pIC50', 'residual_actual_minus_predicted'])
        for i, pred in zip(test, predictions['test']):
            writer.writerow([int(i), metadata[i]['molecule_chembl_id'],
                             float(y[i]), float(pred), float(y[i] - pred)])
    with (OUT / 'model_comparison.csv').open('w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['model', 'dataset', 'n', 'MAE', 'RMSE', 'R2'])
        for name, data in [('Random Forest', rf_report['scores']), ('XGBoost', results)]:
            for split_name in ('train', 'validation', 'test'):
                metric = data[split_name]
                writer.writerow([name, split_name, metric['n'], metric['MAE'],
                                 metric['RMSE'], metric['R2']])

    actual = y[test]
    pred = predictions['test']
    low = min(actual.min(), pred.min())
    high = max(actual.max(), pred.max())
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(actual, pred, alpha=0.35, s=12)
    ax.plot([low, high], [low, high], linestyle='--', linewidth=1)
    ax.set(xlabel='Experimental median pIC50', ylabel='Predicted pIC50',
           title='EGFR XGBoost: scaffold-held-out test')
    fig.tight_layout()
    fig.savefig(OUT / 'test_actual_vs_predicted.png', dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(actual - pred, bins=40)
    ax.set(xlabel='Actual - predicted pIC50', ylabel='Number of molecules',
           title='EGFR XGBoost: test residuals')
    fig.tight_layout()
    fig.savefig(OUT / 'test_residuals.png', dpi=170)
    plt.close(fig)
    print(f'\nSaved model: {MODEL.relative_to(ROOT)}')
    print(f'Saved reports: {OUT.relative_to(ROOT)}')
    print('Step 4 complete. Compare scores, but do not tune repeatedly on the test set.')


if __name__ == '__main__':
    main()
