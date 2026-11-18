"""Step 5: reproducible scaffold-split Random Forest + XGBoost for 4 ChEMBL targets.

Place in the existing project's scripts/ directory. Run from project root:
    .\\.venv\\Scripts\\python.exe scripts/05_train_remaining_targets.py

Read-only: data/raw/*.csv, existing EGFR models/reports. Does not incorporate
BindingDB or generated analog predictions. Targets are *median per-assay pIC50*.
"""
from pathlib import Path
import csv
import json
import warnings

import joblib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import xgboost as xgb

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
REPORTS = ROOT / 'reports'
MODELS = ROOT / 'models'
PROCESSED = ROOT / 'data' / 'processed'
TARGETS = ['VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2']
SEED = 42
FP_BITS = 2048


def metrics(y, predictions):
    return {'n': int(len(y)), 'MAE': float(mean_absolute_error(y, predictions)),
            'RMSE': float(np.sqrt(mean_squared_error(y, predictions))),
            'R2': float(r2_score(y, predictions))}


def load_target(target):
    path = RAW / f'{target}_Final_ML_Dataset.csv'
    if not path.is_file():
        raise FileNotFoundError(f'Missing input: {path}. Check original Step 1 data/raw folder.')
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    needed = {'molecule_chembl_id', 'canonical_smiles', 'pIC50'}
    if not rows or not needed.issubset(rows[0]):
        raise ValueError(f'{target}: empty file or missing columns {needed}')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_BITS)
    X = np.empty((len(rows), FP_BITS), dtype=np.uint8)
    y = np.empty(len(rows), dtype=np.float64)
    scaffolds, ids, smiles = [], [], []
    for i, row in enumerate(rows):
        s = row['canonical_smiles'].strip()
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise ValueError(f'{target}: invalid SMILES at row {i}: {s}')
        fp = generator.GetFingerprintAsNumPy(mol)
        X[i] = fp
        y[i] = float(row['pIC50'])
        ids.append(row['molecule_chembl_id'].strip())
        smiles.append(s)
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffolds.append(Chem.MolToSmiles(core) or '__NO_MURCKO_SCAFFOLD__')
    if not np.isfinite(y).all() or len(set(ids)) != len(ids):
        raise ValueError(f'{target}: nonfinite pIC50 or duplicated molecule IDs; review data')
    # No salts/tautomer canonicalization at this stage, consistent with Steps 2-4.
    return X, y, np.asarray(scaffolds), ids, smiles


def scaffold_split(X, y, scaffolds):
    first = GroupShuffleSplit(n_splits=1, test_size=0.30, random_state=SEED)
    train, hold = next(first.split(X, y, groups=scaffolds))
    second = GroupShuffleSplit(n_splits=1, test_size=0.50, random_state=SEED)
    rel_val, rel_test = next(second.split(X[hold], y[hold], groups=scaffolds[hold]))
    indices = {'train': train, 'validation': hold[rel_val], 'test': hold[rel_test]}
    sets = {name: set(scaffolds[idx]) for name, idx in indices.items()}
    assert len(set(np.concatenate(list(indices.values())))) == len(y), 'Incomplete split'
    assert not (sets['train'] & sets['validation'] or sets['train'] & sets['test']
                or sets['validation'] & sets['test']), 'Scaffold leakage'
    assert all(len(idx) > 1 for idx in indices.values()), 'Split too small'
    return indices


def save_plot(actual, predicted, output, target, model_name):
    lo = min(actual.min(), predicted.min()); hi = max(actual.max(), predicted.max())
    fig, ax = plt.subplots(figsize=(5.7, 5.2))
    ax.scatter(actual, predicted, s=11, alpha=0.4)
    ax.plot([lo, hi], [lo, hi], '--', linewidth=1)
    ax.set(xlabel='Experimental median pIC50', ylabel='Predicted pIC50',
           title=f'{target} {model_name}: scaffold-held-out test')
    fig.tight_layout()
    fig.savefig(output / 'test_actual_vs_predicted.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.7, 4))
    ax.hist(actual - predicted, bins=35)
    ax.set(xlabel='Residual (actual - predicted)', ylabel='Count',
           title=f'{target} {model_name}: test residuals')
    fig.tight_layout()
    fig.savefig(output / 'test_residuals.png', dpi=160)
    plt.close(fig)


def save_result(target, model_name, model, X, y, ids, indices, metrics_by_split,
                output, flagged, extra):
    output.mkdir(parents=True, exist_ok=True)
    model_label = 'RF' if model_name == 'Random Forest' else 'XGB'
    model_path = MODELS / f'{target}_{model_label}_scaffold_baseline.joblib'
    joblib.dump(model, model_path)
    test = indices['test']
    predicted = model.predict(X[test])
    with (output / 'test_predictions.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['fingerprint_row', 'molecule_chembl_id', 'actual_pIC50',
                         'predicted_pIC50', 'residual_actual_minus_predicted',
                         'extreme_review_flag'])
        for i, p in zip(test, predicted):
            writer.writerow([int(i), ids[i], float(y[i]), float(p), float(y[i] - p), bool(flagged[i])])
    save_plot(y[test], predicted, output, target, model_name)
    report = {'target': target, 'model': model_name,
              'target_label': 'ChEMBL median of individual per-measurement pIC50',
              'fingerprints': 'RDKit Morgan ECFP4 radius 2, 2048 bits',
              'split': 'two-stage GroupShuffleSplit grouped by Murcko scaffold; seed=42',
              'metrics': metrics_by_split,
              'test_extreme_review_flags': int(flagged[test].sum()),
              'test_without_extreme_review_flags': (metrics(y[test][~flagged[test]], predicted[~flagged[test]])
                  if (~flagged[test]).sum() >= 2 else None),
              'rdkit_version': rdBase.rdkitVersion,
              'sklearn_version': sklearn_version,
              'xgboost_version': xgb.__version__,
              'extra': extra,
              'limitations': ['No salt/tautomer standardization',
                'Extreme pIC50 values retained; review required',
                'A scaffold split does not prove validity on triptanthrin analogs',
                'Do not tune using the test scores']}
    (output / 'metrics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return model_path


def train_one(target):
    print(f'\n========== {target} ==========', flush=True)
    X, y, scaffolds, ids, smiles = load_target(target)
    indices = scaffold_split(X, y, scaffolds)
    flagged = (y < 3) | (y > 12)  # Existing Step 3 review heuristic, NOT removal.
    print(f'Molecules: {len(y):,}; scaffolds: {len(set(scaffolds)):,}; '
          + ', '.join(f'{name}={len(idx):,}' for name, idx in indices.items()), flush=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PROCESSED / f'{target}_ECFP4_2048.npz', X=X, y=y)
    with (PROCESSED / f'{target}_fingerprint_metadata.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['fingerprint_row', 'molecule_chembl_id', 'canonical_smiles', 'pIC50'])
        writer.writerows((i, ids[i], smiles[i], float(y[i])) for i in range(len(y)))
    split_file = REPORTS / target.lower() / 'split_assignments.csv'
    split_file.parent.mkdir(parents=True, exist_ok=True)
    with split_file.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['fingerprint_row', 'molecule_chembl_id', 'scaffold_smiles', 'split', 'pIC50', 'outlier_review_flag'])
        for group, idx in indices.items():
            writer.writerows((int(i), ids[i], scaffolds[i], group, float(y[i]), bool(flagged[i])) for i in idx)

    train, val, test = (indices[key] for key in ('train', 'validation', 'test'))
    comparisons = []
    baseline = metrics(y[test], np.full(len(test), float(np.median(y[train]))))
    print(f'Test median baseline: MAE={baseline["MAE"]:.4f} RMSE={baseline["RMSE"]:.4f} R2={baseline["R2"]:.4f}')
    print('Training Random Forest (200 trees)...', flush=True)
    rf = RandomForestRegressor(n_estimators=200, min_samples_leaf=2,
                               max_features='sqrt', n_jobs=4, random_state=SEED)
    rf.fit(X[train], y[train])
    rf_scores = {part: metrics(y[idx], rf.predict(X[idx])) for part, idx in indices.items()}
    save_result(target, 'Random Forest', rf, X, y, ids, indices, rf_scores,
                REPORTS / target.lower() / 'random_forest', flagged,
                {'n_estimators': 200, 'min_samples_leaf': 2,
                 'test_training_median_baseline': baseline})
    for part, sc in rf_scores.items():
        comparisons.append([target, 'Random Forest', part, *[sc[k] for k in ('n', 'MAE', 'RMSE', 'R2')]])
    print('Training XGBoost (validation early stopping)...', flush=True)
    gb = xgb.XGBRegressor(objective='reg:squarederror', n_estimators=1200,
        learning_rate=0.04, max_depth=5, min_child_weight=3,
        subsample=0.85, colsample_bytree=0.85, reg_lambda=5.0,
        reg_alpha=0.1, tree_method='hist', n_jobs=4,
        random_state=SEED, eval_metric='rmse', early_stopping_rounds=40)
    gb.fit(X[train], y[train], eval_set=[(X[val], y[val])], verbose=False)
    xgb_scores = {part: metrics(y[idx], gb.predict(X[idx])) for part, idx in indices.items()}
    save_result(target, 'XGBoost', gb, X, y, ids, indices, xgb_scores,
                REPORTS / target.lower() / 'xgboost', flagged,
                {'best_iteration_zero_based': int(gb.best_iteration),
                 'test_training_median_baseline': baseline})
    for part, sc in xgb_scores.items():
        comparisons.append([target, 'XGBoost', part, *[sc[k] for k in ('n', 'MAE', 'RMSE', 'R2')]])
    for name, scores in [('RF', rf_scores), ('XGB', xgb_scores)]:
        r = scores['test']
        print(f'{name:3s} TEST n={r["n"]:5d} MAE={r["MAE"]:.4f} RMSE={r["RMSE"]:.4f} R2={r["R2"]:.4f}', flush=True)
    print(f'Test extreme-value flags: {int(flagged[test].sum())} (not deleted)', flush=True)
    return comparisons


def add_existing_egfr(rows):
    rf_path = REPORTS / 'egfr_random_forest' / 'metrics.json'
    xgb_path = REPORTS / 'egfr_xgboost' / 'metrics.json'
    for path, label, key in [(rf_path, 'Random Forest', 'scores'),
                             (xgb_path, 'XGBoost', 'metrics')]:
        if not path.exists():
            print(f'EGFR {label} report not found: leaving it out of combined CSV')
            continue
        data = json.loads(path.read_text(encoding='utf-8'))
        for part in ('train', 'validation', 'test'):
            m = data[key][part]
            rows.append(['EGFR', label, part, int(m['n']), float(m['MAE']),
                         float(m['RMSE']), float(m['R2'])])


def main():
    MODELS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    rows = []
    add_existing_egfr(rows)
    for target in TARGETS:
        rows.extend(train_one(target))
    dest = REPORTS / 'all_targets_model_comparison.csv'
    with dest.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['target', 'model', 'split', 'n', 'MAE', 'RMSE', 'R2'])
        writer.writerows(rows)
    print(f'\nStep 5 complete. Combined results: {dest.relative_to(ROOT)}')
    print('Test metrics are preliminary; avoid repeated model selection based on test results.')


if __name__ == '__main__':
    main()
