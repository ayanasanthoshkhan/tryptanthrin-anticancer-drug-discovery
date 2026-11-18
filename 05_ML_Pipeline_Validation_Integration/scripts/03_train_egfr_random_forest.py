"""Step 3: EGFR Random Forest baseline with scaffold-disjoint holdout.

Run from anywhere: python scripts/03_train_egfr_random_forest.py
Inputs come from Step 2; original bioactivity datasets are never modified.
"""
from pathlib import Path
import csv
import json
import numpy as np
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
NPZ = ROOT / 'data' / 'processed' / 'EGFR_ECFP4_2048.npz'
META = ROOT / 'data' / 'processed' / 'EGFR_fingerprint_metadata.csv'
OUT = ROOT / 'reports' / 'egfr_random_forest'
MODELS = ROOT / 'models'
SEED = 42


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    with np.load(NPZ) as arr:
        X, y = arr['X'], arr['y']
    with META.open(encoding='utf-8-sig', newline='') as f:
        metadata = list(csv.DictReader(f))
    assert len(metadata) == len(y) == len(X), 'Feature/metadata rows do not align'
    assert X.shape[1] == 2048, 'Expect 2048-bit Step 2 fingerprints'
    assert all(int(r['fingerprint_row']) == i for i, r in enumerate(metadata)), 'Fingerprint row mapping mismatch'
    assert np.isfinite(y).all(), 'Nonfinite labels'

    scaffolds = []
    for row in metadata:
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is None:
            raise ValueError('Invalid input SMILES found')
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffolds.append(Chem.MolToSmiles(core) or '__NO_MURCKO_SCAFFOLD__')
    scaffolds = np.asarray(scaffolds)
    all_idx = np.arange(len(y))

    # ~70/15/15 partition by molecular scaffold; test held out throughout training.
    first = GroupShuffleSplit(n_splits=1, test_size=0.30, random_state=SEED)
    train_idx, hold_idx = next(first.split(X, y, groups=scaffolds))
    second = GroupShuffleSplit(n_splits=1, test_size=0.50, random_state=SEED)
    relative_val, relative_test = next(second.split(X[hold_idx], y[hold_idx], groups=scaffolds[hold_idx]))
    val_idx, test_idx = hold_idx[relative_val], hold_idx[relative_test]

    train_sc = set(scaffolds[train_idx]); val_sc = set(scaffolds[val_idx]); test_sc = set(scaffolds[test_idx])
    assert not (train_sc & val_sc or train_sc & test_sc or val_sc & test_sc), 'Scaffold leakage'
    assert len(set(train_idx) | set(val_idx) | set(test_idx)) == len(y)

    print(f'EGFR samples: {len(y):,} | distinct scaffolds: {len(set(scaffolds)):,}')
    print(f'Scaffold-disjoint split: train={len(train_idx):,}, validation={len(val_idx):,}, test={len(test_idx):,}')
    print('Training Random Forest (200 trees; this may take some time)...', flush=True)
    model = RandomForestRegressor(
        n_estimators=200, min_samples_leaf=2, max_features='sqrt',
        n_jobs=-1, random_state=SEED
    )
    model.fit(X[train_idx], y[train_idx])

    def metrics(indices, predictions):
        actual = y[indices]
        return {'n': len(indices), 'MAE': float(mean_absolute_error(actual, predictions)),
                'RMSE': float(np.sqrt(mean_squared_error(actual, predictions))),
                'R2': float(r2_score(actual, predictions))}

    scores = {}
    preds = {}
    for name, inds in [('train', train_idx), ('validation', val_idx), ('test', test_idx)]:
        pred = model.predict(X[inds])
        preds[name] = pred
        scores[name] = metrics(inds, pred)
        print(f'{name:10s} n={len(inds):5d} MAE={scores[name]["MAE"]:.4f} RMSE={scores[name]["RMSE"]:.4f} R2={scores[name]["R2"]:.4f}')

    baseline = np.full(len(test_idx), float(np.median(y[train_idx])))
    baseline_scores = metrics(test_idx, baseline)
    flagged = np.array([(float(row['pIC50']) < 3 or float(row['pIC50']) > 12) for row in metadata])
    non_outlier_test = ~flagged[test_idx]
    sensitivity = None
    if non_outlier_test.sum() >= 2:
        sensitivity = metrics(test_idx[non_outlier_test], preds['test'][non_outlier_test])
    print(f'Test median baseline: MAE={baseline_scores["MAE"]:.4f} RMSE={baseline_scores["RMSE"]:.4f} R2={baseline_scores["R2"]:.4f}')
    print(f'Outlier review flags in test: {int(flagged[test_idx].sum())} (NOT removed during training)')

    with (OUT / 'split_assignments.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['fingerprint_row', 'molecule_chembl_id', 'scaffold_smiles', 'split', 'pIC50', 'outlier_review_flag'])
        for split, inds in [('train', train_idx), ('validation', val_idx), ('test', test_idx)]:
            for i in inds:
                w.writerow([int(i), metadata[i]['molecule_chembl_id'], scaffolds[i], split, float(y[i]), bool(flagged[i])])
    with (OUT / 'test_predictions.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['fingerprint_row', 'molecule_chembl_id', 'actual_pIC50', 'predicted_pIC50', 'residual_actual_minus_predicted', 'outlier_review_flag'])
        for i, pred in zip(test_idx, preds['test']):
            w.writerow([int(i), metadata[i]['molecule_chembl_id'], float(y[i]), float(pred), float(y[i]-pred), bool(flagged[i])])

    low, high = min(y[test_idx].min(), preds['test'].min()), max(y[test_idx].max(), preds['test'].max())
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(y[test_idx], preds['test'], alpha=0.35, s=12)
    ax.plot([low, high], [low, high], linestyle='--', linewidth=1)
    ax.set(xlabel='Experimental median pIC50', ylabel='Predicted pIC50', title='EGFR: scaffold-holdout test predictions')
    fig.tight_layout(); fig.savefig(OUT / 'test_actual_vs_predicted.png', dpi=170); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.hist(y[test_idx] - preds['test'], bins=40)
    ax.set(xlabel='Residual (actual - predicted)', ylabel='Molecules', title='EGFR test residuals')
    fig.tight_layout(); fig.savefig(OUT / 'test_residuals.png', dpi=170); plt.close(fig)

    model_path = MODELS / 'EGFR_RF_scaffold_baseline.joblib'
    joblib.dump(model, model_path)
    report = {
        'target':'EGFR', 'target_variable':'stored molecule-level median pIC50 (provisional pending source verification)',
        'features':'RDKit Morgan ECFP4 radius=2, 2048 bits', 'split_method':'two-stage GroupShuffleSplit on Bemis-Murcko scaffolds',
        'random_seed':SEED, 'rdkit_version':rdBase.rdkitVersion, 'sklearn_version':sklearn_version,
        'model':'RandomForestRegressor(n_estimators=200, min_samples_leaf=2, max_features=sqrt)',
        'counts':{'all':len(y),'train':len(train_idx),'validation':len(val_idx),'test':len(test_idx),
                  'unique_scaffolds':len(set(scaffolds))},
        'scores': scores, 'test_training_median_baseline':baseline_scores,
        'test_excluding_extreme_flags_sensitivity_only':sensitivity,
        'test_extreme_review_flags':int(flagged[test_idx].sum()),
        'caveats':['Original input SMILES not standardized for salts/tautomers',
                   'Training and evaluation include extreme flags; review before final claims',
                   'Scaffold holdout tests generalization across scaffolds but cannot guarantee tryptanthrin-domain applicability',
                   'No tuning on test set; scaffold group split proportions may differ from nominal fractions']
    }
    (OUT / 'metrics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'Outputs: {OUT.relative_to(ROOT)}; {model_path.relative_to(ROOT)}')
    print('Step 3 complete. Test scores are preliminary, not final publication-quality performance.')

if __name__ == '__main__':
    main()
