r"""Step 12 — fingerprint-overlap sensitivity of FIXED held-out test predictions.

Windows PowerShell (from project root):
    .\.venv\Scripts\python.exe scripts/12_fingerprint_collision_sensitivity.py
Optional target-specific analysis:
    .\.venv\Scripts\python.exe scripts/12_fingerprint_collision_sensitivity.py --targets EGFR TDO2

Read-only with respect to all pre-existing project files. This analysis DOES NOT
retrain, select, or tune models. It measures how many held-out test MOLECULES have
an identical input ECFP4 vector anywhere in the training partition and recomputes
metrics with those test examples excluded. Exclusion is a post-hoc diagnostic,
NOT a corrected unbiased test estimate. Outputs only under reports/final_validation/.
"""
from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
OUT = ROOT / 'reports' / 'final_validation'


def load_csv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing input: {path}. Run the corresponding prior steps first.')
    with path.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def save_csv(filename, rows, columns):
    with (OUT / filename).open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def evaluate(y, pred):
    n = len(y)
    if n == 0:
        return dict(n=0, MAE=None, RMSE=None, R2=None)
    return dict(n=n, MAE=float(mean_absolute_error(y, pred)),
                RMSE=float(np.sqrt(mean_squared_error(y, pred))),
                R2=float(r2_score(y, pred)) if n >= 2 and np.ptp(y) > 0 else None)


def compact(v):
    return 'NA' if v is None else f'{v:.4f}'


def path_for_predictions(target, model):
    if target == 'EGFR':
        folder = 'egfr_random_forest' if model == 'RF' else 'egfr_xgboost'
    else:
        folder = f'{target.lower()}/{"random_forest" if model == "RF" else "xgboost"}'
    return ROOT / 'reports' / folder / 'test_predictions.csv'


def inspect_target(target):
    base = ROOT / 'data' / 'processed' / f'{target}_ECFP4_2048.npz'
    if not base.is_file():
        raise FileNotFoundError(f'Missing fingerprint file: {base}')
    with np.load(base, allow_pickle=False) as npz:
        X = npz['X']
        y = npz['y'].astype(np.float64)
    if X.shape != (len(y), 2048) or not np.all(np.isin(X, [0, 1])):
        raise ValueError(f'{target}: unexpected fingerprint matrix shape or values: {X.shape}')
    if not np.isfinite(y).all():
        raise ValueError(f'{target}: nonfinite target labels')

    split_folder = 'egfr_random_forest' if target == 'EGFR' else target.lower()
    assignments = load_csv(ROOT / 'reports' / split_folder / 'split_assignments.csv')
    if len(assignments) != len(y):
        raise ValueError(f'{target}: split assignment length does not match fingerprint rows')
    indexed = {}
    for rec in assignments:
        idx = int(rec['fingerprint_row'])
        if idx in indexed or not 0 <= idx < len(y) or rec['split'] not in ('train','validation','test'):
            raise ValueError(f'{target}: invalid/duplicate split index {idx}')
        indexed[idx] = rec
        if abs(float(rec['pIC50']) - y[idx]) > 1e-7:
            raise ValueError(f'{target}: split label mismatch at {idx}')
    if set(indexed) != set(range(len(y))):
        raise ValueError(f'{target}: missing split indices')

    # Pack to 256-byte keys (2048 bits) to avoid large Python byte objects and
    # to ensure comparison is exact on the input features actually seen by ML.
    packed = np.packbits(X.astype(np.uint8), axis=1)
    del X
    train = sorted(i for i in indexed if indexed[i]['split'] == 'train')
    validation = sorted(i for i in indexed if indexed[i]['split'] == 'validation')
    test = sorted(i for i in indexed if indexed[i]['split'] == 'test')
    if not train or not validation or not test:
        raise ValueError(f'{target}: empty partition')
    train_representatives = defaultdict(list)
    for i in train:
        train_representatives[packed[i].tobytes()].append(i)

    matching_test = {}
    for i in test:
        matches = train_representatives.get(packed[i].tobytes())
        if matches:
            matching_test[i] = matches
    collision_set = set(matching_test)
    collision_ids = set(packed[i].tobytes() for i in collision_set)
    total_unique_overlaps = len(collision_ids)
    rows_details = []
    for i in sorted(matching_test):
        matches = matching_test[i]
        # First matching training molecule is a traceable representative, NOT
        # necessarily the closest activity label or the only matched molecule.
        rep = matches[0]
        rows_details.append({
            'target':target, 'test_fingerprint_row':i,
            'test_molecule_id':indexed[i]['molecule_chembl_id'],
            'test_pIC50':float(y[i]),
            'training_match_count':len(matches),
            'representative_training_fingerprint_row':rep,
            'representative_training_molecule_id':indexed[rep]['molecule_chembl_id'],
            'representative_training_pIC50':float(y[rep]),
            'abs_test_minus_representative_train_pIC50':abs(float(y[i]-y[rep])),
            'mean_training_match_pIC50':float(np.mean(y[matches])),
            'max_abs_test_minus_training_match_pIC50':float(np.max(np.abs(y[matches]-y[i]))),
        })
    output = []
    for model in ('RF','XGB'):
        path = path_for_predictions(target, model)
        predictions = load_csv(path)
        pred_by_idx = {}
        expected = set(test)
        for rec in predictions:
            i = int(rec['fingerprint_row'])
            if i in pred_by_idx or i not in expected:
                raise ValueError(f'{target} {model}: unexpected duplicate/non-test prediction index {i}')
            if rec['molecule_chembl_id'] != indexed[i]['molecule_chembl_id']:
                raise ValueError(f'{target} {model}: prediction molecule-ID mismatch at {i}')
            if not np.isclose(float(rec['actual_pIC50']), y[i], rtol=0, atol=1e-7):
                raise ValueError(f'{target} {model}: prediction actual-pIC50 mismatch at {i}')
            v = float(rec['predicted_pIC50'])
            if not np.isfinite(v):
                raise ValueError(f'{target} {model}: nonfinite prediction')
            pred_by_idx[i] = v
        if set(pred_by_idx) != expected:
            raise ValueError(f'{target} {model}: missing test predictions')
        for subset, ids in [('full_test', test),
                            ('no_train_identical_fingerprint', [i for i in test if i not in collision_set]),
                            ('identical_fingerprint_only', [i for i in test if i in collision_set])]:
            y_true = y[ids]
            y_pred = np.asarray([pred_by_idx[i] for i in ids], dtype=float)
            output.append({'target':target,'model':model,'subset':subset,
                           'test_total':len(test),
                           'test_molecules_matching_train_fingerprint':len(collision_set),
                           'test_matching_fraction':len(collision_set)/len(test),
                           'distinct_shared_train_test_fingerprints':total_unique_overlaps,
                           **evaluate(y_true, y_pred)})
    print(f'{target:<10} test={len(test):>5,} | affected test molecules={len(collision_set):>4,} '
          f'({len(collision_set)/len(test):.2%}) | distinct shared fingerprints={total_unique_overlaps}')
    for model in ('RF','XGB'):
        r = {a['subset']:a for a in output if a['model'] == model}
        print(f'  {model:<3} R2 full={compact(r["full_test"]["R2"])}  '
              f'without-overlap={compact(r["no_train_identical_fingerprint"]["R2"])} | '
              f'MAE full={compact(r["full_test"]["MAE"])}  '
              f'without-overlap={compact(r["no_train_identical_fingerprint"]["MAE"])}')
    return output, rows_details


def main():
    parser = argparse.ArgumentParser(description='Fixed-test fingerprint collision sensitivity; no retraining.')
    parser.add_argument('--targets',nargs='+',choices=TARGETS,default=list(TARGETS))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    summary, details = [], []
    print('Step 12 — test fingerprint-collision sensitivity (NO retraining)\n', flush=True)
    for target in args.targets:
        metrics, affected = inspect_target(target)
        summary.extend(metrics)
        details.extend(affected)
    save_csv('fingerprint_collision_sensitivity_metrics.csv', summary,
             ['target','model','subset','test_total','test_molecules_matching_train_fingerprint',
              'test_matching_fraction','distinct_shared_train_test_fingerprints','n','MAE','RMSE','R2'])
    save_csv('fingerprint_collision_affected_test_molecules.csv', details,
             ['target','test_fingerprint_row','test_molecule_id','test_pIC50',
              'training_match_count','representative_training_fingerprint_row',
              'representative_training_molecule_id','representative_training_pIC50',
              'abs_test_minus_representative_train_pIC50','mean_training_match_pIC50',
              'max_abs_test_minus_training_match_pIC50'])
    lines = [
        '# Step 12 — Fingerprint Collision Sensitivity', '',
        'This is a read-only, post-hoc diagnostic using existing fixed-test predictions.',
        'No models were retrained, no predictions were recalculated, and no source data were modified.', '',
        '## Results', '',
        '| Target | Model | Affected test molecules | Original test R2 | Without exact train-FP overlap R2 | Original MAE | Without overlap MAE |',
        '|---|---|---:|---:|---:|---:|---:|'
    ]
    for r in summary:
        if r['subset'] != 'full_test':
            continue
        clean = next(x for x in summary if x['target']==r['target'] and x['model']==r['model']
                     and x['subset']=='no_train_identical_fingerprint')
        lines.append(f'| {r["target"]} | {r["model"]} | {r["test_molecules_matching_train_fingerprint"]}/{r["test_total"]} '
                     f'| {compact(r["R2"])} | {compact(clean["R2"])} | {compact(r["MAE"])} | {compact(clean["MAE"])} |')
    lines.extend(['','## Interpretation','',
       '- Counts of affected **test molecules** can exceed counts of distinct overlapping fingerprint patterns.',
       '- Equal 2048-bit ECFP4 input vectors can occur for distinct structures; this is not automatically a canonical-structure duplicate.',
       '- The collision-excluded test subset has a different composition and may have different activity variance. A changed R2 does **not** by itself quantify bias caused by overlap.',
       '- The collision-only subset may be very small. Treat its scores as descriptive, not robust.',
       '- The original fixed test metrics remain the primary reported baseline, accompanied by this sensitivity check and scaffold-disjointness limitations.',
       '- This analysis does not assess stereochemistry-aware fingerprints, tautomer/salt standardization, assay heterogeneity, calibration, or performance on externally measured triptanthrin analogs.',
       '- Avoid selecting or tuning models against these post-hoc test metrics.', ''])
    (OUT/'fingerprint_collision_sensitivity_report.md').write_text('\n'.join(lines),encoding='utf-8')
    print('\nStep 12 complete. Saved under reports\\final_validation\\:')
    print('  fingerprint_collision_sensitivity_metrics.csv')
    print('  fingerprint_collision_affected_test_molecules.csv')
    print('  fingerprint_collision_sensitivity_report.md')
    print('No training, test, model, or shortlist files changed.')

if __name__=='__main__':
    main()
