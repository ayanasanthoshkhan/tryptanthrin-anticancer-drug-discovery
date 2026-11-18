"""Step 7: predict filtered tryptanthrin analogs for five targets.

Run from anywhere:
    python scripts/07_predict_analog_activity.py

Uses BOTH pre-trained RF and XGBoost baselines, preserving earlier artifacts.
Predictions are computational screening estimates, not observed activity or
validated uncertainty intervals. AD similarity cutoffs are exploratory.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
ANALOGS = ROOT / 'data' / 'analogs' / 'filtered_analog_library.csv'
AD_FILE = ROOT / 'reports' / 'applicability_domain' / 'analog_target_similarity.csv'
PROCESSED = ROOT / 'data' / 'processed'
MODELS = ROOT / 'models'
OUT = ROOT / 'reports' / 'analog_predictions'
TARGETS = ['EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2']
BITS = 2048
BATCH = 512


def read_rows(path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing file: {path}')
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def require(rows, fields, source):
    if not rows:
        raise ValueError(f'Empty CSV: {source}')
    missing = set(fields) - set(rows[0])
    if missing:
        raise ValueError(f'Missing columns in {source}: {sorted(missing)}')


def make_analog_features(rows):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=BITS)
    X = np.empty((len(rows), BITS), dtype=np.uint8)
    for i, row in enumerate(rows):
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is None:
            raise ValueError(f'Invalid SMILES: {row["Analog_ID"]}')
        X[i] = np.asarray(generator.GetFingerprintAsNumPy(mol), dtype=np.uint8)
    return X


def predict_batched(model, X):
    result = np.empty(len(X), dtype=np.float64)
    for start in range(0, len(X), BATCH):
        end = min(start + BATCH, len(X))
        result[start:end] = model.predict(X[start:end])
    if not np.isfinite(result).all():
        raise ValueError('Non-finite prediction found')
    return result


def get_training_range(target):
    feature_path = PROCESSED / f'{target}_ECFP4_2048.npz'
    split_path = (ROOT / 'reports' / 'egfr_random_forest' / 'split_assignments.csv'
                  if target == 'EGFR' else ROOT / 'reports' / target.lower() / 'split_assignments.csv')
    metadata_path = PROCESSED / f'{target}_fingerprint_metadata.csv'
    split_rows = read_rows(split_path)
    metadata = read_rows(metadata_path)
    require(split_rows, ['fingerprint_row', 'molecule_chembl_id', 'split'], split_path)
    require(metadata, ['fingerprint_row', 'molecule_chembl_id'], metadata_path)
    with np.load(feature_path) as archive:
        y = archive['y']
    if len(y) != len(split_rows) or len(y) != len(metadata):
        raise ValueError(f'{target}: inconsistent feature/metadata/split sizes')
    seen = set()
    train_indices = []
    for row in split_rows:
        i = int(row['fingerprint_row'])
        if i in seen or i >= len(y) or i < 0 or row['molecule_chembl_id'] != metadata[i]['molecule_chembl_id']:
            raise ValueError(f'{target}: invalid split row {i}')
        seen.add(i)
        if row['split'] == 'train':
            train_indices.append(i)
        elif row['split'] not in {'validation', 'test'}:
            raise ValueError(f'{target}: unknown split')
    if len(seen) != len(y) or not train_indices:
        raise ValueError(f'{target}: missing split rows')
    ytrain = y[train_indices]
    return float(np.min(ytrain)), float(np.max(ytrain)), len(train_indices)


def main():
    analogs = read_rows(ANALOGS)
    require(analogs, ['Analog_ID', 'SMILES', 'Decision', 'SA_score', 'QED'], ANALOGS)
    analog_ids = [r['Analog_ID'] for r in analogs]
    if len(set(analog_ids)) != len(analog_ids):
        raise ValueError('Duplicate Analog_ID in analog library')
    ad = read_rows(AD_FILE)
    require(ad, ['Analog_ID', 'Target', 'Max_train_Tanimoto', 'Similarity_band', 'Nearest_train_chembl_id'], AD_FILE)
    ad_lookup = {}
    for row in ad:
        key = (row['Target'], row['Analog_ID'])
        if key in ad_lookup:
            raise ValueError(f'Duplicate AD record: {key}')
        ad_lookup[key] = row
    if len(ad_lookup) != len(analogs) * len(TARGETS):
        raise ValueError(f'Unexpected AD row count: {len(ad_lookup)}')
    if set(ad_lookup) != {(target, aid) for target in TARGETS for aid in analog_ids}:
        raise ValueError('AD IDs or target names disagree with analog library')

    print(f'Step 7 | analogs={len(analogs):,} | RDKit={rdBase.rdkitVersion}', flush=True)
    print('Generating analog Morgan ECFP4 fingerprints once...', flush=True)
    X = make_analog_features(analogs)
    OUT.mkdir(parents=True, exist_ok=True)
    columns = ['Analog_ID','Target','SMILES','Category','Decision','QED','SA_score',
               'Max_train_Tanimoto','Similarity_band','Nearest_train_chembl_id',
               'Pred_RF_pIC50','Pred_XGB_pIC50','Model_disagreement_abs',
               'RF_outside_training_pIC50_range','XGB_outside_training_pIC50_range',
               'Priority_tier']
    per_target = []
    target_file = OUT / 'all_target_predictions.csv'
    with target_file.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for target in TARGETS:
            paths = [MODELS / f'{target}_{name}_scaffold_baseline.joblib' for name in ('RF','XGB')]
            for path in paths:
                if not path.is_file():
                    raise FileNotFoundError(f'Missing trained model: {path}')
            rf = joblib.load(paths[0])
            xgb = joblib.load(paths[1])
            print(f'\n{target}: predicting RF and XGBoost...', flush=True)
            rf_pred = predict_batched(rf, X)
            xgb_pred = predict_batched(xgb, X)
            train_min, train_max, train_size = get_training_range(target)
            counts = Counter()
            for i, analog in enumerate(analogs):
                similarity = ad_lookup[target, analog['Analog_ID']]
                band = similarity['Similarity_band']
                # This is a structural-review category, NOT prediction confidence.
                if band == 'HIGH_SIMILARITY':
                    tier = 'HIGH_SIMILARITY_SCREEN'
                elif band == 'MODERATE_SIMILARITY':
                    tier = 'MODERATE_SIMILARITY_SCREEN'
                elif band == 'LOW_SIMILARITY_REVIEW':
                    tier = 'LOW_SIMILARITY_EXPLORATORY_ONLY'
                else:
                    raise ValueError(f'Unexpected similarity band: {band}')
                counts[tier] += 1
                writer.writerow({
                    'Analog_ID': analog['Analog_ID'], 'Target': target,
                    'SMILES': analog['SMILES'], 'Category': analog.get('Category',''),
                    'Decision': analog['Decision'], 'QED': analog.get('QED',''),
                    'SA_score': analog.get('SA_score',''),
                    'Max_train_Tanimoto': similarity['Max_train_Tanimoto'],
                    'Similarity_band': band,
                    'Nearest_train_chembl_id': similarity['Nearest_train_chembl_id'],
                    'Pred_RF_pIC50': f'{rf_pred[i]:.5f}',
                    'Pred_XGB_pIC50': f'{xgb_pred[i]:.5f}',
                    'Model_disagreement_abs': f'{abs(rf_pred[i]-xgb_pred[i]):.5f}',
                    'RF_outside_training_pIC50_range': bool(rf_pred[i] < train_min or rf_pred[i] > train_max),
                    'XGB_outside_training_pIC50_range': bool(xgb_pred[i] < train_min or xgb_pred[i] > train_max),
                    'Priority_tier': tier,
                })
            row = {'Target': target, 'Analogs': len(analogs),
                   'Train_reference_molecules': train_size,
                   'Train_pIC50_min': round(train_min, 4),
                   'Train_pIC50_max': round(train_max, 4),
                   **dict(counts),
                   'Mean_RF_prediction': round(float(rf_pred.mean()), 4),
                   'Mean_XGB_prediction': round(float(xgb_pred.mean()), 4),
                   'Mean_abs_model_disagreement': round(float(np.mean(np.abs(rf_pred-xgb_pred))), 4)}
            per_target.append(row)
            print('  tiers:', dict(counts), '| mean absolute RF–XGB difference:',
                  f'{row["Mean_abs_model_disagreement"]:.3f}', flush=True)

    with (OUT / 'prediction_summary.csv').open('w', newline='', encoding='utf-8') as f:
        columns_summary = ['Target','Analogs','Train_reference_molecules','Train_pIC50_min','Train_pIC50_max',
                           'HIGH_SIMILARITY_SCREEN','MODERATE_SIMILARITY_SCREEN','LOW_SIMILARITY_EXPLORATORY_ONLY',
                           'Mean_RF_prediction','Mean_XGB_prediction','Mean_abs_model_disagreement']
        w = csv.DictWriter(f, fieldnames=columns_summary)
        w.writeheader()
        w.writerows(per_target)
    (OUT / 'methodology.json').write_text(json.dumps({
        'models': 'Existing trained Step 3–5 RF and XGBoost models; NO new training',
        'fingerprint': 'RDKit Morgan radius 2, 2048 bits, SMILES as provided',
        'source': str(ANALOGS.relative_to(ROOT)),
        'applicability': 'Step 6 nearest TRAIN-only Tanimoto similarity',
        'labels': 'Stored median individual-measurement pIC50 values',
        'warning': ['No experimental activities were measured for generated analogs by this script.',
                    'Low similarity indicates extrapolation; do NOT promote as a validated hit.',
                    'Similarity tiers are heuristic and NOT calibrated uncertainty.',
                    'The models are trained across broad chemical space; predictions are not target-engagement validation.',
                    'Both model predictions are reported; no model chosen using test results.',
                    'Training label-range flags do not capture model extrapolation fully.',
                    'Further chemical standardization, analogue overlap checks, and evaluation of stability are needed.']
    }, indent=2), encoding='utf-8')
    print(f'\nStep 7 complete. Wrote {len(analogs)*len(TARGETS):,} prediction rows.', flush=True)
    print(f'Outputs: {target_file.relative_to(ROOT)} and prediction_summary.csv', flush=True)
    print('Do NOT interpret low-similarity predictions as confirmed anticancer hits.', flush=True)

if __name__ == '__main__':
    main()
