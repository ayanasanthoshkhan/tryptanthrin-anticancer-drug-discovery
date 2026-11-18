"""Step 6: structure-based applicability screening for the tryptanthrin analog library.

Run from the existing capstone project root:
    .\\.venv\\Scripts\\python.exe scripts/06_applicability_domain.py

Uses the EXACT train-only molecule indices saved by Steps 3 and 5.
Does NOT use experimental activity labels, models, validation, or test compounds
as the reference for analog similarity. No files from earlier steps are changed.

Tanimoto thresholds below are *exploratory*, not statistically calibrated
applicability-domain guarantees. Fingerprints follow the prior steps (Morgan
radius=2, fpSize=2048) and use SMILES AS PROVIDED, including salts if present.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
ANALOG_PATH = ROOT / 'data' / 'analogs' / 'filtered_analog_library.csv'
PROCESSED = ROOT / 'data' / 'processed'
REPORTS = ROOT / 'reports'
OUT_DIR = REPORTS / 'applicability_domain'
TARGETS = ['EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2']
RADIUS = 2
FP_BITS = 2048
# Operational bins only; NOT empirically validated thresholds.
HIGH_SIM = 0.60
MODERATE_SIM = 0.40


def read_csv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'Required file not found: {path}')
    with path.open('r', newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def require_columns(rows, needed, source):
    if not rows:
        raise ValueError(f'Empty dataset: {source}')
    missing = set(needed) - set(rows[0])
    if missing:
        raise ValueError(f'{source} is missing columns: {sorted(missing)}')


def load_training_rows(target: str):
    meta_path = PROCESSED / f'{target}_fingerprint_metadata.csv'
    split_path = (REPORTS / 'egfr_random_forest' / 'split_assignments.csv'
                  if target == 'EGFR' else REPORTS / target.lower() / 'split_assignments.csv')
    meta = read_csv(meta_path)
    assignments = read_csv(split_path)
    require_columns(meta, ['fingerprint_row', 'molecule_chembl_id', 'canonical_smiles'], meta_path)
    require_columns(assignments, ['fingerprint_row', 'molecule_chembl_id', 'split'], split_path)
    if len(assignments) != len(meta):
        raise ValueError(f'{target}: split assignments and metadata differ in length')
    if any(int(row['fingerprint_row']) != i for i, row in enumerate(meta)):
        raise ValueError(f'{target}: metadata is not in fingerprint_row order')
    observed = set()
    training = []
    split_counts = Counter()
    for row in assignments:
        i = int(row['fingerprint_row'])
        if i < 0 or i >= len(meta) or i in observed:
            raise ValueError(f'{target}: duplicate or invalid fingerprint row index {i}')
        observed.add(i)
        if row['molecule_chembl_id'] != meta[i]['molecule_chembl_id']:
            raise ValueError(f'{target}: split metadata mismatch at row {i}')
        split = row['split']
        if split not in ('train', 'validation', 'test'):
            raise ValueError(f'{target}: unknown split {split!r}')
        split_counts[split] += 1
        if split == 'train':
            training.append(meta[i])
    if len(observed) != len(meta) or not all(split_counts[x] > 0 for x in ('train','validation','test')):
        raise ValueError(f'{target}: incomplete train/validation/test assignments')
    return training, dict(split_counts)


def load_analogs(generator):
    rows = read_csv(ANALOG_PATH)
    require_columns(rows, ['Analog_ID', 'SMILES', 'Decision'], ANALOG_PATH)
    ids = [r['Analog_ID'] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('Analog_ID values are not unique')
    mol_fps = []
    invalid = []
    for i, row in enumerate(rows):
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is None:
            invalid.append((i, row['Analog_ID']))
            mol_fps.append(None)
        else:
            mol_fps.append(generator.GetFingerprint(mol))
    if invalid:
        raise ValueError(f'{len(invalid)} invalid analog SMILES; examples: {invalid[:5]}')
    return rows, mol_fps


def similarity_band(score):
    if score >= HIGH_SIM:
        return 'HIGH_SIMILARITY'
    if score >= MODERATE_SIM:
        return 'MODERATE_SIMILARITY'
    return 'LOW_SIMILARITY_REVIEW'


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=FP_BITS)
    analogs, analog_fps = load_analogs(generator)
    print(f'Analogs: {len(analogs):,} | RDKit: {rdBase.rdkitVersion}', flush=True)
    print('Using TRAIN-only neighbors; no validation/test structures in the reference.', flush=True)

    output_csv = OUT_DIR / 'analog_target_similarity.csv'
    columns = ['Analog_ID','Target','Max_train_Tanimoto','Similarity_band',
               'Nearest_train_chembl_id','Exact_training_fingerprint_match',
               'Decision','SA_score','Category']
    summary = {}
    with output_csv.open('w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        for target in TARGETS:
            training, split_counts = load_training_rows(target)
            train_fps = []
            train_ids = []
            for r in training:
                mol = Chem.MolFromSmiles(r['canonical_smiles'])
                if mol is None:
                    raise ValueError(f'{target}: malformed training SMILES: {r["molecule_chembl_id"]}')
                train_fps.append(generator.GetFingerprint(mol))
                train_ids.append(r['molecule_chembl_id'])
            bands = Counter()
            total_sim = 0.0
            exact_matches = 0
            print(f'\n{target}: {len(train_fps):,} train molecules, '
                  f'{split_counts["validation"]:,} validation, '
                  f'{split_counts["test"]:,} test', flush=True)
            for i, (analog, fp) in enumerate(zip(analogs, analog_fps), start=1):
                # BulkTanimotoSimilarity runs C++ similarity against all TRAIN fingerprints.
                similarities = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
                nearest = max(range(len(similarities)), key=similarities.__getitem__)
                best = float(similarities[nearest])
                band = similarity_band(best)
                bands[band] += 1
                total_sim += best
                exact = best >= 1.0 - 1e-12
                exact_matches += int(exact)
                writer.writerow({
                    'Analog_ID': analog['Analog_ID'], 'Target': target,
                    'Max_train_Tanimoto': f'{best:.6f}',
                    'Similarity_band': band,
                    'Nearest_train_chembl_id': train_ids[nearest],
                    'Exact_training_fingerprint_match': exact,
                    'Decision': analog.get('Decision',''),
                    'SA_score': analog.get('SA_score',''),
                    'Category': analog.get('Category','')
                })
                if i % 2000 == 0 or i == len(analogs):
                    print(f'  processed {i:,}/{len(analogs):,} analogs', flush=True)
            summary[target] = {
                'training_molecules': len(train_fps),
                'validation_molecules_excluded_from_reference': split_counts['validation'],
                'test_molecules_excluded_from_reference': split_counts['test'],
                'analog_count': len(analogs),
                'mean_nearest_training_tanimoto': round(total_sim/len(analogs), 4),
                'exact_fingerprint_matches': exact_matches,
                'similarity_band_counts': dict(bands),
            }
            print(f'  mean nearest similarity: {total_sim/len(analogs):.3f} | '
                  f'low similarity: {bands["LOW_SIMILARITY_REVIEW"]:,} | '
                  f'exact matches: {exact_matches:,}', flush=True)

    explanation = {
        'scope': 'Molecular fingerprint nearest-neighbor similarity, NOT predictive uncertainty calibration',
        'fingerprint': f'RDKit Morgan radius {RADIUS}, {FP_BITS} bits (ECFP4-like)',
        'reference': 'TRAIN split only for each target (as saved by Step 3/5)',
        'thresholds_exploratory_only': {
            'high_similarity': f'>= {HIGH_SIM}',
            'moderate_similarity': f'>= {MODERATE_SIM} and < {HIGH_SIM}',
            'low_similarity_review': f'< {MODERATE_SIM}',
        },
        'warnings': [
            'Thresholds are heuristic and have not been calibrated against prediction errors.',
            'Similarity cannot establish that a model predicts the correct activity.',
            'Exact fingerprint match is not necessarily identical chemical structure (hash collisions or stereochemistry).',
            'Training/reference SMILES may include salts; no additional standardization performed.',
            'No pIC50 predictions are generated in this step.',
            'Analogs flagged by AD should not be automatically deleted.'
        ],
        'per_target': summary,
    }
    (OUT_DIR / 'similarity_summary.json').write_text(json.dumps(explanation, indent=2), encoding='utf-8')
    with (OUT_DIR / 'similarity_summary.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['Target','Training_molecules','Analogs','Mean_nearest_Tanimoto',
                    'High_similarity','Moderate_similarity','Low_similarity_review','Exact_fingerprint_matches'])
        for target, s in summary.items():
            c = s['similarity_band_counts']
            w.writerow([target,s['training_molecules'],s['analog_count'],
                        s['mean_nearest_training_tanimoto'],c.get('HIGH_SIMILARITY',0),
                        c.get('MODERATE_SIMILARITY',0),c.get('LOW_SIMILARITY_REVIEW',0),
                        s['exact_fingerprint_matches']])
    print('\nStep 6 complete.', flush=True)
    print(f'Outputs: {output_csv.relative_to(ROOT)} and similarity_summary.csv/.json', flush=True)
    print('Similarity bands are exploratory; do not interpret as validated model confidence.', flush=True)


if __name__ == '__main__':
    main()
