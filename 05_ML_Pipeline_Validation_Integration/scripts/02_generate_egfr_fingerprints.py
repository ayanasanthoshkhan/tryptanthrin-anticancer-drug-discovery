"""Step 2: create reproducible EGFR ECFP4 fingerprints, without model training.

Inputs: data/raw/EGFR_Final_ML_Dataset.csv
Outputs: data/processed/EGFR_ECFP4_2048.npz,
         data/processed/EGFR_fingerprint_metadata.csv,
         reports/egfr_fingerprint_summary.txt

Input rows and stored pIC50 are retained unchanged. No split or normalization.
"""
from pathlib import Path
import csv
import math
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data' / 'raw' / 'EGFR_Final_ML_Dataset.csv'
OUTDIR = ROOT / 'data' / 'processed'
REPORTS = ROOT / 'reports'
BITS = 2048
RADIUS = 2


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=BITS)
    x_rows, labels, records = [], [], []
    failed = []
    with SOURCE.open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for source_row, row in enumerate(reader, start=2):
            smiles = (row.get('canonical_smiles') or '').strip()
            compound = (row.get('molecule_chembl_id') or '').strip()
            try:
                y = float(row['pIC50'])
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                if mol is None or not math.isfinite(y):
                    raise ValueError('invalid SMILES or non-finite pIC50')
                fp = generator.GetFingerprintAsNumPy(mol).astype(np.uint8)
            except (ValueError, TypeError) as exc:
                failed.append((source_row, compound, str(exc)))
                continue
            x_rows.append(fp)
            labels.append(y)
            records.append({
                'fingerprint_row': len(x_rows) - 1,
                'source_csv_row': source_row,
                'molecule_chembl_id': compound,
                'canonical_smiles': smiles,
                'pIC50': row['pIC50'],
                'measurement_count': row.get('measurement_count', ''),
                'review_extreme_pIC50': str(y < 3 or y > 12),
            })
    if failed:
        print(f'WARNING: {len(failed)} rows failed; no files written. First errors: {failed[:5]}')
        raise SystemExit(1)
    x = np.stack(x_rows, axis=0)
    y = np.asarray(labels, dtype=np.float64)
    assert x.shape == (len(records), BITS) and len(y) == len(records)
    assert np.isin(x, [0, 1]).all()
    fp_file = OUTDIR / 'EGFR_ECFP4_2048.npz'
    np.savez_compressed(fp_file, X=x, y=y)
    meta_file = OUTDIR / 'EGFR_fingerprint_metadata.csv'
    with meta_file.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    report = [
        'EGFR ECFP4 / Morgan fingerprint generation',
        f'Input: {SOURCE.relative_to(ROOT)}',
        f'RDKit version: {rdBase.rdkitVersion}',
        f'Molecules: {len(records):,}',
        f'X shape: {x.shape}',
        f'y shape: {y.shape}',
        f'Radius: {RADIUS}, bits: {BITS}',
        f'Unique molecular IDs: {len({r["molecule_chembl_id"] for r in records}):,}',
        f'Extreme pIC50 review flags: {sum(r["review_extreme_pIC50"] == "True" for r in records)}',
        f'Fraction of set bits: {x.mean():.4f}',
        'Stored median pIC50 labels copied unchanged; this is not a train/test split.',
        'Multi-component SMILES (e.g. salts) are fingerprinted as received; standardization remains a later decision.',
        'No machine-learning model is trained at this stage.',
        f'Outputs: {fp_file.relative_to(ROOT)}, {meta_file.relative_to(ROOT)}',
    ]
    (REPORTS / 'egfr_fingerprint_summary.txt').write_text('\n'.join(report) + '\n', encoding='utf-8')
    print('\n'.join(report))
    print('Step 2 complete.')

if __name__ == '__main__':
    main()
