r"""Step 11 — non-destructive final ML validation audit, with optional split stability.

Run from project root (Windows PowerShell):
    .\.venv\Scripts\python.exe scripts/11_final_ml_validation_audit.py
Optional exploratory RF repeat-scaffold-split check (retrains NEW temporary RFs):
    .\.venv\Scripts\python.exe scripts/11_final_ml_validation_audit.py --repeat-splits 3

Reads original data, metadata, split assignments, test predictions; NEVER rewrites
existing models, splits, predictions or source data. Saves under reports/final_validation/.
Repeated-split RF metrics are exploratory stability estimates, not new model
selection results, and do not re-evaluate the original fixed test as a tuning loop.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
OUT = ROOT / 'reports' / 'final_validation'
SPLITS = ('train', 'validation', 'test')
FP_SIZE = 2048


def read_csv(path):
    if not path.is_file():
        raise FileNotFoundError(f'Required file missing: {path}\nRun Steps 2–5 first; do not modify input filenames.')
    with path.open('r', newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def save_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def scores(y, pred):
    if len(y) < 2:
        return {'n': len(y), 'MAE': None, 'RMSE': None, 'R2': None}
    return {'n': len(y), 'MAE': float(mean_absolute_error(y, pred)),
            'RMSE': float(np.sqrt(mean_squared_error(y, pred))),
            'R2': float(r2_score(y, pred)) if len(set(y)) > 1 else None}


def brief(n):
    return f'{n:,}'


def read_target(t):
    raw = read_csv(ROOT / 'data' / 'raw' / f'{t}_Final_ML_Dataset.csv')
    metadata = read_csv(ROOT / 'data' / 'processed' / f'{t}_fingerprint_metadata.csv')
    split_dir = 'egfr_random_forest' if t == 'EGFR' else t.lower()
    splits = read_csv(ROOT / 'reports' / split_dir / 'split_assignments.csv')
    with np.load(ROOT / 'data' / 'processed' / f'{t}_ECFP4_2048.npz') as data:
        X = data['X'].copy()
        y = data['y'].astype(np.float64).copy()
    if X.shape != (len(raw), FP_SIZE) or len(y) != len(raw) or len(metadata) != len(raw):
        raise ValueError(f'{t}: length/fingerprint mismatch; raw={len(raw)}, metadata={len(metadata)}, X={X.shape}')
    if not np.isfinite(y).all():
        raise ValueError(f'{t}: nonfinite activity labels')
    by_idx = {}
    for row in splits:
        idx = int(row['fingerprint_row'])
        if idx in by_idx:
            raise ValueError(f'{t}: duplicate split entry for index {idx}')
        if idx < 0 or idx >= len(raw) or row['split'] not in SPLITS:
            raise ValueError(f'{t}: invalid split assignment {row}')
        by_idx[idx] = row
    if set(by_idx) != set(range(len(raw))):
        raise ValueError(f'{t}: incomplete split assignments; missing {len(raw) - len(by_idx)}')
    for i, (r, m) in enumerate(zip(raw, metadata)):
        if int(m['fingerprint_row']) != i or by_idx[i]['molecule_chembl_id'] != r['molecule_chembl_id'] or m['molecule_chembl_id'] != r['molecule_chembl_id']:
            raise ValueError(f'{t}: split/metadata/raw molecule ID mismatch at row {i}')
        if not np.isclose(y[i], float(r['pIC50']), rtol=0, atol=1e-8):
            raise ValueError(f'{t}: fingerprint labels differ from raw at row {i}')
    return raw, metadata, by_idx, X, y


def bin_label(v):
    if v < 5:
        return '<5'
    if v < 6:
        return '5–<6'
    if v < 7:
        return '6–<7'
    if v < 8:
        return '7–<8'
    if v < 9:
        return '8–<9'
    return '>=9'


def check_target(t):
    raw, metadata, by_idx, X, y = read_target(t)
    idxs = {name: np.asarray([i for i in range(len(y)) if by_idx[i]['split'] == name], dtype=int) for name in SPLITS}
    if any(len(x) < 2 for x in idxs.values()):
        raise ValueError(f'{t}: an empty or singleton partition')
    # Canonicalize all molecules for exact chemical-structure overlap audit.
    canonical = []
    scaffold_recomputed = []
    for i, row in enumerate(raw):
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is None:
            raise ValueError(f'{t}: invalid SMILES in raw input row {i}')
        canonical.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold = Chem.MolToSmiles(core) or '__NO_MURCKO_SCAFFOLD__'
        scaffold_recomputed.append(scaffold)
        if scaffold != by_idx[i]['scaffold_smiles']:
            raise ValueError(f'{t}: scaffold mismatch at {i}; perhaps input changed since training')
    def overlaps(items):
        sets = {p: set(items[i] for i in ids) for p, ids in idxs.items()}
        return {f'{a}_vs_{b}': len(sets[a] & sets[b]) for a,b in
                [('train','validation'), ('train','test'), ('validation','test')]}
    scaffold_overlap = overlaps(scaffold_recomputed)
    exact_overlap = overlaps(canonical)
    # Count equal fingerprint BIT ARRAYS across partitions as additional representation leakage,
    # recognizing hash collisions / non-unique fingerprints are possible.
    fp_bytes = [X[i].tobytes() for i in range(len(y))]
    fingerprint_overlap = overlaps(fp_bytes)
    status = ('FAIL_STRUCTURE_OR_SCAFFOLD' if any(scaffold_overlap.values()) or any(exact_overlap.values())
              else 'REVIEW_FINGERPRINT_COLLISIONS' if any(fingerprint_overlap.values())
              else 'PASS')
    # Distinct molecules may generate the same folded ECFP4 bit vector.
    # Save a small, traceable set of across-split examples for scientific review.
    representatives = {}
    for part, ids in idxs.items():
        for i in ids:
            representatives.setdefault((part, fp_bytes[i]), int(i))
    fp_examples = []
    for a, b in [('train', 'validation'), ('train', 'test'), ('validation', 'test')]:
        left = {fp: i for (part, fp), i in representatives.items() if part == a}
        right = {fp: i for (part, fp), i in representatives.items() if part == b}
        for fp in left.keys() & right.keys():
            i, j = left[fp], right[fp]
            fp_examples.append({'target':t, 'split_a':a, 'split_b':b,
                'index_a':i, 'index_b':j,
                'id_a':raw[i]['molecule_chembl_id'], 'id_b':raw[j]['molecule_chembl_id'],
                'canonical_smiles_a':canonical[i], 'canonical_smiles_b':canonical[j],
                'exact_structure_match':canonical[i] == canonical[j],
                'scaffold_match':scaffold_recomputed[i] == scaffold_recomputed[j]})
    # Check molecule ID duplication as a distinct data integrity issue.
    id_counts = Counter(r['molecule_chembl_id'] for r in raw)
    split_counts = {p: len(ids) for p, ids in idxs.items()}
    train_y = y[idxs['train']]
    test_y = y[idxs['test']]
    train_lo, train_hi = float(np.min(train_y)), float(np.max(train_y))
    outside = int(((test_y < train_lo) | (test_y > train_hi)).sum())
    summary = {
        'target':t,'samples':len(y),'scaffolds':len(set(scaffold_recomputed)),
        'train':split_counts['train'],'validation':split_counts['validation'],'test':split_counts['test'],
        'split_leakage_status':status,
        'overlapping_scaffolds_train_test':scaffold_overlap['train_vs_test'],
        'overlapping_canonical_smiles_train_test':exact_overlap['train_vs_test'],
        'overlapping_fingerprints_train_test':fingerprint_overlap['train_vs_test'],
        'max_scaffold_overlap_any_pair':max(scaffold_overlap.values()),
        'max_exact_structure_overlap_any_pair':max(exact_overlap.values()),
        'max_fingerprint_overlap_any_pair':max(fingerprint_overlap.values()),
        'duplicated_molecule_ids':sum(n-1 for n in id_counts.values()),
        'test_outside_training_activity_range':outside,
        'test_extreme_review_values':int(((test_y<3)|(test_y>12)).sum()),
        'train_min_pIC50':train_lo,'train_max_pIC50':train_hi,
        'train_median_pIC50':float(np.median(train_y)),
        'test_median_pIC50':float(np.median(test_y))
    }
    bins = []
    for part, ids in idxs.items():
        counts = Counter(bin_label(float(v)) for v in y[ids])
        for bin_name in ('<5','5–<6','6–<7','7–<8','8–<9','>=9'):
            bins.append({'target':t, 'split':part, 'pIC50_bin':bin_name,
                         'count':counts[bin_name], 'fraction':counts[bin_name]/len(ids)})
    sensitivity = []
    paths = [
        ('Random Forest', ROOT/'reports'/('egfr_random_forest' if t=='EGFR' else t.lower()+'/random_forest')/'test_predictions.csv'),
        ('XGBoost', ROOT/'reports'/('egfr_xgboost' if t=='EGFR' else t.lower()+'/xgboost')/'test_predictions.csv')]
    for model_name, path in paths:
        preds = read_csv(path)
        test_set = set(map(int, idxs['test']))
        mapped = {}
        for row in preds:
            idx = int(row['fingerprint_row'])
            if idx in mapped:
                raise ValueError(f'{t} {model_name}: repeated test prediction index {idx}')
            mapped[idx] = float(row['predicted_pIC50'])
            if idx not in test_set or row['molecule_chembl_id'] != raw[idx]['molecule_chembl_id']:
                raise ValueError(f'{t} {model_name}: prediction split/ID mismatch at {idx}')
            if not np.isclose(float(row['actual_pIC50']), y[idx], atol=1e-7, rtol=0):
                raise ValueError(f'{t} {model_name}: prediction labels differ at {idx}')
        if set(mapped) != test_set:
            raise ValueError(f'{t} {model_name}: missing test predictions')
        ordered = idxs['test']
        true = y[ordered]
        guess = np.asarray([mapped[int(i)] for i in ordered])
        extreme = (true < 3) | (true > 12)
        groups = {
            'full_test': np.ones(len(ordered), dtype=bool),
            'without_extreme_review_flags': ~extreme,
            'below_train_pIC50_median': true < float(np.median(train_y)),
            'at_or_above_train_pIC50_median': true >= float(np.median(train_y))
        }
        for group, mask in groups.items():
            sc = scores(true[mask], guess[mask])
            sensitivity.append({'target':t,'model':model_name,'subgroup':group, **sc})
    # Compact activity distribution plot (aggregated by bin/split).
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for part in SPLITS:
        entries = [b for b in bins if b['split']==part]
        ax.plot([b['pIC50_bin'] for b in entries], [b['fraction'] for b in entries], marker='o',label=part)
    ax.set(xlabel='pIC50 interval',ylabel='Fraction within split',title=f'{t}: activity distribution by scaffold split')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/f'{t}_activity_distribution.png', dpi=160)
    plt.close(fig)
    return summary,bins,sensitivity,fp_examples,scaffold_recomputed,X,y


def repeat_rf_stability(t, X, y, scaffolds, count):
    # Repeat NEW scaffold splits; do not alter existing baseline model or split.
    # This is a Monte Carlo stability check, not a completely independent external validation.
    rows=[]
    s=np.asarray(scaffolds)
    for seed in range(101,101+count):
        train,hold=next(GroupShuffleSplit(n_splits=1,test_size=0.30,random_state=seed).split(X,y,groups=s))
        vrel,trel=next(GroupShuffleSplit(n_splits=1,test_size=0.50,random_state=seed).split(X[hold],y[hold],groups=s[hold]))
        test=hold[trel]
        if len(set(s[train]) & set(s[test])):
            raise ValueError(f'{t}: scaffold overlap in repeated split seed {seed}')
        # Same RF configuration as Step 5 (EGFR differs only in n_jobs).
        model=RandomForestRegressor(n_estimators=200,min_samples_leaf=2,max_features='sqrt',n_jobs=4,random_state=42)
        model.fit(X[train],y[train])
        sc=scores(y[test],model.predict(X[test]))
        rows.append({'target':t,'seed':seed,'train_n':len(train),'test_n':len(test),**sc})
        print(f'  RF repeat seed={seed} test_n={len(test):,} R2={sc["R2"]:.4f} MAE={sc["MAE"]:.4f}',flush=True)
    return rows


def main():
    parser=argparse.ArgumentParser(description='Audit Step 1–9 dataset/split integrity and subgroup performance.')
    parser.add_argument('--targets',nargs='+',choices=TARGETS,default=list(TARGETS))
    parser.add_argument('--repeat-splits',type=int,default=0,help='Optional RF refits on alternative scaffold splits, e.g. 3 (slow).')
    args=parser.parse_args()
    if not 0<=args.repeat_splits<=10:
        parser.error('--repeat-splits must be 0..10')
    OUT.mkdir(parents=True, exist_ok=True)
    summary=[]; bins=[]; sensitivity=[]; examples=[]; repeats=[]
    print('Step 11 — read-only baseline validation audit',flush=True)
    for t in args.targets:
        print(f'\n========== {t} ==========',flush=True)
        report, b, sens, ex, scaffolds, X, y=check_target(t)
        summary.append(report); bins.extend(b); sensitivity.extend(sens); examples.extend(ex)
        print(f'  samples={brief(report["samples"])} | split leakage={report["split_leakage_status"]}',flush=True)
        print(f'  train/test shared scaffolds={report["overlapping_scaffolds_train_test"]}, '
              f'structures={report["overlapping_canonical_smiles_train_test"]}, '
              f'fingerprints={report["overlapping_fingerprints_train_test"]}',flush=True)
        print(f'  test extreme flags={report["test_extreme_review_values"]}, '
              f'test labels outside train activity range={report["test_outside_training_activity_range"]}',flush=True)
        if args.repeat_splits:
            print(f'  Running {args.repeat_splits} alternative scaffold splits (RF only)...',flush=True)
            repeats.extend(repeat_rf_stability(t,X,y,scaffolds,args.repeat_splits))
    save_csv(OUT/'leakage_and_dataset_audit.csv',summary,list(summary[0]))
    save_csv(OUT/'activity_distribution_bins.csv',bins,list(bins[0]))
    save_csv(OUT/'test_subgroup_performance.csv',sensitivity,list(sensitivity[0]))
    if examples:
        save_csv(OUT/'cross_split_fingerprint_collision_examples.csv',examples,list(examples[0]))
    if repeats:
        save_csv(OUT/'repeat_scaffold_split_rf.csv',repeats,list(repeats[0]))
    bad=[s['target'] for s in summary if s['split_leakage_status']!='PASS' or s['duplicated_molecule_ids']]
    lines=['# Step 11 — ML Validation Audit','',
           'This is a **non-destructive** audit of pre-existing model inputs, scaffold assignments, and test predictions.',
           'No models were tuned or selected based on the original fixed test results.','',
           '## Audit checks','',
           '- Verified raw dataset, fingerprint matrix, metadata, activity labels, and split assignments align by row and molecule ID.',
           '- Recomputed Bemis–Murcko scaffolds; checked overlap between every pair of train, validation and test partitions.',
           '- Checked exact RDKit canonical structure overlap and identical 2048-bit Morgan fingerprint overlap across partitions. **Different molecules can have identical folded fingerprints**: this is a representation-level collision, not necessarily exact chemical-structure leakage.',
           '- Checked pIC50 distribution, test activity extrapolation, extreme-value flags, and performance by activity subgroup.',
           '- No salts or tautomers were standardized: alternative representations may remain a source of unmeasured chemical equivalence.',
           '', '## Summary', '',
           '| Target | Audit | Train/test scaffold overlap | Structure overlap | FP overlap | Test extremes | Test outside train label range |',
           '|---|---|---:|---:|---:|---:|---:|']
    for r in summary:
        lines.append(f'| {r["target"]} | {r["split_leakage_status"]} | {r["overlapping_scaffolds_train_test"]} | {r["overlapping_canonical_smiles_train_test"]} | {r["overlapping_fingerprints_train_test"]} | {r["test_extreme_review_values"]} | {r["test_outside_training_activity_range"]} |')
    lines.extend(['', '## Interpretation','',
                 '- Continuous pIC50 prediction is regression; **classification class imbalance** is not the appropriate metric. The activity-bin table quantifies rare-label representation instead.',
                 '- Cross-target structural overlap is not in itself leakage: models have different endpoint targets. Contamination is checked within each target.',
                 '- **REVIEW_FINGERPRINT_COLLISIONS** indicates identical model inputs for chemically distinct molecules across splits. Review the collision-examples CSV; a chemically scaffold-disjoint split may still contain representation-level duplication.',
                 '- A clean split does not prove transferability to generated triptanthrin analogs; use the Step 6–9 structural applicability warnings.',
                 '- The pIC50 value is a curated median of per-measurement transformed activity. Raw IC50 median may differ under this nonlinear transform.',
                 '- Model agreement does not quantify calibrated uncertainty. No external validation is claimed.',
                 '- Test-subgroup slices are descriptive diagnostics, not independent evaluation sets.',
                 '', '## Next checks','',
                 '- Confirm fifth vs sixth project target scope with teammates.',
                 '- Review source-unit consistency, molecule standardization, assay heterogeneity, and extreme activity records.',
                 '- Where feasible, use a separately sourced external experimental dataset for prospective validity.',
                 '- Document any selection or model hyperparameter tuning that used the same fixed test, and avoid further test-driven tuning.',
                 ''])
    if repeats:
        lines.extend(['## Alternative scaffold-split RF stability (exploratory)', '',
                      f'- Number of additional independently generated scaffold partitions per target: {args.repeat_splits}.',
                      '- Same Random Forest configuration as the original training pipeline. **XGBoost not repeated**.',
                      '- Since test compounds differ across runs, these estimates assess sensitivity to split composition, not the uncertainty of a fixed-sample estimate.',
                      '', '| Target | Mean test R² | SD test R² | Min | Max |', '|---|---:|---:|---:|---:|'])
        for t in args.targets:
            vals=np.asarray([r['R2'] for r in repeats if r['target']==t],dtype=float)
            lines.append(f'| {t} | {vals.mean():.4f} | {vals.std(ddof=1) if len(vals)>1 else float("nan"):.4f} | {vals.min():.4f} | {vals.max():.4f} |')
        lines.append('')
    (OUT/'validation_audit_report.md').write_text('\n'.join(lines),encoding='utf-8')
    print(f'\nStep 11 audit complete. Targets needing further leakage/ID review: {bad if bad else "none"}')
    print('Saved: reports\\final_validation\\leakage_and_dataset_audit.csv')
    print('       reports\\final_validation\\activity_distribution_bins.csv')
    print('       reports\\final_validation\\test_subgroup_performance.csv')
    print('       reports\\final_validation\\validation_audit_report.md')
    if examples:
        print('       reports\\final_validation\\cross_split_fingerprint_collision_examples.csv')
    if args.repeat_splits:
        print('       reports\\final_validation\\repeat_scaffold_split_rf.csv')
    print('Existing model and prediction files were not modified.')

if __name__=='__main__':
    main()
