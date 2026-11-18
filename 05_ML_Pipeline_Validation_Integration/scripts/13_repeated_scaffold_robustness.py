"""Step 13 — repeated scaffold-disjoint RF robustness validation (read-only).

Place this single file in your project's scripts/ folder, then run from project root:
    .\\.venv\\Scripts\\python.exe scripts/13_repeated_scaffold_robustness.py

Default: five targets x three NEW scaffold splits (seeds 13, 23, 37),
200-tree RandomForest per run, same RF hyperparameters as Steps 3/5.
No changes to existing seed-42 model artifacts, original splits, predictions or shortlists.

Optional quick check:
    .\\.venv\\Scripts\\python.exe scripts/13_repeated_scaffold_robustness.py --targets TDO2 --seeds 13 --trees 10
Do not use the quick check's reduced-tree metrics as the scientific result.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from datetime import datetime, timezone

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'data' / 'processed'
OUTPUT = ROOT / 'reports' / 'final_validation' / 'repeated_scaffold_robustness'
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
DEFAULT_SEEDS = (13, 23, 37)  # exclude original seed=42


def load_target(target: str):
    arrpath = INPUT / f'{target}_ECFP4_2048.npz'
    metapath = INPUT / f'{target}_fingerprint_metadata.csv'
    if not arrpath.is_file() or not metapath.is_file():
        raise FileNotFoundError(f'{target}: missing {arrpath} or {metapath}. Run Steps 2/5 first.')
    with np.load(arrpath) as arr:
        X = arr['X']
        y = arr['y'].astype(np.float64)
    with metapath.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    if X.ndim != 2 or X.shape[1] != 2048 or not (len(X) == len(y) == len(rows)):
        raise ValueError(f'{target}: fingerprint/label/metadata shapes do not match')
    if not np.isfinite(y).all():
        raise ValueError(f'{target}: non-finite pIC50 values')
    ids, scaffolds = [], []
    for i, row in enumerate(rows):
        if int(row['fingerprint_row']) != i:
            raise ValueError(f'{target}: fingerprint_row misaligned at {i}')
        if abs(float(row['pIC50']) - float(y[i])) > 1e-6:
            raise ValueError(f'{target}: metadata pIC50 misaligned at {i}')
        mol = Chem.MolFromSmiles(row['canonical_smiles'])
        if mol is None:
            raise ValueError(f'{target}: SMILES invalid at row {i}')
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffolds.append(Chem.MolToSmiles(core) or '__NO_MURCKO_SCAFFOLD__')
        ids.append(row['molecule_chembl_id'])
    if len(set(ids)) != len(ids):
        raise ValueError(f'{target}: duplicate molecular IDs')
    return X, y, np.asarray(scaffolds)


def make_split(y: np.ndarray, groups: np.ndarray, seed: int):
    # Exactly the two-stage GroupShuffleSplit strategy used in Steps 3 and 5.
    dummy = np.arange(len(y))
    first = GroupShuffleSplit(n_splits=1, test_size=0.30, random_state=seed)
    train, hold = next(first.split(dummy, y, groups=groups))
    second = GroupShuffleSplit(n_splits=1, test_size=0.50, random_state=seed)
    val_rel, test_rel = next(second.split(hold, y[hold], groups=groups[hold]))
    val, test = hold[val_rel], hold[test_rel]
    parts = (train, val, test)
    if not all(len(s) >= 2 for s in parts):
        raise ValueError(f'Seed {seed} created an excessively small partition')
    scaffold_sets = [set(groups[s]) for s in parts]
    if (scaffold_sets[0] & scaffold_sets[1] or scaffold_sets[0] & scaffold_sets[2]
            or scaffold_sets[1] & scaffold_sets[2]):
        raise ValueError(f'Scaffold overlap at seed {seed}')
    if not np.array_equal(np.sort(np.concatenate(parts)), np.arange(len(y))):
        raise ValueError(f'Incomplete split at seed {seed}')
    return train, val, test


def scores(y_true, y_pred):
    return {'MAE': float(mean_absolute_error(y_true, y_pred)),
            'RMSE': float(np.sqrt(mean_squared_error(y_true, y_pred))),
            'R2': float(r2_score(y_true, y_pred))}


def get_collision_count(X, train, test):
    # Count affected TEST MOLECULES, not matching distinct fingerprints or pairs.
    # np.packbits avoids constructing large (n_train*n_test) matrices.
    train_set = set(map(bytes, np.packbits(X[train], axis=1)))
    return sum(bytes(fp) in train_set for fp in np.packbits(X[test], axis=1))


def save_csv(path: Path, records: list[dict], fields: list[str]):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--targets', nargs='+', choices=TARGETS, default=list(TARGETS))
    parser.add_argument('--seeds', type=int, nargs='+', default=list(DEFAULT_SEEDS))
    parser.add_argument('--trees', type=int, default=200, help='Use 200 for final result; smaller only for smoke tests')
    args = parser.parse_args()
    if len(set(args.seeds)) != len(args.seeds):
        parser.error('Seeds must be distinct')
    if 42 in args.seeds:
        parser.error('Seed 42 is the original fixed split; select only NEW seeds')
    if args.trees < 1:
        parser.error('--trees must be >= 1')
    return args


def main():
    args = parse_arguments()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    full_run = args.trees == 200 and args.seeds == list(DEFAULT_SEEDS) and list(args.targets) == list(TARGETS)
    if not full_run:
        print('NOTE: Non-default configuration. For project report, run the full default command.', flush=True)
    print('Step 13 — repeated scaffold split RF robustness; preserves existing models', flush=True)
    print(f'Targets={args.targets} | seeds={args.seeds} | trees={args.trees}', flush=True)
    results = []
    summaries = []
    for target in args.targets:
        X, y, scaffolds = load_target(target)
        print(f'\n========== {target} | molecules={len(y):,} | scaffolds={len(set(scaffolds)):,} ==========', flush=True)
        for seed in args.seeds:
            train, val, test = make_split(y, scaffolds, seed)
            model = RandomForestRegressor(n_estimators=args.trees, min_samples_leaf=2,
                                          max_features='sqrt', n_jobs=4, random_state=42)
            model.fit(X[train], y[train])
            y_pred = model.predict(X[test])
            m = scores(y[test], y_pred)
            base = scores(y[test], np.full(len(test), np.median(y[train])))
            overlap = get_collision_count(X, train, test)
            record = {
                'target': target, 'seed': seed, 'n_train': len(train),
                'n_validation': len(val), 'n_test': len(test),
                'scaffolds_train': len(set(scaffolds[train])),
                'scaffolds_validation': len(set(scaffolds[val])),
                'scaffolds_test': len(set(scaffolds[test])),
                'test_affected_fingerprint_overlap': overlap,
                'test_affected_fingerprint_overlap_pct': 100 * overlap / len(test),
                'test_extreme_pIC50_flags': int(((y[test] < 3) | (y[test] > 12)).sum()),
                'test_below_5_pct': 100 * np.mean(y[test] < 5),
                'test_above_equal_9_pct': 100 * np.mean(y[test] >= 9),
                'MAE': m['MAE'], 'RMSE': m['RMSE'], 'R2': m['R2'],
                'baseline_MAE': base['MAE'], 'baseline_RMSE': base['RMSE'],
                'baseline_R2': base['R2'],
            }
            results.append(record)
            print(f'  seed={seed:<3} train={len(train):5,} val={len(val):5,} test={len(test):5,}'
                  f' | MAE={m["MAE"]:.4f} RMSE={m["RMSE"]:.4f} R2={m["R2"]:.4f}'
                  f' | FP-overlap test={overlap}', flush=True)
            del model
        rs = [r for r in results if r['target'] == target]
        summary = {'target': target, 'repeats': len(rs), 'trees': args.trees}
        for metric in ('MAE', 'RMSE', 'R2'):
            vals = np.asarray([r[metric] for r in rs], dtype=float)
            summary.update({f'{metric}_mean': float(np.mean(vals)),
                            f'{metric}_std': float(np.std(vals, ddof=1)) if len(vals) > 1 else float('nan'),
                            f'{metric}_min': float(vals.min()),
                            f'{metric}_max': float(vals.max())})
        summaries.append(summary)
        print(f'  SUMMARY mean R2={summary["R2_mean"]:.4f}'
              f' | std={summary["R2_std"]:.4f}'
              f' | range=[{summary["R2_min"]:.4f}, {summary["R2_max"]:.4f}]', flush=True)
        del X, y, scaffolds

    suffix = '' if full_run else '_custom_run'
    detail_path = OUTPUT / f'rf_repeated_split_metrics{suffix}.csv'
    summary_path = OUTPUT / f'rf_robustness_summary{suffix}.csv'
    save_csv(detail_path, results, list(results[0].keys()))
    save_csv(summary_path, summaries, list(summaries[0].keys()))
    fig, ax = plt.subplots(figsize=(9, 5))
    target_order = list(args.targets)
    means = [next(r for r in summaries if r['target'] == t)['R2_mean'] for t in target_order]
    stds = [next(r for r in summaries if r['target'] == t)['R2_std'] for t in target_order]
    stds = [v if np.isfinite(v) else 0 for v in stds]
    x = np.arange(len(target_order))
    ax.errorbar(x, means, yerr=stds, fmt='o', capsize=5, label='mean ± 1 sample SD')
    for i, target in enumerate(target_order):
        vals = [r['R2'] for r in results if r['target'] == target]
        ax.scatter([i] * len(vals), vals, s=20, alpha=0.65)
    ax.set_xticks(x, target_order)
    ax.set_ylabel('R² on new scaffold-held-out test split')
    ax.set_title('RF scaffold split robustness across seeds (new refits)')
    ax.axhline(0, linestyle='--', linewidth=1, alpha=0.6)
    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(OUTPUT / f'rf_robustness_r2{suffix}.png', dpi=170)
    plt.close(fig)
    report = [
        '# Step 13 — Repeated scaffold-split robustness (Random Forest)', '',
        f'Run UTC: {datetime.now(timezone.utc).isoformat()}',
        f'RDKit: {rdBase.rdkitVersion}; scikit-learn: {sklearn_version}',
        f'Targets: {", ".join(args.targets)}',
        f'New seeds: {args.seeds}; n_estimators: {args.trees}',
        'Same RF settings as baseline: min_samples_leaf=2, max_features=sqrt, random_state=42; CPU n_jobs=4.',
        'Scaffold-disjoint split: 70/15/15 approximate, two-stage GroupShuffleSplit by Murcko scaffold; exact sample proportions can vary.',
        'Each model fits ONLY its own new train partition; validation/test labels are not used in fitting.',
        'Original seed-42 saved models, test predictions, candidate rankings, and docking handover were untouched.',
        '', '## Findings', '',
        '| Target | Repeats | Mean test MAE ± SD | Mean test RMSE ± SD | Mean test R² ± SD | Test R² range |',
        '|---|---:|---:|---:|---:|---:|',
    ]
    def fm(value):
        return f'{value:.4f}' if np.isfinite(value) else 'N/A (one repeat)'
    for r in summaries:
        report.append(f'| {r["target"]} | {r["repeats"]} | {fm(r["MAE_mean"])} ± {fm(r["MAE_std"])}'
                      f' | {fm(r["RMSE_mean"])} ± {fm(r["RMSE_std"])}'
                      f' | {fm(r["R2_mean"])} ± {fm(r["R2_std"])}'
                      f' | {fm(r["R2_min"])}–{fm(r["R2_max"])} |')
    report.extend(['', '## Interpretation and limitations', '',
        '- These are robustness refits on alternative scaffold splits, not external or prospective experimental validation.',
        '- Only Random Forest was repeated; do not report XGBoost as robustly revalidated from this script.',
        '- Different runs have different test structures, sizes and activity distributions, so variability includes split composition.',
        '- Identical fingerprint vectors may occur across distinct structures even in scaffold-disjoint splits; overlap counts are reported.',
        '- The model is trained exclusively with target-specific ChEMBL labels; assay heterogeneity and salts/tautomers are not resolved.',
        '- Three repeats illustrate variation but do not yield a calibrated uncertainty interval.',
        '- Retain original seed-42 metrics as baseline; avoid model tuning or ranking based on any new held-out test scores.',
        '- Structural applicability to generated triptanthrin molecules remains a separate, major concern.',
    ])
    if not full_run:
        report.extend(['', '**Non-default/smoke-test configuration — not the final report output.**'])
    (OUTPUT / f'rf_repeated_split_report{suffix}.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    print(f'\nStep 13 complete. Results saved under {OUTPUT.relative_to(ROOT)}', flush=True)
    print('Original models/predictions untouched. Review variation across seeds before making final claims.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'Step 13 failed: {type(exc).__name__}: {exc}', file=sys.stderr, flush=True)
        raise
