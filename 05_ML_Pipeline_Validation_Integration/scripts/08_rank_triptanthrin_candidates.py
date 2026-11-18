"""Step 8 — transparent, target-specific triptanthrin prioritization.

Run from your existing project:
    python scripts/08_rank_triptanthrin_candidates.py

No model training is performed. Separate target rankings are NOT comparable
potency claims; predictions and structural similarity are uncalibrated.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'reports' / 'analog_predictions' / 'all_target_predictions.csv'
OUTPUT = ROOT / 'reports' / 'candidate_prioritization'
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
TIER_ORDER = {'HIGH_SIMILARITY_SCREEN': 0, 'MODERATE_SIMILARITY_SCREEN': 1,
              'LOW_SIMILARITY_EXPLORATORY_ONLY': 2}
TOP_N = 25
# Heuristic review cutoff, NOT a calibrated uncertainty interval.
MAX_MODEL_DIFFERENCE = 0.50


def load_csv(path):
    if not path.exists():
        raise FileNotFoundError(f'Input not found: {path}\nRun Step 7 first.')
    with path.open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f'Empty input: {path}')
    required = {'Analog_ID', 'Target', 'SMILES', 'Decision', 'Category',
                'QED', 'SA_score', 'Max_train_Tanimoto', 'Similarity_band',
                'Pred_RF_pIC50', 'Pred_XGB_pIC50', 'Model_disagreement_abs',
                'Priority_tier', 'RF_outside_training_pIC50_range',
                'XGB_outside_training_pIC50_range'}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    return rows


def asfloat(row, key):
    try:
        return float(row[key])
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{row["Analog_ID"]}, {row["Target"]}: invalid {key}') from exc


def boolean(s):
    if s in ('True', 'true', '1'): return True
    if s in ('False', 'false', '0'): return False
    raise ValueError(f'Invalid boolean value {s!r}')


def save_csv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = load_csv(INPUT)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    seen = set()
    by_target = {t: [] for t in TARGETS}
    for row in rows:
        target = row['Target']
        if target not in by_target:
            raise ValueError(f'Unexpected target: {target}')
        identity = (target, row['Analog_ID'])
        if identity in seen:
            raise ValueError(f'Duplicate analog–target combination: {identity}')
        seen.add(identity)
        if row['Priority_tier'] not in TIER_ORDER:
            raise ValueError(f'Unknown tier: {row["Priority_tier"]}')
        if row['Decision'] not in {'RETAIN', 'FLAG'}:
            raise ValueError(f'Unexpected filter Decision: {row["Decision"]}')
        rf = asfloat(row, 'Pred_RF_pIC50')
        xgb = asfloat(row, 'Pred_XGB_pIC50')
        difference = abs(rf - xgb)
        if abs(difference - asfloat(row, 'Model_disagreement_abs')) > 0.00005:
            raise ValueError(f'Model disagreement mismatch: {identity}')
        row['Consensus_pIC50'] = round((rf + xgb) / 2, 5)
        row['Model_disagreement_abs'] = round(difference, 5)
        row['QED'] = asfloat(row, 'QED')
        row['SA_score'] = asfloat(row, 'SA_score')
        row['Max_train_Tanimoto'] = asfloat(row, 'Max_train_Tanimoto')
        outside = (boolean(row['RF_outside_training_pIC50_range']) or
                   boolean(row['XGB_outside_training_pIC50_range']))
        row['Review_note'] = ('OUTSIDE_TRAINING_LABEL_RANGE; ' if outside else '') + (
            'RF_XGB_DISAGREE' if difference > MAX_MODEL_DIFFERENCE else 'Model_agreement_within_heuristic_cutoff')
        # Explicitly disallow low-similarity and out-of-range compounds from screening shortlist.
        row['Screen_eligible'] = (row['Priority_tier'] != 'LOW_SIMILARITY_EXPLORATORY_ONLY'
                                  and difference <= MAX_MODEL_DIFFERENCE and not outside)
        row['Screen_eligible'] = str(row['Screen_eligible'])
        by_target[target].append(row)

    counts = {t: len(group) for t, group in by_target.items()}
    if len(set(counts.values())) != 1:
        raise ValueError(f'Unequal analog counts by target: {counts}')
    num_analogs = next(iter(counts.values()))
    if len(rows) != num_analogs * len(TARGETS):
        raise ValueError('Unexpected target coverage')
    for target in TARGETS:
        if {r['Analog_ID'] for r in by_target[target]} != {r['Analog_ID'] for r in by_target[TARGETS[0]]}:
            raise ValueError('Different analog sets across targets')

    columns = ['Rank', 'Analog_ID', 'Target', 'SMILES', 'Category', 'Decision',
               'Priority_tier', 'Max_train_Tanimoto', 'Consensus_pIC50',
               'Pred_RF_pIC50', 'Pred_XGB_pIC50', 'Model_disagreement_abs',
               'QED', 'SA_score', 'Nearest_train_chembl_id', 'Review_note', 'Screen_eligible']
    all_candidates = []
    all_exploratory = []
    summary = []
    for target in TARGETS:
        group = by_target[target]
        # Never pool pIC50 scores across proteins. Within target, favor retained
        # molecules, then high similarity, then predicted potency and agreement.
        candidates = [r for r in group if r['Screen_eligible'] == 'True']
        candidates.sort(key=lambda r: (r['Decision'] != 'RETAIN',
                                        TIER_ORDER[r['Priority_tier']],
                                        -r['Consensus_pIC50'],
                                        r['Model_disagreement_abs'],
                                        -r['QED'], r['SA_score'], r['Analog_ID']))
        exploratory = [r for r in group if r['Screen_eligible'] != 'True']
        exploratory.sort(key=lambda r: (TIER_ORDER[r['Priority_tier']],
                                        r['Decision'] != 'RETAIN',
                                        -r['Consensus_pIC50'], r['Analog_ID']))
        for i, r in enumerate(candidates, 1): r['Rank'] = i
        for i, r in enumerate(exploratory, 1): r['Rank'] = i
        save_csv(OUTPUT / f'{target}_screen_candidates.csv', candidates, columns)
        save_csv(OUTPUT / f'{target}_top_{TOP_N}.csv', candidates[:TOP_N], columns)
        save_csv(OUTPUT / f'{target}_exploratory_review.csv', exploratory, columns)
        all_candidates.extend(candidates[:TOP_N])
        all_exploratory.extend(exploratory[:TOP_N])
        bands = Counter(r['Priority_tier'] for r in group)
        notes = Counter(r['Decision'] for r in candidates)
        summary.append({'Target': target, 'All_predictions': len(group),
                        'High_similarity': bands['HIGH_SIMILARITY_SCREEN'],
                        'Moderate_similarity': bands['MODERATE_SIMILARITY_SCREEN'],
                        'Low_similarity': bands['LOW_SIMILARITY_EXPLORATORY_ONLY'],
                        'Eligible_for_screening': len(candidates),
                        'Eligible_RETAIN': notes['RETAIN'],
                        'Eligible_FLAG': notes['FLAG'],
                        'Exploratory_or_review': len(exploratory),
                        'Top_file': f'{target}_top_{TOP_N}.csv'})
        print(f'{target:10s} eligible={len(candidates):5,d} | exploratory/review={len(exploratory):5,d} | '
              f'high={bands["HIGH_SIMILARITY_SCREEN"]:4,d} moderate={bands["MODERATE_SIMILARITY_SCREEN"]:5,d} '
              f'low={bands["LOW_SIMILARITY_EXPLORATORY_ONLY"]:5,d}', flush=True)
    save_csv(OUTPUT / 'top_candidates_all_targets.csv', all_candidates, columns)
    save_csv(OUTPUT / 'top_exploratory_for_review.csv', all_exploratory, columns)
    save_csv(OUTPUT / 'screening_summary.csv', summary, list(summary[0]))
    (OUTPUT / 'methodology.json').write_text(json.dumps({
        'purpose': 'Computational candidate prioritization, not proof of anticancer activity',
        'source': str(INPUT.relative_to(ROOT)),
        'per_target_only': 'Do not interpret numerical pIC50 values as cross-target selectivity or clinical potency',
        'consensus': 'Arithmetic mean of RF and XGBoost pIC50 predictions; not an ensemble with validated calibration',
        'screen_eligibility': ['At least MODERATE_SIMILARITY_SCREEN',
                               f'Absolute RF-XGB difference <= {MAX_MODEL_DIFFERENCE} (heuristic)',
                               'Neither prediction outside observed training label range'],
        'rank_order': ['Drug-likeness Decision RETAIN before FLAG',
                       'HIGH before MODERATE similarity',
                       'Higher RF/XGB mean predicted pIC50',
                       'Lower model disagreement', 'Higher QED', 'Lower SA score'],
        'warning': ['Similarity cutoffs and model agreement cutoffs are arbitrary review heuristics, not validated probabilities',
                    'Even screened candidates may be out of domain; compare to training-neighbor similarity distributions',
                    'SA and QED are secondary annotations, not synthetic feasibility or toxicity guarantees',
                    'No docking, experimental testing or retrosynthesis performed',
                    'Review structural-alert flags and prioritize chemical diversity before wet-lab or docking selection',
                    'Do not use repeated test-set comparisons to tune the selection rules']
    }, indent=2), encoding='utf-8')
    print(f'\nStep 8 complete: {len(rows):,} target–analog rows evaluated.')
    print(f'Outputs: {OUTPUT.relative_to(ROOT)}')
    print('Interpret target-wise shortlists as computational hypotheses only.')

if __name__ == '__main__':
    main()
