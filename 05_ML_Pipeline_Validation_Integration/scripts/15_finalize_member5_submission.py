r"""Step 15 — Finish the Member 5 report and verify the team handover.

Place in scripts/ and run from the project root on Windows:
    .\.venv\Scripts\python.exe scripts/15_finalize_member5_submission.py

Inputs are your existing Step 12/14 outputs. No training or source data is changed.
Writes NEW files to reports/member5_submission/ and a notebook review checklist.
Requires matplotlib and numpy (already installed).
"""
from __future__ import annotations
import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports'
FINAL = REPORTS / 'member5_final'
OUT = REPORTS / 'member5_submission'
VALID = REPORTS / 'final_validation'
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
DISPLAY = {'PI3Kalpha': 'PI3Kα'}


def read_csv(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def asfloat(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')


def fmt(value, decimals=4):
    x = asfloat(value)
    return f'{x:.{decimals}f}' if math.isfinite(x) else 'N/A'


def collision_section(rows):
    look = {(r['target'], r['model'], r['subset']): r for r in rows}
    lines = [
        '### Numerical sensitivity results (saved Step 12 predictions)', '',
        '| Target | Affected test molecules | RF full R² | RF excluded R² | XGB full R² | XGB excluded R² |',
        '|---|---:|---:|---:|---:|---:|',
    ]
    changes = []
    for t in TARGETS:
        a = look[(t, 'RF', 'full_test')]
        b = look[(t, 'RF', 'no_train_identical_fingerprint')]
        c = look[(t, 'XGB', 'full_test')]
        d = look[(t, 'XGB', 'no_train_identical_fingerprint')]
        n = int(float(a['test_molecules_matching_train_fingerprint']))
        fraction = asfloat(a['test_matching_fraction']) * 100
        lines.append(f"| {DISPLAY.get(t,t)} | {n} ({fraction:.2f}%) | {fmt(a['R2'])} | {fmt(b['R2'])} | {fmt(c['R2'])} | {fmt(d['R2'])} |")
        changes.extend([abs(asfloat(b['R2']) - asfloat(a['R2'])), abs(asfloat(d['R2']) - asfloat(c['R2']))])
    largest = max(changes)
    lines += [
        '', f'Largest absolute change in R² after excluding matching-fingerprint test molecules: **{largest:.4f}**.',
        'The excluded subset differs in composition from the original test set; changes in R² do not directly quantify leakage-induced bias.',
        'This is a descriptive fixed-model sensitivity analysis, not independent prospective validation.', ''
    ]
    return '\n'.join(lines), largest


def similarity_counts():
    summary = REPORTS / 'applicability_domain' / 'similarity_summary.csv'
    if summary.exists():
        rows = read_csv(summary)
        if rows:
            cols = {k.lower(): k for k in rows[0].keys()}
            target_col = cols.get('target')
            variants = [
                ('high_similarity', 'moderate_similarity', 'low_similarity_review'),
                ('high_similarity', 'moderate_similarity', 'low_similarity'),
                ('high_count', 'moderate_count', 'low_count'),
                ('high_similarity_count', 'moderate_similarity_count', 'low_similarity_count'),
                ('HIGH_SIMILARITY_SCREEN', 'MODERATE_SIMILARITY_SCREEN', 'LOW_SIMILARITY_EXPLORATORY_ONLY'),
            ]
            if target_col:
                for variant in variants:
                    if all(s.lower() in cols for s in variant):
                        got = {r[target_col]: [int(float(r[cols[s.lower()]])) for s in variant] for r in rows}
                        if all(t in got for t in TARGETS):
                            return got
    # Fall back to counting existing predictions, without loading them all into RAM.
    path = REPORTS / 'analog_predictions' / 'all_target_predictions.csv'
    if path.exists():
        counts = Counter()
        with path.open(encoding='utf-8-sig', newline='') as f:
            reader = csv.DictReader(f)
            cols = {k.lower(): k for k in (reader.fieldnames or [])}
            target = cols.get('target')
            tier = next((cols[k] for k in ('priority_tier', 'applicability_tier', 'similarity_tier') if k in cols), None)
            if target and tier:
                for row in reader:
                    counts[row[target], row[tier]] += 1
                return {t: [counts[t,'HIGH_SIMILARITY_SCREEN'], counts[t,'MODERATE_SIMILARITY_SCREEN'],
                            counts[t,'LOW_SIMILARITY_EXPLORATORY_ONLY']] for t in TARGETS}
    return None


def fix_similarity_figure():
    counts = similarity_counts()
    if counts is None:
        return False
    values = np.asarray([counts[t] for t in TARGETS], dtype=float)
    fractions = 100 * values / np.maximum(values.sum(axis=1, keepdims=True), 1)
    fig, ax = plt.subplots(figsize=(10.6, 5.1))
    fig.subplots_adjust(right=0.79, top=0.88, bottom=0.13)
    bottom = np.zeros(len(TARGETS))
    for i, (name, color) in enumerate((('High', '#26734d'), ('Moderate', '#e5b155'), ('Low', '#bd6868'))):
        ax.bar(np.arange(len(TARGETS)), fractions[:, i], bottom=bottom, color=color, label=name)
        bottom += fractions[:, i]
    ax.set_xticks(np.arange(len(TARGETS)), [DISPLAY.get(t,t) for t in TARGETS])
    ax.set_ylim(0, 100)
    ax.set_ylabel('Analogs (%)')
    ax.set_title('Applicability-domain structural similarity tiers')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), frameon=True)
    ax.grid(axis='y', alpha=0.12)
    fig.savefig(OUT / 'fig03_applicability_coverage.png', dpi=180)
    plt.close(fig)
    return True


def notebook_check():
    notebook = ROOT / 'notebooks' / 'Member5_results_review.ipynb'
    if not notebook.exists():
        return 'MISSING: notebooks/Member5_results_review.ipynb', 0
    try:
        doc = json.loads(notebook.read_text(encoding='utf-8'))
        cells = doc.get('cells', [])
        assert doc.get('nbformat') == 4 and isinstance(cells, list)
        code_cells = sum(c.get('cell_type') == 'code' for c in cells)
        return f'Valid Jupyter structure: {len(cells)} cells, {code_cells} code cells; EXECUTION NOT VERIFIED', code_cells
    except (ValueError, AssertionError, TypeError) as ex:
        return f'NOTEBOOK INVALID: {ex}', 0


def main():
    original = FINAL / 'Member5_ML_Final_Report.md'
    collision = VALID / 'fingerprint_collision_sensitivity_metrics.csv'
    required = [original, collision, FINAL / 'final_metrics_summary.csv']
    absent = [str(p.relative_to(ROOT)) for p in required if not p.is_file()]
    if absent:
        raise SystemExit('Required inputs missing. Run Step 12 and Step 14 first:\n' + '\n'.join('  - ' + x for x in absent))
    data = read_csv(collision)
    required_cols = {'target', 'model', 'subset', 'R2', 'test_molecules_matching_train_fingerprint', 'test_matching_fraction'}
    if not data or not required_cols.issubset(data[0]):
        raise SystemExit('Step 12 metrics have unexpected columns. Please upload that CSV rather than changing it manually.')

    OUT.mkdir(parents=True, exist_ok=True)
    report = original.read_text(encoding='utf-8-sig')
    section, max_delta = collision_section(data)
    marker = 'Source: `reports/final_validation/fingerprint_collision_sensitivity_metrics.csv`.'
    if marker not in report:
        raise SystemExit('Could not locate Section 5 insertion point in Step 14 report; refusing an unsafe edit.')
    report = report.replace(marker, marker + '\n\n' + section, 1)
    report = report.replace('## 10. Reproducibility and handover checklist',
        '## 10. Reproducibility and handover checklist')
    # Clearly differentiate validation and unverified items, not a claim of full release readiness.
    report += ('\n\n## 11. Submission verification and pending decisions\n\n'
               '- The final report was assembled from previously saved outputs; no models were retrained.\n'
               '- **Team decision pending:** confirm whether the sixth proposed target was formally excluded.\n'
               '- **Scientific review pending:** confirm published-reference identities and parent-compound activity comparisons.\n'
               '- **Notebook:** confirm review cells execute with the project virtual environment.\n'
               '- **Reproducibility:** preserve requirements, source data versions, model files and processing scripts.\n'
               '- **Uncertainty:** no calibrated confidence intervals or external experimental validation are claimed.\n')
    (OUT / 'Member5_ML_Final_Report_REVIEWED.md').write_text(report, encoding='utf-8')

    shutil.copy2(FINAL / 'final_metrics_summary.csv', OUT / 'final_metrics_summary.csv')
    for name in ('fig01_original_rf_vs_xgb.png', 'fig02_rf_scaffold_robustness.png'):
        src = FINAL / name
        if src.exists():
            shutil.copy2(src, OUT / name)
    fixed = fix_similarity_figure()
    if not fixed and (FINAL / 'fig03_applicability_coverage.png').exists():
        shutil.copy2(FINAL / 'fig03_applicability_coverage.png', OUT / 'fig03_applicability_coverage.png')

    status, ncode = notebook_check()
    deliverables = {
        'ML scripts': ROOT / 'scripts',
        'saved trained models': ROOT / 'models',
        'original target metrics': REPORTS / 'all_targets_model_comparison.csv',
        'applicability summary': REPORTS / 'applicability_domain' / 'similarity_summary.csv',
        'full analog predictions': REPORTS / 'analog_predictions' / 'all_target_predictions.csv',
        'docking shortlist': REPORTS / 'diversity_shortlist' / 'master_docking_handover.csv',
        'Step 11 leakage audit': VALID / 'leakage_and_dataset_audit.csv',
        'Step 12 collision sensitivity': collision,
        'Step 13 robustness': VALID / 'repeated_scaffold_robustness' / 'rf_robustness_summary.csv',
        'Step 14 notebook': ROOT / 'notebooks' / 'Member5_results_review.ipynb',
    }
    checks = ['# Member 5 — Submission checklist', '',
        'Run Step 15 from project root; this is a file-integrity and completeness check, not a scientific certification.', '',
        f'- Step 12 maximum absolute ΔR² (collision-excluded minus full): **{max_delta:.4f}**',
        f'- Notebook: **{status}**',
        f'- Applicability plot regenerated with external legend: **{"YES" if fixed else "NO (original retained if present)"}**', '',
        '## File inventory', '']
    for name, path in deliverables.items():
        checks.append(f'- [{"x" if path.exists() else " "}] {name}: `{path.relative_to(ROOT)}`')
    checks += ['', '## Human checks before handover', '',
        '- [ ] Open and run each code cell of `notebooks/Member5_results_review.ipynb` in the correct Python environment.',
        '- [ ] Verify scientific terminology, figure captions, citations and assay metadata with Member 4.',
        '- [ ] Confirm five versus six targets with the team and guide.',
        '- [ ] Confirm the published references and parent triptanthrin activity comparison; do not claim experimental hits.',
        '- [ ] Share report, figures, saved source code, full predictions and audited 50-entry handover.',
        '- [ ] Do not include virtual environment folders or large raw data archives in GitHub by default.', '']
    (OUT / 'Member5_submission_checklist.md').write_text('\n'.join(checks), encoding='utf-8')
    print('Step 15 complete — final report review package generated. NO retraining.')
    print(' Report   :', (OUT / 'Member5_ML_Final_Report_REVIEWED.md').relative_to(ROOT))
    print(' Checklist:', (OUT / 'Member5_submission_checklist.md').relative_to(ROOT))
    print(' Plot     :', (OUT / 'fig03_applicability_coverage.png').relative_to(ROOT) if (OUT / 'fig03_applicability_coverage.png').exists() else 'Not available')
    print(f' Collision sensitivity maximum |ΔR²|: {max_delta:.4f}')
    print(' Notebook :', status)
    print(f' Missing  : {sum(not p.exists() for p in deliverables.values())} of {len(deliverables)} inventory items')
    print('Original reports, models, predictions and shortlists unchanged.')

if __name__ == '__main__':
    main()
