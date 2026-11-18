"""Step 14 — Assemble reproducible Member 5 ML report, figures and review notebook.

Run from the project root:
    .\\.venv\\Scripts\\python.exe scripts/14_generate_member5_final_report.py

READ ONLY input files under reports/. Writes only to reports/member5_final/ and
notebooks/Member5_results_review.ipynb; never retrains or modifies models.
Requires numpy + matplotlib (already installed in earlier steps).
"""
from __future__ import annotations
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports'
OUT = REPORTS / 'member5_final'
NOTEBOOKS = ROOT / 'notebooks'
TARGETS = ('EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2')
LABELS = {'PI3Kalpha': 'PI3Kα'}

INPUT = {
    'model_metrics': REPORTS / 'all_targets_model_comparison.csv',
    'leakage': REPORTS / 'final_validation' / 'leakage_and_dataset_audit.csv',
    'collision': REPORTS / 'final_validation' / 'fingerprint_collision_sensitivity_metrics.csv',
    'robustness': REPORTS / 'final_validation' / 'repeated_scaffold_robustness' / 'rf_robustness_summary.csv',
    'ad_summary': REPORTS / 'applicability_domain' / 'similarity_summary.csv',
    'prediction_summary': REPORTS / 'analog_predictions' / 'prediction_summary.csv',
    'screening': REPORTS / 'candidate_prioritization' / 'screening_summary.csv',
    'diversity': REPORTS / 'diversity_shortlist' / 'diversity_summary.csv',
    'handover': REPORTS / 'diversity_shortlist' / 'master_docking_handover.csv',
}

def read_csv(path):
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def number(x, default=float('nan')):
    try:
        v = float(x)
        return v if math.isfinite(v) else default
    except (TypeError, ValueError):
        return default

def label(t):
    return LABELS.get(t, t)

def fmt(v, nd=4):
    v = number(v)
    return f'{v:.{nd}f}' if math.isfinite(v) else 'N/A'

def pick(rows, target, model, split='test'):
    for r in rows:
        if r.get('target') == target and r.get('model','').lower().replace(' ','') == model.lower().replace(' ','') and r.get('split') == split:
            return r
    raise ValueError(f'Missing results for {target}, {model}, {split}')

def indexed(rows, key='target'):
    return {r[key]: r for r in rows}

def save_csv(path, rows):
    if not rows:
        return
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def make_comparison_chart(metrics):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    x = np.arange(len(TARGETS)); width = .36
    for ax, metric, title in [(axes[0], 'R2', 'Test R² (higher is better)'),
                              (axes[1], 'MAE', 'Test MAE (lower is better)')]:
        for shift, model, display in [(-width/2, 'Random Forest', 'Random Forest'),
                                      (width/2, 'XGBoost', 'XGBoost')]:
            values = [number(pick(metrics, t, model)[metric]) for t in TARGETS]
            ax.bar(x + shift, values, width=width, label=display)
        ax.set_xticks(x, [label(t) for t in TARGETS])
        ax.set_title(title)
        ax.set_ylabel(metric + (' (pIC50 units)' if metric == 'MAE' else ''))
        ax.grid(axis='y', alpha=.18)
    axes[0].legend(fontsize=9)
    fig.suptitle('Original scaffold-disjoint held-out test — one split per target', fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / 'fig01_original_rf_vs_xgb.png', dpi=180)
    plt.close(fig)

def make_robustness_chart(metrics, robustness):
    rt = indexed(robustness)
    x = np.arange(len(TARGETS)); orig = [number(pick(metrics,t,'Random Forest')['R2']) for t in TARGETS]
    means = [number(rt[t]['R2_mean']) for t in TARGETS]
    std = [number(rt[t]['R2_std'], 0) for t in TARGETS]
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.scatter(x-.11, orig, s=75, marker='s', label='Original RF test R²', zorder=4)
    ax.errorbar(x+.11, means, yerr=std, fmt='o', capsize=5,
                label='RF new splits: mean ± 1 sample SD', zorder=4)
    ax.set_xticks(x, [label(t) for t in TARGETS]); ax.set_ylabel('Test R²')
    ax.set_title('Original performance vs three new scaffold split refits')
    ax.grid(axis='y', alpha=.2); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / 'fig02_rf_scaffold_robustness.png', dpi=180); plt.close(fig)

def make_similarity_chart(sim):
    keys = {k.lower():k for k in sim[0]}
    target_key = keys.get('target', 'Target')
    # Step 6 summary format may vary. Prefer the exact integer tiers when available;
    # otherwise derive tier counts from Step 7 long prediction output.
    def tier_counts_from_summary(r):
        variants = [
            ('high_similarity', 'moderate_similarity', 'low_similarity_review'),
            ('high_similarity', 'moderate_similarity', 'low_similarity'),
            ('high_count', 'moderate_count', 'low_count'),
            ('high_similarity_count', 'moderate_similarity_count', 'low_similarity_count'),
            ('HIGH_SIMILARITY_SCREEN','MODERATE_SIMILARITY_SCREEN','LOW_SIMILARITY_EXPLORATORY_ONLY'),
        ]
        lower = {k.lower():k for k in r}
        for v in variants:
            if all(k.lower() in lower for k in v):
                return [number(r[lower[k.lower()]],0) for k in v]
        return None
    by_target={r[target_key]:r for r in sim}
    tiers={t:tier_counts_from_summary(by_target[t]) for t in TARGETS if t in by_target}
    if any(x is None for x in tiers.values()) or len(tiers)!=5:
        longpath=REPORTS/'analog_predictions'/'all_target_predictions.csv'
        if not longpath.exists():
            return False
        with longpath.open(encoding='utf-8-sig',newline='') as f:
            rdr=csv.DictReader(f)
            tkey=next((x for x in rdr.fieldnames if x.lower()=='target'),None)
            tierkey=next((x for x in rdr.fieldnames if x.lower() in ('priority_tier','applicability_tier','similarity_tier')),None)
            if not (tkey and tierkey):
                return False
            d=Counter((r[tkey],r[tierkey]) for r in rdr)
        tiers={t:[d[t,'HIGH_SIMILARITY_SCREEN'],d[t,'MODERATE_SIMILARITY_SCREEN'],d[t,'LOW_SIMILARITY_EXPLORATORY_ONLY']] for t in TARGETS}
    x=np.arange(len(TARGETS)); total=np.array([tiers[t] for t in TARGETS],dtype=float)
    den=np.maximum(total.sum(axis=1),1); pc=100*total/den[:,None]
    fig,ax=plt.subplots(figsize=(9.4,5))
    bottoms=np.zeros(len(TARGETS))
    for i,(name,color) in enumerate([('High','#26734d'),('Moderate','#e5b155'),('Low','#bd6868')]):
        ax.bar(x,pc[:,i],bottom=bottoms,label=name,color=color)
        bottoms+=pc[:,i]
    ax.set_xticks(x,[label(t) for t in TARGETS]); ax.set_ylabel('Analogs (%)')
    ax.set_ylim(0,100); ax.legend(loc='upper right'); ax.set_title('Applicability-domain structural similarity tiers')
    fig.tight_layout(); fig.savefig(OUT/'fig03_applicability_coverage.png',dpi=180);plt.close(fig)
    return True

def make_notebook():
    """Review notebook: renders existing results and charts; NO refitting."""
    md = lambda s: {'cell_type':'markdown','metadata':{},'source':[line+'\n' for line in s.strip().split('\n')]}
    code = lambda s: {'cell_type':'code','metadata':{},'execution_count':None,'outputs':[],
                      'source':[line+'\n' for line in s.strip().split('\n')]}
    cells=[
      md('''# Member 5 — Tryptanthrin anticancer ML results review
This notebook **reviews saved outputs only**. It does not train, tune, or select models.\n\nRun from the project root with the existing `.venv`. For full methodology see `reports/member5_final/Member5_ML_Final_Report.md`.'''),
      code('''from pathlib import Path
import csv
import matplotlib.pyplot as plt
from IPython.display import display, Markdown, Image
root = Path.cwd()
if not (root / 'reports' / 'all_targets_model_comparison.csv').exists():
    root = Path.cwd().parent
reports = root / 'reports'
assert (reports / 'all_targets_model_comparison.csv').exists(), 'Open notebook from project root or notebooks/ folder'
def rows(filename):
    with open(filename, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
'''),
      md('''## Fixed scaffold-split model comparison\nThe reported test split is held out by Bemis–Murcko scaffolds. Scores are preliminary estimates, not external validation.'''),
      code('''results = [r for r in rows(reports / 'all_targets_model_comparison.csv') if r['split'] == 'test']
for r in results:
    print(f"{r['target']:10} {r['model']:14} MAE={float(r['MAE']):.4f} RMSE={float(r['RMSE']):.4f} R2={float(r['R2']):.4f}")
display(Image(filename=str(reports / 'member5_final' / 'fig01_original_rf_vs_xgb.png')))'''),
      md('''## Leakage and fingerprint collision sensitivity\nDistinct molecules can share identical non-chiral folded fingerprints. Collision exclusion is a diagnostic, not an independent test set.'''),
      code('''checks = rows(reports / 'final_validation' / 'leakage_and_dataset_audit.csv')
for r in checks:
    print(r['target'], 'scaffold overlap:', r['overlapping_scaffolds_train_test'],
          '| exact structures:', r['overlapping_canonical_smiles_train_test'],
          '| shared fingerprints:', r['overlapping_fingerprints_train_test'])'''),
      md('''## Alternative scaffold split robustness (Random Forest only)\nThree new RF refits per target (seeds 13, 23, 37). XGBoost was **not** repeated.'''),
      code('''rob = rows(reports / 'final_validation' / 'repeated_scaffold_robustness' / 'rf_robustness_summary.csv')
for r in rob:
    print(r['target'], 'R2 mean:', round(float(r['R2_mean']),4), 'SD:',round(float(r['R2_std']),4))
display(Image(filename=str(reports / 'member5_final' / 'fig02_rf_scaffold_robustness.png')))'''),
      md('''## Analog screening and final handover\nCandidates are in-silico hypotheses; source datasets and assay comparability impose significant limitations.'''),
      code('''handover = rows(reports / 'diversity_shortlist' / 'master_docking_handover.csv')
print('Target–analog entries:',len(handover))
print('Unique Analog_IDs:',len({r['Analog_ID'] for r in handover}))
for t in sorted({r['Target'] for r in handover}):
    print(t, sum(r['Target']==t for r in handover))
if (reports/'member5_final'/'fig03_applicability_coverage.png').exists():
    display(Image(filename=str(reports/'member5_final'/'fig03_applicability_coverage.png')))'''),
      md('''## Conclusions and limitations\n- Regression on target-specific curated median pIC50, not binary class classification.\n- Structural applicability is poor for EGFR, VEGFR2 and PI3Kα generated analog predictions.\n- TDO2 has weaker performance under alternative scaffold splits; similarity support is not predictive confidence.\n- Test results are from one original scaffold split; only Random Forest has a three-repeat robustness assessment.\n- Parent compound predicted pIC50 and reference identities require independent checks; no experimentally validated anticancer hits are claimed.\n- Team should confirm whether the original six-target scope was formally reduced to five.''')
    ]
    doc={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},
                                  'language_info':{'name':'python'}},'nbformat':4,'nbformat_minor':5}
    NOTEBOOKS.mkdir(parents=True,exist_ok=True)
    (NOTEBOOKS/'Member5_results_review.ipynb').write_text(json.dumps(doc,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')

def main():
    required=['model_metrics','leakage','collision','robustness','diversity','handover']
    missing=[str(INPUT[k].relative_to(ROOT)) for k in required if not INPUT[k].is_file()]
    if missing:
        raise SystemExit('Missing required outputs. Run their earlier steps first:\n'+'\n'.join(' - '+s for s in missing))
    data={k:read_csv(p) for k,p in INPUT.items() if p.is_file()}
    if any(not data[k] for k in required):
        raise SystemExit('A required report CSV is empty; cannot generate a verified report.')
    metrics=data['model_metrics']; leak=indexed(data['leakage']); robust=indexed(data['robustness']); handover=data['handover']; div=indexed(data['diversity'], 'Target')
    for t in TARGETS:
        if t not in leak or t not in robust or t not in div:
            raise SystemExit(f'Missing audit/robustness/diversity result for {t}')
        for model in ('Random Forest','XGBoost'):
            pick(metrics,t,model)
    if len(handover) != 50:
        print(f'WARNING: shortlist has {len(handover)} rows, not the previously expected 50.')
    OUT.mkdir(parents=True,exist_ok=True)
    comparison=[]
    for t in TARGETS:
        a=pick(metrics,t,'Random Forest');b=pick(metrics,t,'XGBoost');rr=robust[t]
        comparison.append({'target':t,'test_n':a['n'],'RF_test_MAE':a['MAE'],'RF_test_RMSE':a['RMSE'],
                           'RF_test_R2':a['R2'],'XGB_test_MAE':b['MAE'],'XGB_test_RMSE':b['RMSE'],
                           'XGB_test_R2':b['R2'],'RF_repeat3_mean_R2':rr['R2_mean'],
                           'RF_repeat3_sd_R2':rr['R2_std'], 'RF_repeat3_min_R2':rr['R2_min'],
                           'RF_repeat3_max_R2':rr['R2_max']})
    save_csv(OUT/'final_metrics_summary.csv',comparison)
    make_comparison_chart(metrics); make_robustness_chart(metrics, data['robustness'])
    sim_chart=False
    if 'ad_summary' in data:
        try:
            sim_chart=make_similarity_chart(data['ad_summary'])
        except (KeyError,ValueError,TypeError) as e:
            print('Applicability coverage figure skipped:',e)
    summary=[]
    counts=Counter(r['Target'] for r in handover)
    unique_analogs=len({r['Analog_ID'] for r in handover})
    summary += ['# Member 5 — Final ML Evaluation and Integration Report','',
        f'Generated (UTC): {datetime.now(timezone.utc).isoformat()}',
        '**Project:** Tryptanthrin analog anticancer target-wise bioactivity screening',
        '**Scope:** Five modelled targets (EGFR, VEGFR2, PI3Kα, IDO1, TDO2). The original six-target plan requires team confirmation.','',
        '## 1. Role and objectives','',
        'Member 5 developed target-specific ECFP4/Morgan fingerprint features, trained Random Forest and XGBoost pIC50 regression models, evaluated scaffold-disjoint test performance, audited fingerprints and train/test leakage, measured analog structural applicability, and integrated in-silico candidate predictions.','',
        '## 2. Data and modelling approach','',
        '- Input labels are target-specific curated median **pIC50** from exact IC50 measurements (ChEMBL); they are **not** binary active/inactive classifications. BindingDB source data were not automatically merged into the trained datasets.',
        '- Features: RDKit Morgan radius 2 (ECFP4-like), 2048-bit folded fingerprints; chirality not explicitly encoded.',
        '- Train/validation/test: two-stage Bemis–Murcko scaffold GroupShuffleSplit, random seed 42 for original models; observed molecule proportions depend on scaffold group sizes.',
        '- Random Forest: 200 estimators, min_samples_leaf=2, max_features=sqrt, random_state=42. XGBoost: histogram tree method with validation-based early stopping.',
        '- Evaluation: MAE, RMSE and R² on a fixed scaffold-held-out test partition, plus training-median baseline. Higher pIC50 corresponds to lower IC50 within the defined endpoint, but cross-assay comparability has not been established.','',
        '## 3. Original held-out test performance','',
        '| Target | Test n | RF MAE | RF RMSE | RF R² | XGB MAE | XGB RMSE | XGB R² |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in comparison:
        summary.append(f"| {label(r['target'])} | {r['test_n']} | {fmt(r['RF_test_MAE'])} | {fmt(r['RF_test_RMSE'])} | {fmt(r['RF_test_R2'])} | {fmt(r['XGB_test_MAE'])} | {fmt(r['XGB_test_RMSE'])} | {fmt(r['XGB_test_R2'])} |")
    summary += ['', '![Original fixed-split model comparison](fig01_original_rf_vs_xgb.png)','',
        '**Interpretation:** Test performance is preliminary and represents one scaffold split per target. Differences between RF and XGBoost do not establish a universally superior algorithm. Do not repeatedly use the original test set for hyperparameter tuning.','',
        '## 4. Data integrity and fingerprint leakage audit','',
        '| Target | Training | Test | Shared train/test scaffolds | Shared exact structures | Shared fingerprints |',
        '|---|---:|---:|---:|---:|---:|']
    for t in TARGETS:
        l=leak[t]
        summary.append(f"| {label(t)} | {l['train']} | {l['test']} | {l['overlapping_scaffolds_train_test']} | {l['overlapping_canonical_smiles_train_test']} | {l['overlapping_fingerprints_train_test']} |")
    summary += ['', 'Distinct structures can produce identical 2048-bit Morgan fingerprints, including stereochemical differences or folding collisions. This is representation-level overlap, not proof of exact-molecule duplication. Chemical standardization for salt/tautomer equivalence remains incomplete.', '',
        '## 5. Fingerprint collision sensitivity','',
        'The existing test predictions were re-evaluated after excluding test molecules whose fingerprints matched training fingerprints. The analysis did not retrain models. Results are descriptive and cannot replace independently designed validation.','']
    coll=data['collision']
    if coll:
        summary.append('Source: `reports/final_validation/fingerprint_collision_sensitivity_metrics.csv`.')
    summary += ['', '## 6. Repeated scaffold-split robustness — Random Forest only','',
        'Three new scaffold-disjoint splits (seeds 13, 23, 37) were each used to refit the baseline RF configuration. Variability reflects changed held-out scaffolds and activity distributions. **No analogous XGBoost repeated-split experiment was performed.**','',
        '| Target | Original RF R² | RF repeat mean R² | Sample SD | Range |',
        '|---|---:|---:|---:|---:|']
    for t in TARGETS:
        r=robust[t]; orig=pick(metrics,t,'Random Forest')
        summary.append(f"| {label(t)} | {fmt(orig['R2'])} | {fmt(r['R2_mean'])} | {fmt(r['R2_std'])} | {fmt(r['R2_min'])}–{fmt(r['R2_max'])} |")
    summary += ['', '![Random Forest repeated scaffold split robustness](fig02_rf_scaffold_robustness.png)', '',
        '**TDO2 caveat:** Its relatively small dataset and lower repeated RF performance warrant particular caution, even if generated analogs are often structurally similar to TDO2 training compounds.', '',
        '## 7. Applicability and compound predictions','',
        '- The filtered generated library has 11,355 analogs. Predictions were generated for all five targets, yielding 56,775 target–analog rows. Both RF and XGBoost provide predictions; their agreement is not a calibrated uncertainty estimate.',
        '- Training-only nearest-neighbor Morgan Tanimoto is used as an exploratory structural applicability indicator. Similarity thresholds are heuristic, not validated confidence bounds.',
        '- EGFR, VEGFR2 and PI3Kα generated analog predictions mostly fall in the low-similarity category; IDO1 and TDO2 have better structural coverage. Better coverage does not establish experimental efficacy or robust model extrapolation.', '',
        '![Applicability-domain tier coverage](fig03_applicability_coverage.png)' if sim_chart else 'The applicability-domain tier chart was skipped because matching tier-count columns were unavailable. Original results remain under `reports/applicability_domain/`.', '',
        '## 8. Integrated diversity shortlist','',
        f'- {len(handover)} target–analog shortlist entries correspond to {unique_analogs} distinct Analog IDs. This is a computational docking handover, **not** a set of experimentally validated hits.',
        '- Ranking considers predicted activity, RF/XGB disagreement, applicability tier, drug-likeness and synthetic accessibility; within-target diversity uses pairwise fingerprint Tanimoto ≤ 0.75.',
        '- Chemically diverse selection was restricted to the top 500 eligible rows for IDO1 and TDO2; the final candidate list is dependent on this pre-ranking.',
        '- Parent triptanthrin comparison is structural unless parent predicted pIC50 and relevant assay context are separately verified.', '',
        '| Target | Selected | Maximum pairwise selected Tanimoto |', '|---|---:|---:|']
    for t in TARGETS:
        d=div[t]
        summary.append(f"| {label(t)} | {counts[t]} | {fmt(d['Max_pairwise_Tanimoto'])} |")
    summary += ['', '## 9. Methodological limitations and project decisions','',
        '- No direct in-vitro target inhibition, cellular anticancer activity or external prospective validation is established by these predictions.',
        '- Exact ChEMBL endpoint curation and median label construction are inherited from Member 4. Mixed assay conditions, biological heterogeneity, potential stereochemical conflation, salt/tautomer standardization and outlier records require documented review.',
        '- Scaffold isolation reduces one form of leakage; distinct structures with equal fingerprints remain possible; overlap sensitivity was assessed separately.',
        '- Three repeated RF splits are insufficient to calculate a calibrated generalization confidence interval. The original and repeated results should both be shown.',
        '- Since pIC50 is continuous, classification class-imbalance metrics do not directly apply; examine underrepresented extremes in the response distribution.',
        '- Candidate similarity tiers and RF/XGB model agreement are prioritization heuristics, not calibrated probability or prediction intervals.',
        '- Target selectivity cannot be inferred solely by ranking pIC50 from heterogeneous target-specific datasets; harmonized experimental validation is necessary.',
        '- A sixth target was described in the original team work division but is not included in the current five-target modelling workflow: document the agreed scope.', '',
        '## 10. Reproducibility and handover checklist','',
        '1. Preserve original `data/raw/` target datasets and filtered analog library, along with Member 4/Member 3 version details.',
        '2. Preserve scripts `01`–`09` and `11`–`14`, saved model artifacts and environment package versions.',
        '3. Store `reports/all_targets_model_comparison.csv`, Step 11/12 audit tables, Step 13 repeat-split results and Step 6 applicability results.',
        '4. Share `reports/diversity_shortlist/master_docking_handover.csv` and the independently audited handover with the docking team; verify unique structures and parent/reference identities.',
        '5. Ensure figures include captions stating that all activity results are model-derived estimates; do not call candidates confirmed inhibitors.',
        '6. Open `notebooks/Member5_results_review.ipynb` to inspect final outputs without retraining models.', '',
        '## Source files (relative to project root)','']
    for k,p in INPUT.items():
        summary.append(f"- `{p.relative_to(ROOT)}`" + (' (available)' if p.exists() else ' (optional / unavailable)'))
    (OUT/'Member5_ML_Final_Report.md').write_text('\n'.join(summary)+'\n',encoding='utf-8')
    make_notebook()
    print('Step 14 complete — final Member 5 report compiled from YOUR saved project files.')
    print(' Report : reports/member5_final/Member5_ML_Final_Report.md')
    print(' Table  : reports/member5_final/final_metrics_summary.csv')
    print(' Figures: reports/member5_final/fig01_original_rf_vs_xgb.png')
    print('          reports/member5_final/fig02_rf_scaffold_robustness.png')
    if sim_chart:print('          reports/member5_final/fig03_applicability_coverage.png')
    print(' Notebook: notebooks/Member5_results_review.ipynb')
    print('No models, training datasets, predictions or shortlist files modified.')

if __name__=='__main__':
    main()
