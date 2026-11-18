"""Audit Member 4's ChEMBL IC50 datasets before ML. No values are removed.
Run from the project root: python scripts/01_audit_datasets.py
"""
from pathlib import Path
from collections import Counter,defaultdict
import csv, math
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
REPORT = ROOT / 'reports'
REPORT.mkdir(exist_ok=True)
EXPECTED = ['EGFR', 'VEGFR2', 'PI3Kalpha', 'IDO1', 'TDO2']
FIELDNAMES = ['target','rows','missing_smiles','invalid_smiles','duplicate_input_smiles','duplicate_canonical_structures','missing_pic50','invalid_ic50','pic50_mismatch_gt_0_01','pic50_outside_3_to_12','measurements_gt_1','pic50_std_missing','min_pic50','median_pic50','max_pic50']
summary=[]
flags=[]
for target in EXPECTED:
    filename=RAW / f'{target}_Final_ML_Dataset.csv'
    if not filename.exists():
        raise FileNotFoundError(f'Missing dataset: {filename}')
    with filename.open(newline='', encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        required={'canonical_smiles','IC50_nM','pIC50','molecule_chembl_id','measurement_count'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{filename}: missing columns {required-set(reader.fieldnames or [])}')
        rows=list(reader)
    tally=Counter(); seen_raw=Counter(); seen_canon=Counter(); y=[]
    for line,row in enumerate(rows, start=2):
        smiles=(row.get('canonical_smiles') or '').strip()
        identity=row.get('molecule_chembl_id') or f'line_{line}'
        if not smiles:
            tally['missing_smiles']+=1;flags.append([target,identity,'missing_smiles',''])
        else:
            seen_raw[smiles]+=1
            mol=Chem.MolFromSmiles(smiles)
            if mol is None:
                tally['invalid_smiles']+=1; flags.append([target,identity,'invalid_smiles',smiles])
            else:
                seen_canon[Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)]+=1
        try:
            ic50=float(row['IC50_nM'])
            if not math.isfinite(ic50) or ic50<=0: raise ValueError()
        except (ValueError,TypeError):
            ic50=None;tally['invalid_ic50']+=1;flags.append([target,identity,'invalid_ic50',str(row.get('IC50_nM'))])
        try:
            pic50=float(row['pIC50'])
            if not math.isfinite(pic50): raise ValueError()
            y.append(pic50)
        except (ValueError,TypeError):
            pic50=None;tally['missing_pic50']+=1;flags.append([target,identity,'invalid_or_missing_pic50',str(row.get('pIC50'))])
        if pic50 is not None:
            if pic50<3 or pic50>12:
                tally['pic50_outside_3_to_12']+=1
                flags.append([target,identity,'review_extreme_pic50',f'{pic50:.6g}'])
            if ic50 is not None:
                expected=9-math.log10(ic50)
                if abs(pic50-expected)>0.01:
                    tally['pic50_mismatch_gt_0_01']+=1
                    flags.append([target,identity,'pic50_ic50_mismatch',f'actual={pic50:.5f},expected={expected:.5f}'])
        try:
            if int(float(row.get('measurement_count') or '0'))>1: tally['measurements_gt_1']+=1
        except ValueError: flags.append([target,identity,'invalid_measurement_count',str(row.get('measurement_count'))])
        if not (row.get('pIC50_std') or '').strip(): tally['pic50_std_missing']+=1
    y.sort()
    summary.append({'target':target,'rows':len(rows),'missing_smiles':tally['missing_smiles'],
        'invalid_smiles':tally['invalid_smiles'],
        'duplicate_input_smiles':sum(n-1 for n in seen_raw.values() if n>1),
        'duplicate_canonical_structures':sum(n-1 for n in seen_canon.values() if n>1),
        'missing_pic50':tally['missing_pic50'],'invalid_ic50':tally['invalid_ic50'],
        'pic50_mismatch_gt_0_01':tally['pic50_mismatch_gt_0_01'],
        'pic50_outside_3_to_12':tally['pic50_outside_3_to_12'],
        'measurements_gt_1':tally['measurements_gt_1'],'pic50_std_missing':tally['pic50_std_missing'],
        'min_pic50':round(y[0],4) if y else '',
        'median_pic50':round(y[len(y)//2],4) if y else '',
        'max_pic50':round(y[-1],4) if y else ''})
    print(f"{target:10s} rows={len(rows):5d} | invalid SMILES={tally['invalid_smiles']} | outlier flags={tally['pic50_outside_3_to_12']}")
with (REPORT/'dataset_audit_summary.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=FIELDNAMES);w.writeheader();w.writerows(summary)
with (REPORT/'rows_for_review.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['target','molecule_chembl_id','issue','details']);w.writerows(flags)
print('Audit complete. Check reports/dataset_audit_summary.csv and reports/rows_for_review.csv')
print('Outliers are REVIEW flags only, never automatically deleted.')
