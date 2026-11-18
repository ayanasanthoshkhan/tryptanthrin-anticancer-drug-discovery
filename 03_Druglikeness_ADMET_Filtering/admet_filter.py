import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")
IN = "all_analogs_with_SA.csv"   # upload this to Colab first
OUT = ""

df = pd.read_csv(IN)

def catalog(*cats):
    p = FilterCatalogParams()
    for c in cats:
        p.AddCatalog(c)
    return FilterCatalog(p)

PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
BRENK = catalog(FilterCatalogParams.FilterCatalogs.BRENK)
NIH = catalog(FilterCatalogParams.FilterCatalogs.NIH)

def alerts(cat, m):
    return "; ".join(sorted({e.GetDescription() for e in cat.GetMatches(m)}))

rows = []
for smi in df["SMILES"]:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        rows.append(None)
        continue
    rows.append(dict(
        MW=Descriptors.MolWt(m),
        LogP=Crippen.MolLogP(m),
        TPSA=rdMolDescriptors.CalcTPSA(m),
        HBD=Lipinski.NumHDonors(m),
        HBA=Lipinski.NumHAcceptors(m),
        QED=QED.qed(m),
        RotB=rdMolDescriptors.CalcNumRotatableBonds(m),
        HeavyAtoms=m.GetNumHeavyAtoms(),
        Rings=rdMolDescriptors.CalcNumRings(m),
        AromRings=rdMolDescriptors.CalcNumAromaticRings(m),
        FractionCSP3=rdMolDescriptors.CalcFractionCSP3(m),
        MR=Crippen.MolMR(m),
        PAINS_alerts=alerts(PAINS, m),
        Brenk_alerts=alerts(BRENK, m),
        NIH_alerts=alerts(NIH, m),
    ))

invalid = [i for i, r in enumerate(rows) if r is None]
desc = pd.DataFrame([r if r else {} for r in rows])
out = pd.concat([df, desc], axis=1)
out["RDKit_parse_ok"] = [r is not None for r in rows]

# ---- Lipinski / Veber / Egan / lead-likeness ----
out["Lipinski_violations"] = ((out.MW > 500).astype(int) + (out.LogP > 5).astype(int)
                              + (out.HBD > 5).astype(int) + (out.HBA > 10).astype(int))
out["Lipinski_pass"] = out.Lipinski_violations <= 1          # standard "Ro5: max 1 violation"
out["Veber_pass"] = (out.RotB <= 10) & (out.TPSA <= 140)
out["Egan_pass"] = (out.LogP <= 5.88) & (out.TPSA <= 131.6)
out["Lead_like"] = out.MW.between(200, 350) & (out.LogP <= 3.5) & (out.RotB <= 7)

# ---- Rule-based ADMET proxies (no ML model: transparent, reproducible) ----
out["GI_absorption_pred"] = np.where(out.Egan_pass, "High", "Low")
out["BBB_permeant_pred"] = np.where((out.TPSA <= 90) & (out.MW <= 450) & (out.HBD <= 3), "Likely", "Unlikely")
out["Pfizer_3_75_risk"] = (out.LogP > 3) & (out.TPSA < 75)   # higher in-vivo tox / promiscuity risk
out["GSK_4_400_pass"] = (out.MW <= 400) & (out.LogP <= 4)
out["Solubility_risk"] = (out.LogP > 4) | (out.AromRings >= 4)

# ---- Structural / toxicity alerts ----
out["PAINS_flag"] = out.PAINS_alerts.fillna("") != ""
out["Brenk_flag"] = out.Brenk_alerts.fillna("") != ""
out["NIH_flag"] = out.NIH_alerts.fillna("") != ""
out["Tox_alert_count"] = (out.Brenk_alerts.fillna("").apply(lambda s: len([x for x in s.split("; ") if x]))
                         + out.NIH_alerts.fillna("").apply(lambda s: len([x for x in s.split("; ") if x])))

# ---- Decision logic: HARD EXCLUDE / FLAG / RETAIN ----
def decide(r):
    hard, soft = [], []
    if not r.RDKit_parse_ok:
        return "EXCLUDE", "RDKit could not parse SMILES"
    if r.PAINS_flag:
        hard.append("PAINS alert: " + r.PAINS_alerts)
    if r.Lipinski_violations >= 2:
        hard.append(f"Lipinski: {int(r.Lipinski_violations)} violations")
    if not r.Veber_pass:
        hard.append("Veber fail (RotB>10 or TPSA>140)")
    if r.QED < 0.20:
        hard.append(f"Very low QED ({r.QED:.2f})")
    if r.NIH_flag:
        soft.append("NIH alert: " + r.NIH_alerts)
    if r.Brenk_flag:
        soft.append("Brenk alert: " + r.Brenk_alerts)
    if r.Lipinski_violations == 1:
        soft.append("Lipinski: 1 violation")
    if r.QED < 0.35:
        soft.append(f"Low QED ({r.QED:.2f})")
    if not r.Egan_pass:
        soft.append("Predicted low GI absorption (Egan)")
    if r.Pfizer_3_75_risk:
        soft.append("Pfizer 3/75 (LogP>3 & TPSA<75)")
    if r.Solubility_risk:
        soft.append("Solubility risk (LogP>4 or >=4 aromatic rings)")
    if hard:
        return "EXCLUDE", " | ".join(hard + soft)
    if soft:
        return "FLAG", " | ".join(soft)
    return "RETAIN", "Passes all filters"

res = out.apply(decide, axis=1, result_type="expand")
out["Decision"], out["Reason"] = res[0], res[1]

# ---- Save ----
cols_front = ["Analog_ID", "Category", "SMILES"]
out.to_csv(OUT + "descriptor_table_full.csv", index=False)
keep = out[out.Decision != "EXCLUDE"]
keep.to_csv(OUT + "filtered_analog_library.csv", index=False)
out[out.Decision == "EXCLUDE"].to_csv(OUT + "excluded_analogs_with_reasons.csv", index=False)
out[out.Decision == "FLAG"].to_csv(OUT + "flagged_analogs_with_reasons.csv", index=False)

# ---- Report stats ----
n = len(out)
stats = dict(
    n=n, invalid=len(invalid),
    retain=(out.Decision == "RETAIN").sum(), flag=(out.Decision == "FLAG").sum(),
    exclude=(out.Decision == "EXCLUDE").sum(),
    lip0=(out.Lipinski_violations == 0).sum(), lip1=(out.Lipinski_violations == 1).sum(),
    lip2=(out.Lipinski_violations >= 2).sum(),
    veber=(~out.Veber_pass).sum(), egan=(~out.Egan_pass).sum(),
    pains=out.PAINS_flag.sum(), brenk=out.Brenk_flag.sum(), nih=out.NIH_flag.sum(),
    lowqed=(out.QED < 0.35).sum(), pfizer=out.Pfizer_3_75_risk.sum(),
    bbb=(out.BBB_permeant_pred == "Likely").sum(),
)
print(stats)
print(out.groupby(["Category", "Decision"]).size().unstack(fill_value=0))
print(out[["MW", "LogP", "TPSA", "HBD", "HBA", "QED", "RotB"]].describe().round(2))

from collections import Counter
for c in ["PAINS_alerts", "Brenk_alerts", "NIH_alerts"]:
    cnt = Counter(x for s in out[c].fillna("") for x in s.split("; ") if x)
    print(c, cnt.most_common(12))

print("Parent / Published_reference decisions:")
print(out[out.Category.isin(["Parent", "Published_reference"])][["Analog_ID", "Category", "Decision", "Reason"]]
      .to_string(max_colwidth=90))
out.to_pickle("out.pkl")
