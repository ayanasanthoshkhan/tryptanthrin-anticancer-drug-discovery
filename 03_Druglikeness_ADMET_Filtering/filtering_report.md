# Drug-Likeness & ADMET Filtering Report (Member 3 – Abhi)

**Input:** `all_analogs_with_SA.csv` – 11,393 unique tryptanthrin-type analogs (Analog_ID and Category kept unchanged).
**Tools:** RDKit 2026.09.1 (descriptors, QED, FilterCatalog), pandas. All 11,393 SMILES parsed successfully.

## 1. Descriptors calculated
MW, LogP (Crippen), TPSA, HBD, HBA, QED, rotatable bonds, plus heavy atoms, ring counts, aromatic rings, Fsp3. SA scores from Nakshathra's file were carried through (all < 3, so no SA-based removal).

| Property | Mean | Min | Max |
|---|---|---|---|
| MW | 333.0 | 248.2 | 464.5 |
| LogP | 2.25 | 0.13 | 4.18 |
| TPSA | 74.8 | 52.0 | 138.2 |
| HBD | 0.59 | 0 | 2 |
| HBA | 4.49 | 3 | 7 |
| QED | 0.54 | 0.27 | 0.63 |
| RotB | 1.29 | 0 | 4 |

## 2. Rules and filters applied
- **Lipinski Ro5:** MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 (max 1 violation allowed).
- **Veber:** RotB ≤ 10 and TPSA ≤ 140.
- **Structural alerts:** PAINS (A/B/C), Brenk (unwanted/reactive/toxic groups), NIH.
- **Rule-based ADMET proxies (not ML):** Egan (GI absorption: TPSA ≤ 131.6, LogP ≤ 5.88), BBB heuristic (TPSA ≤ 90, MW ≤ 450, HBD ≤ 3), Pfizer 3/75 (LogP > 3 and TPSA < 75 = higher in-vivo toxicity risk), GSK 4/400, solubility risk (LogP > 4 or ≥ 4 aromatic rings).

## 3. Decision logic
- **EXCLUDE:** PAINS alert, or ≥ 2 Lipinski violations, or Veber failure, or QED < 0.20, or unparsable SMILES.
- **FLAG (kept, with reason):** Brenk/NIH alert, 1 Lipinski violation, QED < 0.35, low predicted GI absorption, Pfizer 3/75, solubility risk.
- **RETAIN:** passes everything.

## 4. Results
| Decision | Count | % |
|---|---|---|
| RETAIN | 8,369 | 73.5% |
| FLAG | 2,986 | 26.2% |
| EXCLUDE | 38 | 0.3% |

- **Lipinski:** all 11,393 have 0 violations; none failed Veber. The library is uniformly small and rule-of-five compliant, so these rules do not discriminate here.
- **PAINS (38 excluded):** 32 × `anil_di_alk_C` (dialkylamino-aniline, redox-active) and 6 × `catechol_A` (also Brenk catechol). All 38 are Exploratory_design. No Parent or Published_reference compound was excluded.
- **Brenk alerts (1,163 compounds):** aniline 1,101, catechol 42, N–O single bond 17, nitro 17, hydroquinone 14, iodine 7. Aniline is the amine substituent on the core; it is a metabolic-activation concern, not a hard failure, so these are flagged only.
- **Pfizer 3/75:** 1,803 flagged (LogP > 3 with TPSA < 75).
- **Egan:** 29 predicted low GI absorption (TPSA > 131.6).
- **Solubility risk:** 84. **Low QED (< 0.35):** 2. **NIH alerts:** 0.
- **Lead-like** (MW 200–350, LogP ≤ 3.5, RotB ≤ 7): 7,905. **BBB likely:** 8,902.
- **Parent (tryptanthrin):** RETAIN. **Published references (37):** 32 retained, 5 flagged (iodine, aniline, nitro/N–O, Pfizer 3/75); none excluded.

## 5. Why molecules are retained or flagged
Every row in `descriptor_table.csv` has a `Decision` and a plain-text `Reason` column listing each rule that triggered, or "Passes all filters".

## 6. Limitations
- ADMET here is rule-based; it is not a trained model and gives no quantitative hERG, Ames or CYP values. Treat flags as prioritisation hints.
- PAINS/Brenk matches are alerts, not proof of toxicity or assay interference.
- Filters are not applied to SA (already screened upstream).

## 7. Deliverables
1. `descriptor_table.csv`: all 11,393 compounds with descriptors, filters, ADMET proxies, decision and reason.
2. `filtered_analog_library.csv`: 11,355 compounds (RETAIN + FLAG).
3. `strict_clean_analogs.csv`: 8,369 compounds with no flags (RETAIN only).
4. `excluded_analogs_with_reasons.csv` (38) and `flagged_analogs_with_reasons.csv` (2,986).
5. This report.
