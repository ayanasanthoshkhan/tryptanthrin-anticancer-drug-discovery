"""Step 9 — diversity-aware, target-specific triptanthrin shortlist.

Run in existing project:
    .\\.venv\\Scripts\\python.exe scripts/09_diversity_and_docking_shortlist.py

Inputs (from prior steps):
    data/analogs/filtered_analog_library.csv
    reports/candidate_prioritization/{TARGET}_screen_candidates.csv

Outputs:
    reports/diversity_shortlist/{TARGET}_diverse_top10.csv
    reports/diversity_shortlist/master_docking_handover.csv
    reports/diversity_shortlist/parent_comparison.csv
    reports/diversity_shortlist/diversity_summary.csv
    reports/diversity_shortlist/methodology.json

This is computational triage, NOT verified binding, potency, novelty or synthesis.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "data" / "analogs" / "filtered_analog_library.csv"
RANKINGS = ROOT / "reports" / "candidate_prioritization"
OUT = ROOT / "reports" / "diversity_shortlist"
TARGETS = ("EGFR", "VEGFR2", "PI3Kalpha", "IDO1", "TDO2")
N_PER_TARGET = 10
# Diversity threshold is a heuristic, not a calibrated definition of chemical novelty.
MAX_PAIRWISE_SIMILARITY = 0.75
# Inspect a substantial candidate pool, not only the original top 25.
MAX_POOL = 500
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def read_csv(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}. Complete Step 8 and ensure the filtered library is present.")
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"Empty CSV: {path}")
    return rows


def save_csv(path, rows, columns):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def require(rows, keys, label):
    missing = set(keys) - set(rows[0])
    if missing:
        raise ValueError(f"{label}: missing columns: {sorted(missing)}")


def number(row, key):
    try:
        return float(row[key])
    except (ValueError, TypeError, KeyError) as e:
        raise ValueError(f"{row.get('Analog_ID', '?')}: missing/invalid {key}") from e


def mol_fp(smiles, mol_id):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {mol_id}")
    return FP_GEN.GetFingerprint(mol)


def main():
    library_rows = read_csv(LIBRARY)
    require(library_rows, ["Analog_ID", "SMILES", "Category", "Decision", "Reason",
                           "QED", "SA_score", "MW", "LogP", "TPSA", "PAINS_flag",
                           "Brenk_flag", "InChIKey"], "analog library")
    library = {r["Analog_ID"]: r for r in library_rows}
    if len(library) != len(library_rows):
        raise ValueError("Duplicate Analog_ID in library")
    parents = [r for r in library_rows if r["Category"].strip().lower() == "parent"]
    if len(parents) != 1:
        raise ValueError(f"Expected exactly one parent triptanthrin row, found {len(parents)}")
    parent = parents[0]
    parent_id = parent["Analog_ID"]
    fp_cache = {}
    def fp_for(row):
        aid = row["Analog_ID"]
        if aid not in fp_cache:
            fp_cache[aid] = mol_fp(row["SMILES"], aid)
        return fp_cache[aid]
    parent_fp = fp_for(parent)
    OUT.mkdir(parents=True, exist_ok=True)

    master, parent_comparison, summaries = [], [], []
    output_columns = [
        "Target", "Selection_rank", "Original_target_rank", "Analog_ID", "SMILES", "InChIKey",
        "Category", "Decision", "Filter_reason", "Consensus_pIC50", "Pred_RF_pIC50",
        "Pred_XGB_pIC50", "Model_disagreement_abs", "Priority_tier", "Max_train_Tanimoto",
        "Nearest_train_chembl_id", "Exact_training_structure_match", "QED", "SA_score",
        "MW", "LogP", "TPSA", "PAINS_flag", "Brenk_flag",
        "Similarity_to_parent", "Max_similarity_to_earlier_selected", "Diversity_rule",
    ]
    for target in TARGETS:
        candidates = read_csv(RANKINGS / f"{target}_screen_candidates.csv")
        require(candidates, ["Rank", "Analog_ID", "Target", "SMILES", "Screen_eligible",
                             "Consensus_pIC50", "Pred_RF_pIC50", "Pred_XGB_pIC50",
                             "Model_disagreement_abs", "Priority_tier", "Max_train_Tanimoto",
                             "Nearest_train_chembl_id"], f"{target} ranking")
        candidates.sort(key=lambda r: int(r["Rank"]))
        selected = []
        excluded_exact = 0
        missing = 0
        # Only Step 8-screened molecules. Exclude exact training matches from new-lead list.
        pool = []
        for r in candidates:
            if r["Target"] != target or r["Screen_eligible"] != "True":
                raise ValueError(f"Unexpected candidate row for {target}: {r.get('Analog_ID')}")
            aid = r["Analog_ID"]
            if aid not in library:
                missing += 1
                continue
            if r["SMILES"] != library[aid]["SMILES"]:
                raise ValueError(f"SMILES mismatch for {aid}")
            if number(r, "Max_train_Tanimoto") >= 0.999999:
                excluded_exact += 1
                continue
            pool.append(r)
            if len(pool) >= MAX_POOL:
                break
        if missing:
            raise ValueError(f"{target}: {missing} candidate IDs missing in analog library")

        # Greedy MaxMin-style diversity selection: among the early ranked pool,
        # choose the best-ranked molecule that satisfies pairwise diversity.
        for r in pool:
            if len(selected) >= N_PER_TARGET:
                break
            fp = fp_for(library[r["Analog_ID"]])
            sims = [DataStructs.TanimotoSimilarity(fp, fp_for(library[s["Analog_ID"]]))
                    for s in selected]
            if sims and max(sims) > MAX_PAIRWISE_SIMILARITY:
                continue
            selected.append(r)

        # Do not silently relax diversity criteria to fill 10 slots.
        if len(selected) < N_PER_TARGET:
            print(f"  {target}: only {len(selected)} candidates met the strict diversity rule; "
                  "retaining a shorter list.", flush=True)

        rows = []
        for i, r in enumerate(selected, 1):
            aid = r["Analog_ID"]
            lib = library[aid]
            fp = fp_for(lib)
            previous = selected[:i-1]
            max_prev = max((DataStructs.TanimotoSimilarity(
                fp, fp_for(library[p["Analog_ID"]])) for p in previous), default=None)
            parent_sim = DataStructs.TanimotoSimilarity(fp, parent_fp)
            entry = {
                "Target": target, "Selection_rank": i, "Original_target_rank": r["Rank"],
                "Analog_ID": aid, "SMILES": lib["SMILES"], "InChIKey": lib["InChIKey"],
                "Category": lib["Category"], "Decision": lib["Decision"],
                "Filter_reason": lib["Reason"], "Consensus_pIC50": r["Consensus_pIC50"],
                "Pred_RF_pIC50": r["Pred_RF_pIC50"],
                "Pred_XGB_pIC50": r["Pred_XGB_pIC50"],
                "Model_disagreement_abs": r["Model_disagreement_abs"],
                "Priority_tier": r["Priority_tier"],
                "Max_train_Tanimoto": r["Max_train_Tanimoto"],
                "Nearest_train_chembl_id": r["Nearest_train_chembl_id"],
                "Exact_training_structure_match": "False", "QED": lib["QED"],
                "SA_score": lib["SA_score"], "MW": lib["MW"], "LogP": lib["LogP"],
                "TPSA": lib["TPSA"], "PAINS_flag": lib["PAINS_flag"],
                "Brenk_flag": lib["Brenk_flag"],
                "Similarity_to_parent": round(parent_sim, 4),
                "Max_similarity_to_earlier_selected": "" if max_prev is None else round(max_prev, 4),
                "Diversity_rule": f"pairwise_Tanimoto<={MAX_PAIRWISE_SIMILARITY}",
            }
            rows.append(entry)
            master.append(entry)
            parent_comparison.append({"Target": target, "Analog_ID": aid,
                                      "Category": lib["Category"],
                                      "Similarity_to_parent": round(parent_sim, 4),
                                      "Original_target_rank": r["Rank"],
                                      "Consensus_pIC50": r["Consensus_pIC50"]})
        save_csv(OUT / f"{target}_diverse_top10.csv", rows, output_columns)
        pairwise = [DataStructs.TanimotoSimilarity(
            fp_for(library[a["Analog_ID"]]), fp_for(library[b["Analog_ID"]]))
            for i, a in enumerate(selected) for b in selected[i+1:]]
        summary = {
            "Target": target, "Eligible_rows_step8": len(candidates),
            "Nonexact_pool_reviewed": len(pool), "Exact_training_matches_skipped": excluded_exact,
            "Selected": len(selected), "Mean_pairwise_Tanimoto": round(sum(pairwise)/len(pairwise), 4) if pairwise else "",
            "Max_pairwise_Tanimoto": round(max(pairwise), 4) if pairwise else "",
            "Parent_analog_id": parent_id,
            "RETAIN_selected": sum(r["Decision"] == "RETAIN" for r in rows),
            "FLAG_selected": sum(r["Decision"] == "FLAG" for r in rows),
            "Published_reference_selected": sum(r["Category"] == "Published_reference" for r in rows),
        }
        summaries.append(summary)
        print(f"{target:10s} eligible={len(candidates):5,d} | pool={len(pool):3d} | "
              f"selected={len(rows):2d} | max pairwise={summary['Max_pairwise_Tanimoto']} | "
              f"retained={summary['RETAIN_selected']} flagged={summary['FLAG_selected']}", flush=True)

    save_csv(OUT / "master_docking_handover.csv", master, output_columns)
    save_csv(OUT / "parent_comparison.csv", parent_comparison,
             ["Target", "Analog_ID", "Category", "Similarity_to_parent", "Original_target_rank", "Consensus_pIC50"])
    save_csv(OUT / "diversity_summary.csv", summaries, list(summaries[0]))
    methodology = {
        "step": 9, "software": f"RDKit {rdBase.rdkitVersion}",
        "input": "Step 8 eligible candidates + original filtered library",
        "parent_analog_id": parent_id, "parent_smiles": parent["SMILES"],
        "fingerprints": "Morgan/ECFP4, radius 2, 2048-bit, Tanimoto similarity",
        "selection": f"Top-ranked Step 8 screened candidates, inspect first {MAX_POOL} nonexact rows; "
                     f"greedy top-rank-first pairwise similarity <= {MAX_PAIRWISE_SIMILARITY}",
        "target_shortlist_size": N_PER_TARGET,
        "exact_training_match": "Nearest training fingerprint similarity >= 0.999999 excluded from prospective shortlist; "
                                "fingerprint identity is NOT proof of identical chemical identity or novelty",
        "known_compounds": "Published_reference and other categories retained as labels. "
                           "An analog without exact training fingerprint match is NOT necessarily novel.",
        "parent_comparison": "Structural similarity only; parent activity NOT inferred; no potency improvement claimed",
        "unverified": ["No docking", "No experimental target inhibition", "No validated confidence intervals",
                       "No AiZynthFinder retrosynthesis", "No structure-standardization guarantee",
                       "No independent prospective validation"],
        "notes": ["Top-ranked greedy diversity triage, not a globally optimized medoid set",
                  "Results should be manually inspected for chemical feasibility and assay scope",
                  "Cross-target pIC50 model outputs are not proven selectivity or comparable potency",
                  "Retain FLAG, PAINS and Brenk information for medicinal-chemistry review",
                  "Scaffold split and training-neighbor similarity are initial safeguards, not calibrated applicability domain"],
    }
    (OUT / "methodology.json").write_text(json.dumps(methodology, indent=2), encoding="utf-8")
    print(f"\nStep 9 complete. {len(master)} selected target–analog rows (duplicates across targets possible).")
    print(f"Outputs: {OUT.relative_to(ROOT)}")
    print("Shortlists are hypotheses for review and docking, not validated anticancer hits.")

if __name__ == "__main__":
    main()
