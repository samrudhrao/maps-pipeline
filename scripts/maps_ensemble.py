#!/usr/bin/env python3
"""
MAPS Ensemble Prediction - maps_ensemble.py
Step 14b of the MAPS pipeline — RF + XGBoost weighted ensemble prediction.
Combines two AI models for more robust pathogen scoring with SHAP explainability.

Model weights: RF 55% + XGBoost 45%
Composite score = (ensemble probability * 0.40) + (abundance * 0.35) + (pathogen score * 0.25)
Confidence tiers: HIGH >=70%, MODERATE 40-69%, LOW <40%
SHAP values explain top 3 feature drivers for each organism prediction.

Supplemental DB logic (same as full_pipeline.sh):
  - S. cohnii: only added if main DB has <500 reads of known pathogenic Staph
  - B. licheniformis: always added if >=1% abundance

To retrain models: run the model training script in /maps/models/
Model files: maps_ai_model.pkl (RF), maps_xgb_model.pkl (XGB),
             maps_shap_explainer.pkl (RF SHAP), maps_xgb_shap_explainer.pkl (XGB SHAP)

Usage: python3 maps_ensemble.py SAMPLE OUTDIR HUMAN_PCT SPECIMEN_TYPE
Output: OUTDIR/SAMPLE_ensemble_scores.csv
"""

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.expanduser('/maps/scripts'))
try:
    from maps_specimen_scores import get_specimen_adjusted_score, detect_specimen_type
except: 
    get_specimen_adjusted_score = lambda s,o,t: s
    detect_specimen_type = lambda s: "unknown"

import numpy as np
import joblib
import os
import sys

# Load all model files
# RF and XGB models for ensemble prediction
# SHAP explainers for feature importance interpretation
RF_MODEL  = joblib.load(os.path.expanduser("/maps/models/maps_ai_model.pkl"))
XGB_MODEL = joblib.load(os.path.expanduser("/maps/models/maps_xgb_model.pkl"))
LE        = joblib.load(os.path.expanduser("/maps/models/maps_label_encoder.pkl"))
FEATURES  = joblib.load(os.path.expanduser("/maps/models/maps_features.pkl"))
RF_SHAP   = joblib.load(os.path.expanduser("/maps/models/maps_shap_explainer.pkl"))
XGB_SHAP  = joblib.load(os.path.expanduser("/maps/models/maps_xgb_shap_explainer.pkl"))


with open(os.path.expanduser("/maps/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]
def apply_confidence_naming(organism, ai_score):
    """
    >= 70% confidence: report full species name
    < 70% confidence: report genus level only
    """
    if ai_score >= 0.70:
        return organism
    # Extract prefix and genus only
    prefix = ""
    name = organism
    if organism.startswith("[BACT] ") or organism.startswith("[FUNG] "):
        prefix = organism[:7]
        name = organism[7:]
    parts = name.strip().split()
    if len(parts) >= 2:
        return f"{prefix}{parts[0]} sp."
    return organism



PATHOGEN_SCORES = {
    "pseudomonas aeruginosa":1.0,"staphylococcus aureus":1.0,
    "escherichia coli":0.9,"klebsiella pneumoniae":1.0,
    "candida albicans":1.0,"aspergillus fumigatus":0.9,
    "pneumocystis jirovecii":1.0,"histoplasma":1.0,
    "coccidioides":1.0,"mucor mucedo":1.0,"scedosporium":1.0,
    "lomentospora":1.0,"bacillus licheniformis":0.7,
    "microbacterium":0.6,"enterococcus faecium":0.9,
    "citrobacter freundii":0.8,"staphylococcus cohnii":0.8,
    "staphylococcus epidermidis":0.4,"none":0.0,
}

# Ensemble weights and abundance threshold
# RF gets slightly more weight (55%) as it was validated first
# Change weights here if one model consistently outperforms the other
RF_WEIGHT=0.55; XGB_WEIGHT=0.45; THRESHOLD=3.0

def get_pathogen_score(org):
    org=str(org).lower()
    for k,v in PATHOGEN_SCORES.items():
        if k in org: return v
    return 0.4

def get_is_fungal(org):
    return 1 if any(f in org.lower() for f in
        ["aspergillus","candida","fusarium","pneumocystis",
         "histoplasma","mucor","scedosporium","coccidioides","lomentospora"]) else 0

def parse_report(path, min_reads=50):
    rows = {}
    if not os.path.exists(path): return rows
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            pct, reads = float(parts[0]), int(parts[1])
            rank, name = parts[3].strip(), parts[5].strip()
            if rank == "S" and reads >= min_reads:
                if not any(e in name.lower() for e in EXCLUDE):
                    rows[name] = (pct, reads)
    return rows

def get_shap_explanation(X_df):
    """
    Calculate SHAP feature importance for a single prediction.
    Averages RF and XGB SHAP values and returns top 3 drivers.
    Used to explain why the model scored an organism high or low.
    Returns string like 'abundance(0.789) | above_threshold(0.226)'
    """
    try:
        rf_vals  = RF_SHAP.shap_values(X_df)
        xgb_vals = XGB_SHAP.shap_values(X_df)

        # RF SHAP returns different shapes depending on version
        # Handle both (1, 11, 2) and flat array formats
        # Take class 1 (positive class) values
        if isinstance(rf_vals, list):
            rf_vals = rf_vals[1]
        elif hasattr(rf_vals, 'shape') and len(rf_vals.shape) == 3:
            rf_vals = rf_vals[0, :, 1]  # shape (1,11,2) -> (11,)
        else:
            rf_vals = rf_vals[0]

        # XGB SHAP returns (1, 11) shape — flatten to 1D array
        if hasattr(xgb_vals, 'shape') and len(xgb_vals.shape) == 2:
            xgb_vals = xgb_vals[0]
        else:
            xgb_vals = xgb_vals[0]

        avg_shap = [(FEATURES[i], (abs(rf_vals[i]) + abs(xgb_vals[i]))/2)
                    for i in range(len(FEATURES))]
        avg_shap = sorted(avg_shap, key=lambda x: -x[1])
        top3 = [f"{f}({v:.3f})" for f,v in avg_shap[:3] if v > 0.001]
        return " | ".join(top3) if top3 else "N/A"
    except Exception as e:
        return f"N/A"

def predict_ensemble(sample, outdir, human_pct=20.0):
    """
    Run RF + XGBoost ensemble prediction on all detected organisms.
    Reads from Kraken2 reports and merges supplemental DB results.
    Returns sorted list of organism prediction dictionaries with SHAP explanations.
    """
    bact = parse_report(f"{outdir}/bact_results/{sample}_bact_report.txt")
    fung = parse_report(f"{outdir}/fungi_results/{sample}_fungi_report.txt")
    supp = parse_report(f"{outdir}/bact_results/{sample}_extra_report.txt")

    # Supplemental DB merge — same logic as full_pipeline.sh
    # Only S. cohnii and B. licheniformis from supplemental DB
    SUPP_KEEP = ["staphylococcus cohnii", "bacillus licheniformis"]

    # Check main DB for dominant known pathogenic Staphylococcus
    KNOWN_STAPH = ['aureus', 'epidermidis', 'lugdunensis', 'saprophyticus',
                   'haemolyticus', 'capitis', 'warneri', 'caprae',
                   'hominis', 'pasteuri', 'pettenkoferi', 'simulans']
    main_staph_reads = sum(v[1] for k,v in bact.items()
                          if 'staphylococcus' in k.lower()
                          and any(s in k.lower() for s in KNOWN_STAPH))

    for org, (pct, reads) in supp.items():
        if not any(s in org.lower() for s in SUPP_KEEP):
            continue
        if 'cohnii' in org.lower() and main_staph_reads > 500:
            continue  # Main DB has dominant known Staph — skip S. cohnii
        if pct >= 1.0:
            if org not in bact:
                bact[org] = (pct, reads)
            elif pct > bact[org][0] * 2:
                bact[org] = (pct, reads)

    bact_df = pd.Series({k:v[0] for k,v in bact.items()})
    fung_df = pd.Series({k:v[0] for k,v in fung.items()})
    if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
    if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
    bact_df.index = ["[BACT] " + x for x in bact_df.index]
    fung_df.index = ["[FUNG] " + x for x in fung_df.index]
    combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)

    total_reads     = sum(v[1] for v in list(bact.values())) if bact else 100000
    microbial_reads = total_reads * (1 - human_pct/100)

    # Import specimen-type exclusion to filter biologically implausible organisms
    try:
        from maps_specimen_scores import should_exclude
    except:
        should_exclude = lambda o,s: False

    results = []
    for org, pct in combined.head(15).items():
        if should_exclude(org, specimen_type):
            continue
        ps        = get_pathogen_score(org)
        ps        = get_specimen_adjusted_score(ps, org, specimen_type)
        is_fungal = get_is_fungal(org)
        try: stype_enc = LE.transform(["unknown"])[0]
        except: stype_enc = 0

        X = pd.DataFrame([[
            pct, human_pct, microbial_reads, ps,
            1 if pct >= THRESHOLD else 0,
            1 if human_pct > 50 else 0,
            pct * ps, np.log1p(microbial_reads),
            1 if microbial_reads < 5000 else 0,
            is_fungal, stype_enc
        ]], columns=FEATURES)

        # Get probability from each model (class 1 = pathogen)
        rf_prob       = RF_MODEL.predict_proba(X)[0][1]
        xgb_prob      = XGB_MODEL.predict_proba(X)[0][1]
        # Calculate weighted ensemble probability
        ensemble_prob = (rf_prob * RF_WEIGHT) + (xgb_prob * XGB_WEIGHT)
        composite     = (ensemble_prob * 0.40) + (min(pct/100,1.0) * 0.35) + (ps * 0.25)
        spec_ps       = get_specimen_adjusted_score(ps, org, specimen_type)
        adj_composite = (ensemble_prob * 0.40) + (min(pct/100,1.0) * 0.35) + (spec_ps * 0.25)
        tier_raw      = "HIGH" if adj_composite >= 0.70 else "MODERATE" if adj_composite >= 0.40 else "LOW"
        # Apply ANSI color codes for terminal display
        # GREEN=HIGH, CYAN=MODERATE, MAGENTA=LOW
        tier          = f"\033[1;32m{tier_raw}\033[0m" if tier_raw=="HIGH" else \
                        f"\033[1;36m{tier_raw}\033[0m" if tier_raw=="MODERATE" else \
                        f"\033[1;35m{tier_raw}\033[0m"
        shap_exp      = get_shap_explanation(X)

        display_org = apply_confidence_naming(org, ensemble_prob)
        results.append({
            "organism":    display_org,
            "abundance":   round(pct,2),
            "rf_score":    round(rf_prob,3),
            "xgb_score":   round(xgb_prob,3),
            "ensemble":    round(ensemble_prob,3),
            "composite":   round(composite,3),
            "tier":        tier,
            "tier_raw":    tier_raw,
            "above_3pct":  pct >= THRESHOLD,
            "shap_drivers":shap_exp,
        })

    results = sorted(results, key=lambda x: -x["composite"])

    print(f"\n    \033[1;34m{'─'*76}\033[0m")
    print(f"    \033[1;34mMAPS ENSEMBLE PREDICTION — {sample}\033[0m")
    print(f"    \033[1;36mRF weight: {RF_WEIGHT*100:.0f}%  |  XGBoost weight: {XGB_WEIGHT*100:.0f}%\033[0m")
    print(f"    \033[1;34m{'─'*76}\033[0m")
    print(f"    {'Rk':<4} {'Organism':<38} {'Abund':>7} {'RF':>7} {'XGB':>7} {'Ens':>7} {'Tier':>18} {'Top SHAP drivers'}")
    print(f"    {'─'*116}")

    for i, r in enumerate(results[:8], 1):
        flag     = f" \033[1;32m✓\033[0m" if r["above_3pct"] else ""
        org_col  = f"\033[1;36m{r['organism']}\033[0m"
        shap_col = f"\033[0;37m{r['shap_drivers']}\033[0m"
        print(f"    {i:<4} {org_col:<49} {r['abundance']:>6.1f}%"
              f" {r['rf_score']:>6.1%} {r['xgb_score']:>6.1%}"
              f" {r['ensemble']:>6.1%} {r['tier']}{flag}"
              f"  {shap_col}")

    top = results[0] if results else None
    if top and top["above_3pct"]:
        org = top["organism"].replace("[BACT] ","").replace("[FUNG] ","")
        print(f"\n    \033[1;34m→ Most likely pathogen: \033[1;32m{org}\033[0m")
        print(f"      RF: \033[1;36m{top['rf_score']:.1%}\033[0m  "
              f"XGBoost: \033[1;36m{top['xgb_score']:.1%}\033[0m  "
              f"Ensemble: \033[1;36m{top['ensemble']:.1%}\033[0m  {top['tier']}")
        if top["rf_score"] >= 0.7 and top["xgb_score"] >= 0.7:
            print(f"      \033[1;32m✓ Both models agree — HIGH confidence result\033[0m")
        elif top["rf_score"] >= 0.7 or top["xgb_score"] >= 0.7:
            print(f"      \033[1;36m~ Models partially agree — interpret with moderate confidence\033[0m")
        else:
            print(f"      \033[1;35m✗ Models disagree — interpret with caution\033[0m")
    else:
        print(f"\n    \033[1;36m→ No organism above {THRESHOLD}% threshold\033[0m")

    # Save ensemble scores to CSV for use in final report and database
    out_path = f"{outdir}/{sample}_ensemble_scores.csv"
    pd.DataFrame(results).to_csv(out_path, index=False)
    print(f"      \033[1;36m✓ Saved: {sample}_ensemble_scores.csv\033[0m")

    return results

if __name__ == "__main__":
    # Parse command line arguments
    # human_pct read from QC JSON if available, otherwise from argv or default 20%
    sample    = sys.argv[1]
    outdir    = sys.argv[2]
    specimen_type = sys.argv[4] if len(sys.argv) > 4 else "unknown"
    import json
    qc_path = f"{outdir}/{sample}_qc.json"
    if len(sys.argv) > 3:
        human_pct = float(sys.argv[3])
    elif os.path.exists(qc_path):
        with open(qc_path) as qf: human_pct = json.load(qf).get("human_pct", 20.0)
    else:
        human_pct = 20.0
    predict_ensemble(sample, outdir, human_pct)
