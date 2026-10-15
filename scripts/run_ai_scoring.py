#!/usr/bin/env python3
"""
MAPS AI Scoring - run_ai_scoring.py
Step 10 of the MAPS pipeline — Random Forest confidence scoring.
Scores each detected organism above 3% threshold using the trained RF model.

This is a simpler version of maps_predict.py (Step 14) — it runs earlier
in the pipeline to provide preliminary AI scores before the full prediction.
Results saved to OUTDIR/SAMPLE_ai_scores.csv for use in downstream steps.

Primary input: OUTDIR/SAMPLE_abundance.csv (generated in Step 9)
Fallback input: Kraken2 report files if abundance CSV not available

Confidence tiers: HIGH >=70%, MODERATE 40-69%, LOW <40%

Usage: python3 run_ai_scoring.py SAMPLE OUTDIR
Output: OUTDIR/SAMPLE_ai_scores.csv
"""

import joblib
import numpy as np
import pandas as pd
import os

import os as _os
MAPS_HOME = "/maps" if _os.path.exists("/maps/kraken2_db") else _os.path.expanduser("~/maps")
import sys

sample = sys.argv[1]
outdir = sys.argv[2]

# Load contaminant exclusion list
# 137 environmental organisms filtered from results
with open(os.path.expanduser(f"{MAPS_HOME}/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]

# Clinical pathogen prior probability scores
# Must match PATHOGEN_SCORES in maps_final_report.py and maps_predict.py
# Update all three files if adding new organisms
PATHOGEN_SCORES = {
    "pseudomonas aeruginosa": 1.0, "staphylococcus aureus": 1.0,
    "escherichia coli": 0.9,       "klebsiella pneumoniae": 1.0,
    "candida albicans": 1.0,       "aspergillus fumigatus": 0.9,
    "pneumocystis jirovecii": 1.0, "histoplasma": 1.0,
    "coccidioides": 1.0,           "mucor mucedo": 1.0,
    "scedosporium": 1.0,           "lomentospora": 1.0,
    "bacillus licheniformis": 0.7, "microbacterium": 0.6,
    "enterococcus faecium": 0.9,   "citrobacter freundii": 0.8,
    "staphylococcus cohnii": 0.8,  "staphylococcus epidermidis": 0.4,
    "none": 0.0,
}

def get_pathogen_score(org):
    org = str(org).lower()
    for k, v in PATHOGEN_SCORES.items():
        if k in org: return v
    return 0.4

def get_is_fungal(org):
    return 1 if any(f in org.lower() for f in
        ["aspergillus fumigatus","candida","fusarium oxysporum","pneumocystis",
         "histoplasma","mucor mucedo","scedosporium","coccidioides","lomentospora"]) else 0

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

# Load trained Random Forest model
# If model not found pipeline exits gracefully with a warning
model_path = os.path.expanduser(f"{MAPS_HOME}/models/maps_ai_model.pkl")
if not os.path.exists(model_path):
    print("  AI model not found")
    sys.exit(1)

model    = joblib.load(model_path)
le       = joblib.load(os.path.expanduser(f"{MAPS_HOME}/models/maps_label_encoder.pkl"))
features = joblib.load(os.path.expanduser(f"{MAPS_HOME}/models/maps_features.pkl"))

# Use abundance CSV as primary data source
# Contains supplemental DB organisms already merged with correct logic
# Falls back to raw Kraken2 reports if abundance CSV not found
import csv as _csv
abundance_path = f"{outdir}/{sample}_abundance.csv"
bact = {}
fung = {}
if os.path.exists(abundance_path):
    with open(abundance_path) as _f:
        _reader = _csv.DictReader(_f)
        for _row in _reader:
            _org = _row['organism']
            _pct = float(_row['relative_abundance_%'])
            if '[BACT]' in _org:
                bact[_org.replace('[BACT] ','')] = (_pct, int(_pct * 1000))
            elif '[FUNG]' in _org:
                fung[_org.replace('[FUNG] ','')] = (_pct, int(_pct * 1000))
else:
    bact = parse_report(f"{outdir}/bact_results/{sample}_bact_report.txt")
    fung = parse_report(f"{outdir}/fungi_results/{sample}_fungi_report.txt")

bact_df = pd.Series({k: v[0] for k,v in bact.items()})
fung_df = pd.Series({k: v[0] for k,v in fung.items()})
if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
bact_df.index = ["[BACT] " + x for x in bact_df.index]
fung_df.index = ["[FUNG] " + x for x in fung_df.index]
combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
# Only score organisms above 3% threshold
# Below-threshold organisms are tracked in low_abundance_detections column
combined = combined[combined >= 3.0]

# Estimate total reads and human DNA for feature calculation
# human_pct defaults to 20% if QC JSON not available
total_reads = sum(v[1] for v in list(bact.values())[:1]) if bact else 100000
human_pct   = 20.0

print(f"\n  AI Confidence Scores for {sample}:")
print(f"  {'Organism':<45} {'Abund':>7} {'AI Score':>9} {'Tier':>10}")
print("  " + "-"*75)

ai_results = []
for org, pct in combined.head(10).items():
    path_score      = get_pathogen_score(org)
    microbial_reads = total_reads * (1 - human_pct/100)
    above_threshold = 1
    high_human      = 1 if human_pct > 50 else 0
    abund_x_score   = pct * path_score
    log_reads       = np.log1p(microbial_reads)
    low_reads       = 1 if microbial_reads < 5000 else 0
    fung_flag       = get_is_fungal(org)
    try: stype_enc  = le.transform(["unknown"])[0]
    except: stype_enc = 0

    # Build 11-feature vector matching training feature order exactly
    # Order must match maps_features.pkl — do not change without retraining
    vals = [pct, human_pct, microbial_reads, path_score,
            above_threshold, high_human, abund_x_score,
            log_reads, low_reads, fung_flag, stype_enc]

    X    = pd.DataFrame([vals], columns=features)
    prob = model.predict_proba(X)[0][1]
    tier = "HIGH" if prob >= 0.70 else "MODERATE" if prob >= 0.40 else "LOW"
    print(f"  {org:<45} {pct:>6.1f}% {prob:>8.1%} {tier:>10}")
    ai_results.append({"organism": org, "abundance_%": pct,
                        "ai_score": round(prob,3), "tier": tier})

# Save AI scores to CSV for use in Step 13 final report
pd.DataFrame(ai_results).to_csv(f"{outdir}/{sample}_ai_scores.csv", index=False)
print(f"\n  AI scores saved to {outdir}/{sample}_ai_scores.csv")
