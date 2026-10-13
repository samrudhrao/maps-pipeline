#!/usr/bin/env python3
"""
MAPS Final Clinical Report - maps_final_report.py
Generates the ranked pathogen report used for clinical interpretation.
This is the primary output physicians use to make treatment decisions.

Key functions:
  - Merges bacterial, fungal, and supplemental DB results
  - Applies 3% abundance threshold for clinical reporting
  - Applies genus-level naming for organisms below 70% AI confidence
  - Applies specimen-type specific scoring adjustments

Supplemental DB logic (S. cohnii and B. licheniformis):
  - S. cohnii: only added if main DB has <500 reads of known pathogenic Staph
  - B. licheniformis: always added if >=1% abundance and not in main DB

Genus-level naming:
  - Organisms with ensemble AI score <70% reported as "Genus sp."
  - Full species name reported when score >=70%

Usage: python3 maps_final_report.py SAMPLE OUTDIR HUMAN_PCT SPECIMEN_TYPE
Output: OUTDIR/SAMPLE_final_report.csv and OUTDIR/SAMPLE_summary.json
"""

import pandas as pd
import sys
import os
sys.path.insert(0, os.path.expanduser('~/maps/scripts'))
try:
    from maps_specimen_scores import get_specimen_adjusted_score, detect_specimen_type
except: 
    get_specimen_adjusted_score = lambda s,o,t: s
    detect_specimen_type = lambda s: "unknown"

import numpy as np
import os
import joblib
import sys
import json
from datetime import datetime

# Command line arguments
sample  = sys.argv[1]       # Sample name
outdir  = sys.argv[2]       # Output directory
outdir  = sys.argv[2]
specimen_type = sys.argv[3] if len(sys.argv) > 3 else "unknown"  # Clinical specimen type


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

# Load contaminant filter list (137 environmental/reagent organisms to exclude)
# Located at ~/maps/filter_list.txt — add new contaminants here
with open(os.path.expanduser("~/maps/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]

# Clinical pathogen scores — prior probability that organism is a true pathogen
# Scale 0.0-1.0: 1.0 = always pathogenic, 0.4 = common contaminant
# Used to weight AI confidence scores for clinical reporting
# To add new organisms: append to this dictionary with appropriate score
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
    "streptococcus pneumoniae": 1.0,"none": 0.0,
}

def get_pathogen_score(org):
    """
    Look up clinical pathogen score for an organism.
    Returns 0.4 (default) for organisms not in PATHOGEN_SCORES.
    Higher score = higher prior probability of being a true pathogen.
    """
    org = str(org).lower()
    for k, v in PATHOGEN_SCORES.items():
        if k in org: return v
    return 0.4

def get_is_fungal(org):
    """
    Check if organism is a clinically significant fungus.
    Used as a binary feature (0/1) in the AI model.
    Add new fungal pathogens to the list if needed.
    """
    return 1 if any(f in org.lower() for f in
        ["aspergillus fumigatus","candida","fusarium oxysporum","pneumocystis",
         "histoplasma","mucor mucedo","scedosporium","coccidioides","lomentospora"]) else 0

def parse_report(path, min_reads=50):
    """
    Parse Kraken2 report file and return species-level detections.
    Filters to species rank (S) with minimum 50 reads.
    Excludes organisms in the EXCLUDE contaminant list.
    Returns dict: {organism_name: (abundance_pct, read_count)}
    """
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

# Load results from all three Kraken2 databases
# bact = standard bacterial DB, fung = fungal DB, supp_raw = supplemental DB
bact = parse_report(f"{outdir}/bact_results/{sample}_bact_report.txt")
fung = parse_report(f"{outdir}/fungi_results/{sample}_fungi_report.txt")
supp_raw = parse_report(f"{outdir}/bact_results/{sample}_extra_report.txt")

# Check how many reads from known pathogenic Staphylococcus species are in main DB
# This determines whether S. cohnii from supplemental DB should be included
# If main DB already has dominant Staph (>500 reads), S. cohnii is likely noise
KNOWN_STAPH = ['aureus', 'epidermidis', 'lugdunensis', 'saprophyticus',
               'haemolyticus', 'capitis', 'warneri', 'caprae',
               'hominis', 'pasteuri', 'pettenkoferi', 'simulans', 'xylosus']
main_staph_reads = sum(v[1] for k,v in bact.items()
                       if 'staphylococcus' in k.lower() and
                       any(s in k.lower() for s in KNOWN_STAPH))

for org, (pct, reads) in supp_raw.items():
    # Only process S. cohnii and B. licheniformis from supplemental
    if 'cohnii' not in org.lower() and 'licheniformis' not in org.lower():
        continue
    # Skip if below meaningful signal threshold
    if pct < 1.0:
        continue
    # S. cohnii — skip if main DB already has dominant Staphylococcus
    if 'cohnii' in org.lower() and main_staph_reads > 500:
        continue
    # Add to bact if not present or supplemental has more signal
    if org not in bact or bact[org][1] < reads:
        bact[org] = (pct, reads)

# Get human DNA percentage and total read count for quality assessment
# Primary source: QC JSON generated in Step 3b
# Fallback: parse Kraken2 human depletion report directly
human_report = f"{outdir}/human_depleted/{sample}_human_report.txt"
human_pct = 20.0   # Default value if QC data not available
total_raw = 100000  # Default value if QC data not available

# Read from QC JSON first — most accurate source
# Falls back to parsing Kraken2 report if QC JSON not found
qc_path = f"{outdir}/{sample}_qc.json"
if os.path.exists(qc_path):
    with open(qc_path) as f:
        qc = json.load(f)
    human_pct  = qc.get("human_pct", 20.0)
    total_raw  = qc.get("total_reads", 100000)
elif os.path.exists(human_report):
    with open(human_report) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6 and "Homo sapiens" in parts[5]:
                human_pct = float(parts[0])
            if len(parts) >= 4 and parts[3].strip() == "R":
                try: total_raw = int(parts[1])
                except: pass

# Estimate microbial read count from total reads minus human fraction
microbial_reads_est = total_raw * (1 - human_pct/100)

# Normalize bacterial and fungal abundances to relative percentages
# Then add prefixes and combine into single ranked series
bact_df = pd.Series({k: v[0] for k,v in bact.items()})
fung_df = pd.Series({k: v[0] for k,v in fung.items()})
if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
bact_df.index = ["[BACT] " + x for x in bact_df.index]
fung_df.index = ["[FUNG] " + x for x in fung_df.index]
combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
# Apply 3% minimum abundance threshold for clinical reporting
above    = combined[combined >= 3.0]

# Assign quality flag based on microbial read count and human DNA percentage
# CRITICAL: results unreliable, LOW: interpret with caution, PASS: acceptable
if microbial_reads_est < 1000:
    quality_flag = "CRITICAL"
    quality_msg  = f"Fewer than 1,000 microbial reads ({microbial_reads_est:,.0f}). Results unreliable."
elif microbial_reads_est < 5000:
    quality_flag = "LOW"
    quality_msg  = f"Fewer than 5,000 microbial reads ({microbial_reads_est:,.0f}). Interpret with caution."
elif human_pct > 80:
    quality_flag = "LOW"
    quality_msg  = f"{human_pct:.1f}% human DNA. Very high host contamination."
elif human_pct > 50:
    quality_flag = "MODERATE"
    quality_msg  = f"{human_pct:.1f}% human DNA may limit sensitivity."
else:
    quality_flag = "GOOD"
    quality_msg  = f"Sufficient microbial reads ({microbial_reads_est:,.0f}). Human DNA: {human_pct:.1f}%."

# Classify infection type based on number of organisms above threshold
# MONOMICROBIAL: single pathogen, POLYMICROBIAL: multiple pathogens
# DOMINATED: one organism clearly dominant, COMPLEX: truly mixed infection
above_list = [(org, pct) for org, pct in above.items()]
if len(above_list) == 0:
    infection_type = "NO DETECTION"
elif len(above_list) == 1:
    infection_type = "MONOMICROBIAL"
elif len(above_list) == 2:
    infection_type = "POLYMICROBIAL"
else:
    max_pct    = above_list[0][1]
    second_pct = above_list[1][1]
    if max_pct > 60 and second_pct < 15:
        infection_type = "DOMINATED POLYMICROBIAL"
    else:
        infection_type = "COMPLEX POLYMICROBIAL"

# Load AI model and rank organisms by composite score
# Composite = (AI probability * 0.40) + (abundance * 0.35) + (pathogen score * 0.25)
# Top 15 organisms from combined bacterial+fungal results are scored
model_path = os.path.expanduser("~/maps/models/maps_ai_model.pkl")
ranked_results = []

if os.path.exists(model_path):
    model    = joblib.load(model_path)
    le       = joblib.load(os.path.expanduser("~/maps/models/maps_label_encoder.pkl"))
    features = joblib.load(os.path.expanduser("~/maps/models/maps_features.pkl"))

    # Import specimen-type exclusion function
    # Some organisms are biologically implausible in certain specimen types
    # e.g. gut bacteria in synovial fluid are likely contaminants
    try:
        from maps_specimen_scores import should_exclude
    except:
        should_exclude = lambda o,s: False

    for org, pct in combined.head(15).items():
        if should_exclude(org, specimen_type):
            continue
        ps        = get_pathogen_score(org)
        ps        = get_specimen_adjusted_score(ps, org, specimen_type)
        fung_flag = get_is_fungal(org)
        try: stype_enc = le.transform([specimen_type])[0]
        except: stype_enc = 0

        # Build feature vector for AI model prediction
        # Must match exactly the 11 features used during model training
        X = pd.DataFrame([[pct, human_pct, microbial_reads_est, ps,
                            1 if pct >= 3.0 else 0,
                            1 if human_pct > 50 else 0,
                            pct * ps, np.log1p(microbial_reads_est),
                            1 if microbial_reads_est < 5000 else 0,
                            fung_flag, stype_enc]], columns=features)
        prob      = model.predict_proba(X)[0][1]
        composite = (prob * 0.40) + (min(pct/100,1.0) * 0.35) + (ps * 0.25)
        tier      = "HIGH" if prob >= 0.70 else "MODERATE" if prob >= 0.40 else "LOW"

        # Apply genus-level naming if AI confidence below 70%
        display_org = apply_confidence_naming(org, prob)
        ranked_results.append({
            "rank": 0, "organism": display_org,
            "abundance_%": round(pct, 2),
            "pathogen_score": round(ps, 2),
            "ai_score": round(prob, 3),
            "composite_rank": round(composite, 3),
            "tier": tier,
            "above_3pct": "YES" if pct >= 3.0 else "NO",
        })

    ranked_results = sorted(ranked_results, key=lambda x: -x["composite_rank"])
    for i, r in enumerate(ranked_results): r["rank"] = i + 1

# Generate plain-language clinical summary based on top organism
# Used in the final report and stored in summary JSON for database
if ranked_results and ranked_results[0]["above_3pct"] == "YES":
    top      = ranked_results[0]
    org_name = top["organism"].replace("[BACT] ","").replace("[FUNG] ","")
    pct      = top["abundance_%"]
    ai_pct   = top["ai_score"] * 100
    if pct > 50 and top["tier"] == "HIGH":
        clinical_summary = (f"Dominant detection of {org_name} at {pct:.1f}% microbial DNA abundance "
                           f"with {ai_pct:.0f}% AI confidence strongly suggests this is the primary pathogen.")
    elif pct > 20 and top["tier"] == "HIGH":
        clinical_summary = (f"Significant detection of {org_name} at {pct:.1f}% microbial DNA abundance "
                           f"with {ai_pct:.0f}% AI confidence is consistent with active infection.")
    else:
        clinical_summary = (f"{org_name} detected at {pct:.1f}% microbial DNA abundance "
                           f"with {ai_pct:.0f}% AI confidence. Correlate with clinical presentation.")
else:
    clinical_summary = "No pathogen detected above 3% threshold."

# Print formatted clinical report to stdout
# This output is captured by the pipeline and shown in terminal
print(f"\n{'='*70}")
print(f"MAPS FINAL REPORT — {sample}")
print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
print(f"{'='*70}")
print(f"\n[QUALITY] {quality_flag}")
print(f"  {quality_msg}")
print(f"\n[INFECTION TYPE] {infection_type}")
print(f"\n[CLINICAL SUMMARY]")
print(f"  {clinical_summary}")

if ranked_results:
    print(f"\n[RANKED PATHOGENS]")
    print(f"  {'Rk':<4} {'Organism':<42} {'Abund':>7} {'AI%':>7} {'Score':>7} {'Tier':>10} {'>3%':>5}")
    print("  " + "─"*82)
    for r in ranked_results[:8]:
        print(f"  {r['rank']:<4} {r['organism']:<42} {r['abundance_%']:>6.1f}% "
              f"{r['ai_score']:>6.1%} {r['composite_rank']:>7.3f} "
              f"{r['tier']:>10} {r['above_3pct']:>5}")

print(f"\n{'='*70}")

# Save ranked results to CSV and summary JSON
# CSV used by downstream steps (maps_predict.py, maps_ensemble.py, maps_stats.py)
# JSON used by maps_database.py to store results in SQLite
df_out   = pd.DataFrame(ranked_results)
csv_path = f"{outdir}/{sample}_final_report.csv"
df_out.to_csv(csv_path, index=False)

# Build summary dictionary for JSON output
# This is the primary data structure passed to maps_database.py
top_result = ranked_results[0] if ranked_results else {}
summary = {
    "sample":              sample,
    "quality":             quality_flag,
    "human_pct":           round(human_pct,1),
    "microbial_reads":     int(microbial_reads_est),
    "infection_type":      infection_type,
    "organisms_above_3pct":len(above_list),
    "top_organism":        top_result.get("organism","None"),
    "top_ai_score":        top_result.get("ai_score",0),
    "top_composite":       top_result.get("composite_rank",0),
    "clinical_summary":    clinical_summary,
    "timestamp":           datetime.now().strftime('%Y-%m-%d %H:%M'),
}
json_path = f"{outdir}/{sample}_summary.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)

print(f"Saved: {sample}_final_report.csv")
print(f"Saved: {sample}_summary.json")
