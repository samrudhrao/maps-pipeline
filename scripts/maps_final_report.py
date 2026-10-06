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

sample  = sys.argv[1]
outdir  = sys.argv[2]
specimen_type = sys.argv[3] if len(sys.argv) > 3 else "unknown"


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

with open(os.path.expanduser("~/maps/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]

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

bact = parse_report(f"{outdir}/bact_results/{sample}_bact_report.txt")
fung = parse_report(f"{outdir}/fungi_results/{sample}_fungi_report.txt")
supp_raw = parse_report(f"{outdir}/bact_results/{sample}_extra_report.txt")

# Check main DB Staphylococcus read counts
KNOWN_STAPH = ['aureus', 'epidermidis', 'lugdunensis', 'saprophyticus',
               'haemolyticus', 'capitis', 'warneri', 'caprae']
main_staph_reads = sum(v[1] for k,v in bact.items()
                       if 'staphylococcus' in k.lower() and
                       any(s in k.lower() for s in KNOWN_STAPH))

for org, (pct, reads) in supp_raw.items():
    if 'cohnii' not in org.lower():
        continue  # Only S. cohnii from supplemental
    if main_staph_reads > 500:
        continue  # Main DB has dominant Staph — skip
    if "staphylococcus" in org.lower():
        if org not in bact or bact[org][1] < reads:
            bact[org] = (pct, reads)

human_report = f"{outdir}/human_depleted/{sample}_human_report.txt"
human_pct = 20.0
total_raw = 100000

# Read from QC json first — most accurate source
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

microbial_reads_est = total_raw * (1 - human_pct/100)

bact_df = pd.Series({k: v[0] for k,v in bact.items()})
fung_df = pd.Series({k: v[0] for k,v in fung.items()})
if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
bact_df.index = ["[BACT] " + x for x in bact_df.index]
fung_df.index = ["[FUNG] " + x for x in fung_df.index]
combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
above    = combined[combined >= 3.0]

# Quality flag
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

# Infection type
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

# AI ranking
model_path = os.path.expanduser("~/maps/models/maps_ai_model.pkl")
ranked_results = []

if os.path.exists(model_path):
    model    = joblib.load(model_path)
    le       = joblib.load(os.path.expanduser("~/maps/models/maps_label_encoder.pkl"))
    features = joblib.load(os.path.expanduser("~/maps/models/maps_features.pkl"))

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

        X = pd.DataFrame([[pct, human_pct, microbial_reads_est, ps,
                            1 if pct >= 3.0 else 0,
                            1 if human_pct > 50 else 0,
                            pct * ps, np.log1p(microbial_reads_est),
                            1 if microbial_reads_est < 5000 else 0,
                            fung_flag, stype_enc]], columns=features)
        prob      = model.predict_proba(X)[0][1]
        composite = (prob * 0.40) + (min(pct/100,1.0) * 0.35) + (ps * 0.25)
        tier      = "HIGH" if prob >= 0.70 else "MODERATE" if prob >= 0.40 else "LOW"

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

# Clinical summary
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

# Print report
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

# Save
df_out   = pd.DataFrame(ranked_results)
csv_path = f"{outdir}/{sample}_final_report.csv"
df_out.to_csv(csv_path, index=False)

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
