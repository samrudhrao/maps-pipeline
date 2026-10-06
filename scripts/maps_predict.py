import pandas as pd
import sys
import os
sys.path.insert(0, os.path.expanduser('~/maps/scripts'))
try:
    from maps_specimen_scores import get_specimen_adjusted_score
except: get_specimen_adjusted_score = lambda s,o,t: s
import numpy as np
import joblib
import os
import sys

model    = joblib.load(os.path.expanduser("~/maps/models/maps_ai_model.pkl"))
le       = joblib.load(os.path.expanduser("~/maps/models/maps_label_encoder.pkl"))
features = joblib.load(os.path.expanduser("~/maps/models/maps_features.pkl"))


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
    "none": 0.0,
}

THRESHOLD = 3.0

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

def predict(sample, outdir, human_pct=20.0):
    # Use abundance CSV — already has supplemental DB merged correctly by pipeline
    bact = {}
    fung = {}
    abundance_path = f"{outdir}/{sample}_abundance.csv"
    if os.path.exists(abundance_path):
        import csv
        combined = {}
        with open(abundance_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                combined[row['organism']] = float(row['relative_abundance_%'])
        combined = pd.Series(combined)
    else:
        # Fallback to Kraken2 reports if no abundance CSV
        bact = parse_report(f"{outdir}/bact_results/{sample}_bact_report.txt")
        fung = parse_report(f"{outdir}/fungi_results/{sample}_fungi_report.txt")
        bact_df = pd.Series({k: v[0] for k,v in bact.items()})
        fung_df = pd.Series({k: v[0] for k,v in fung.items()})
        if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
        if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
        bact_df.index = ["[BACT] " + x for x in bact_df.index]
        fung_df.index = ["[FUNG] " + x for x in fung_df.index]
        combined = pd.concat([bact_df, fung_df])
    combined = combined.sort_values(ascending=False)
    above    = combined[combined >= THRESHOLD]

    total_reads     = sum(v[1] for v in list(bact.values())) if bact else 100000
    microbial_reads = total_reads * (1 - human_pct/100)

    # Import specimen exclusion check
    try:
        from maps_specimen_scores import should_exclude
    except:
        should_exclude = lambda o,s: False

    results = []
    for org, pct in combined.head(15).items():
        if should_exclude(org, specimen_type):
            continue
        ps        = get_pathogen_score(org)
        specimen  = specimen_type if specimen_type != 'unknown' else detect_specimen_type(sample)
        ps        = get_specimen_adjusted_score(ps, org, specimen)
        is_fungal = get_is_fungal(org)
        try: stype_enc = le.transform(["unknown"])[0]
        except: stype_enc = 0

        X = pd.DataFrame([[
            pct, human_pct, microbial_reads, ps,
            1 if pct >= THRESHOLD else 0,
            1 if human_pct > 50 else 0,
            pct * ps,
            np.log1p(microbial_reads),
            1 if microbial_reads < 5000 else 0,
            is_fungal, stype_enc
        ]], columns=features)

        prob      = model.predict_proba(X)[0][1]
        composite = (prob * 0.40) + (min(pct/100,1.0) * 0.35) + (ps * 0.25)
        tier      = "HIGH" if prob >= 0.70 else "MODERATE" if prob >= 0.40 else "LOW"

        display_org = apply_confidence_naming(org, prob)
        results.append({
            "organism":   display_org,
            "abundance":  round(pct, 2),
            "ai_score":   round(prob, 3),
            "composite":  round(composite, 3),
            "tier":       tier,
            "above_3pct": pct >= THRESHOLD
        })

    results = sorted(results, key=lambda x: -x["composite"])

    print(f"\n{'='*65}")
    print(f"MAPS PREDICTION — {sample}")
    print(f"{'='*65}")
    print(f"{'Rank':<5} {'Organism':<42} {'Abund':>7} {'AI':>7} {'Score':>7} {'Tier':>10}")
    print("─"*65)
    for i, r in enumerate(results[:8], 1):
        flag = " ✓" if r["above_3pct"] else ""
        print(f"{i:<5} {r['organism']:<42} {r['abundance']:>6.1f}% "
              f"{r['ai_score']:>6.1%} {r['composite']:>7.3f} {r['tier']:>10}{flag}")

    top = results[0] if results else None
    if top and top["above_3pct"]:
        org = top["organism"].replace("[BACT] ","").replace("[FUNG] ","")
        print(f"\n→ Most likely pathogen: {org}")
        print(f"  Abundance: {top['abundance']}% | AI: {top['ai_score']:.1%} | {top['tier']}")
    else:
        print(f"\n→ No organism above {THRESHOLD}% threshold")

    return results

def detect_specimen_type(sample):
    sample_lower = sample.lower()
    if "syn" in sample_lower: return "synovial"
    if "pla" in sample_lower: return "plasma"
    if "bal" in sample_lower: return "bal"
    if "pet" in sample_lower: return "peritoneal"
    if "plu" in sample_lower: return "pleural"
    if "abs" in sample_lower: return "abscess"
    return "unknown"

if __name__ == "__main__":
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
    predict(sample, outdir, human_pct)
