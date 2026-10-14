#!/usr/bin/env python3
"""
MAPS Statistical Analysis - maps_stats.py
Step 16 of the MAPS pipeline — confidence intervals and anomaly detection.

Two main functions:
  1. Wilson score confidence intervals for each detected organism
     - Better than normal approximation for proportions near 0 or 100%
     - 95% confidence level by default
     - Precision interpreted as Very High/High/Moderate/Low

  2. Isolation Forest anomaly detection
     - Trained on 26 confirmed clinical samples from BSW
     - Flags samples with unusual microbial profiles
     - Score < -0.15 = HIGH anomaly, < -0.05 = MODERATE anomaly
     - Anomaly does NOT invalidate results — flags for manual review

To retrain anomaly model: update training data in train_anomaly_model()
and delete /maps/models/maps_anomaly_model.pkl to force retraining.

Usage: python3 maps_stats.py SAMPLE OUTDIR [SPECIMEN_TYPE]
Output: OUTDIR/SAMPLE_stats.json and OUTDIR/SAMPLE_confidence_intervals.csv
"""

import numpy as np
import pandas as pd
import joblib
import os
import sys
import json
from scipy.stats import norm
from sklearn.ensemble import IsolationForest

# ── CONFIDENCE INTERVALS ──────────────────────────────────────────
# Wilson score interval is preferred over normal approximation
# because it works correctly for proportions near 0% or 100%
def wilson_confidence_interval(reads_organism, total_reads, confidence=0.95):
    """
    Calculate Wilson score confidence interval for a proportion.
    reads_organism: reads assigned to this organism
    total_reads: total microbial reads in sample
    confidence: confidence level (default 0.95 = 95% CI)
    Returns (lower, upper) as percentages
    """
    """Wilson score interval — best method for proportions"""
    if total_reads == 0:
        return 0.0, 0.0
    p = reads_organism / total_reads
    z = norm.ppf((1 + confidence) / 2)
    denom  = 1 + z**2/total_reads
    center = (p + z**2/(2*total_reads)) / denom
    margin = z * np.sqrt(p*(1-p)/total_reads + z**2/(4*total_reads**2)) / denom
    lower  = max(0, (center - margin) * 100)
    upper  = min(100, (center + margin) * 100)
    return round(lower, 2), round(upper, 2)

def interpret_ci(lower, upper, abundance):
    """
    Interpret the width of a confidence interval as a precision label.
    Narrow CI = high precision (many reads), wide CI = low precision (few reads).
    """
    """Interpret width of confidence interval"""
    width = upper - lower
    if width < 2:
        return "Very high precision"
    elif width < 5:
        return "High precision"
    elif width < 15:
        return "Moderate precision"
    else:
        return "Low precision — interpret with caution"

# ── ANOMALY DETECTION ─────────────────────────────────────────────
# Isolation Forest detects samples that differ significantly from
# the 26 training samples. Anomaly = unusual microbial profile,
# NOT necessarily a failed sample or incorrect result.
def train_anomaly_model():
    """
    Train Isolation Forest anomaly detector on known BSW clinical samples.
    Training data: 26 confirmed samples with known abundance, human_pct, total_reads.
    Model saved to /maps/models/maps_anomaly_model.pkl
    Called automatically if model file not found.
    To retrain: delete the pkl file and run any sample through the pipeline.
    """
    """Train Isolation Forest on known samples"""
    training = [
        (0.01,  8.87,  633193),
        (98.91, 36.93, 568662),
        (78.05, 9.16,  608055),
        (95.99, 14.80, 221738),
        (9.83,  21.55, 579728),
        (9.36,  5.44,  531291),
        (62.07, 7.95,  508839),
        (8.29,  6.60,  231896),
        (30.13, 8.87,  633193),
        (0.17,  36.93, 568662),
        (16.24, 9.16,  608055),
        (0.04,  85.20, 221738),
        (18.67, 14.80, 579728),
        (9.59,  21.55, 531291),
        (0.0,   80.00, 508839),
        (23.89, 5.44,  231896),
        (0.09,  7.95,  633193),
        (0.49,  81.00, 568662),
        (0.0,   84.00, 608055),
        (8.10,  92.05, 221738),
        (0.85,  78.45, 579728),
        (0.0,   93.40, 531291),
        (52.14, 23.02, 633193),
        (59.46, 27.75, 568662),
        (57.30, 11.14, 608055),
        (22.12, 19.35, 221738),
    ]

    df = pd.DataFrame(training, columns=["abundance","human_pct","total_reads"])
    df["microbial_reads"] = df["total_reads"] * (1 - df["human_pct"]/100)
    df["log_reads"]       = np.log1p(df["microbial_reads"])
    df["log_abundance"]   = np.log1p(df["abundance"])

    X = df[["abundance","human_pct","log_reads","log_abundance"]].values

    iso = IsolationForest(n_estimators=200, contamination=0.15, random_state=42)
    iso.fit(X)

    model_path = os.path.expanduser("/maps/models/maps_anomaly_model.pkl")
    joblib.dump(iso, model_path)
    return iso

def detect_anomaly(abundance, human_pct, total_reads):
    """
    Check if sample profile is anomalous compared to training data.
    abundance: top organism abundance percentage
    human_pct: human DNA percentage
    total_reads: total read count
    Returns dict with is_anomaly, score, severity, reason
    """
    """Check if sample looks anomalous vs training data"""
    model_path = os.path.expanduser("/maps/models/maps_anomaly_model.pkl")

    if not os.path.exists(model_path):
        iso = train_anomaly_model()
    else:
        iso = joblib.load(model_path)

    microbial_reads = total_reads * (1 - human_pct/100)
    log_reads       = np.log1p(microbial_reads)
    log_abundance   = np.log1p(abundance)

    X       = np.array([[abundance, human_pct, log_reads, log_abundance]])
    score   = iso.decision_function(X)[0]
    is_anom = iso.predict(X)[0] == -1

    if score < -0.15:
        severity = "HIGH"
        reason   = "Very unusual sample profile"
    elif score < -0.05:
        severity = "MODERATE"
        reason   = "Somewhat unusual sample profile"
    else:
        severity = "LOW"
        reason   = "Sample profile within normal range"

    return {
        "is_anomaly": is_anom,
        "score":      round(float(score), 4),
        "severity":   severity if is_anom else "NONE",
        "reason":     reason if is_anom else "Normal",
    }

# ── MAIN FUNCTION ─────────────────────────────────────────────────
# Runs confidence intervals and anomaly detection for a sample
# Called by full_pipeline.sh at Step 16
def run_stats(sample, outdir):
    with open(os.path.expanduser("/maps/filter_list.txt")) as f:
        EXCLUDE = [line.strip().lower() for line in f if line.strip()]

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

    # Use abundance CSV as primary data source
    # This ensures consistency with all other pipeline steps
    # Falls back to Kraken2 reports if abundance CSV not available
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

    # Apply specimen-type exclusions to remove biologically implausible organisms
    # e.g. gut bacteria should not appear in synovial fluid results
    try:
        import sys as _sys
        _sys.path.insert(0, os.path.expanduser("/maps/scripts"))
        from maps_specimen_scores import should_exclude
        bact = {k:v for k,v in bact.items() if not should_exclude(k, specimen_type)}
        fung = {k:v for k,v in fung.items() if not should_exclude(k, specimen_type)}
    except:
        pass

    bact_df = pd.Series({k: v[0] for k,v in bact.items()})
    fung_df = pd.Series({k: v[0] for k,v in fung.items()})
    if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
    if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
    bact_df.index = ["[BACT] " + x for x in bact_df.index]
    fung_df.index = ["[FUNG] " + x for x in fung_df.index]
    combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
    above    = combined[combined >= 3.0]

    total_reads = sum(v[1] for v in list(bact.values())) if bact else 100000

    # Get human DNA percentage for anomaly detection
    # Read from Kraken2 human depletion report
    human_pct   = 20.0
    human_report = f"{outdir}/human_depleted/{sample}_human_report.txt"
    if os.path.exists(human_report):
        with open(human_report) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6 and "Homo sapiens" in parts[5]:
                    human_pct = float(parts[0])

    print(f"\n{'='*70}")
    print(f"MAPS STATISTICAL ANALYSIS — {sample}")
    print(f"{'='*70}")

    # Calculate Wilson score confidence intervals for each organism above 3% threshold
    # Wider CI = fewer reads supporting the detection = less precise estimate
    print(f"\n[CONFIDENCE INTERVALS — 95%]")
    print(f"  {'Organism':<42} {'Abundance':>10} {'95% CI':>20} {'Precision':>25}")
    print("  " + "─"*100)

    ci_results = []
    for org, pct in above.items():
        reads_org   = int(pct/100 * total_reads)
        lower, upper = wilson_confidence_interval(reads_org, total_reads)
        precision   = interpret_ci(lower, upper, pct)
        print(f"  {org:<42} {pct:>9.1f}% "
              f"  ({lower:.1f}% - {upper:.1f}%){precision:>25}")
        ci_results.append({
            "organism":   org,
            "abundance":  round(pct, 2),
            "ci_lower":   lower,
            "ci_upper":   upper,
            "ci_width":   round(upper-lower, 2),
            "precision":  precision,
        })

    if not ci_results:
        print("  No organisms above 3% threshold")

    # Run Isolation Forest anomaly detection on top organism abundance
    # Anomaly flag does not invalidate results — flags for physician review
    print(f"\n[ANOMALY DETECTION]")
    top_abundance = above.iloc[0] if len(above) > 0 else 0.0
    anomaly = detect_anomaly(top_abundance, human_pct, total_reads)

    status = "⚠ ANOMALY DETECTED" if anomaly["is_anomaly"] else "✓ NORMAL"
    print(f"  Status:   {status}")
    print(f"  Score:    {anomaly['score']:.4f} (more negative = more anomalous)")
    print(f"  Severity: {anomaly['severity']}")
    print(f"  Reason:   {anomaly['reason']}")

    if anomaly["is_anomaly"]:
        print(f"\n  ⚠ This sample has an unusual profile compared to your")
        print(f"    26 training samples. Consider:")
        print(f"    - Checking for sequencing run issues")
        print(f"    - Verifying sample preparation")
        print(f"    - Reviewing raw data quality in FastQC report")

    # Save statistical results to JSON and CSV for downstream use
    # JSON included in pipeline output, CSV for data analysis
    stats_output = {
        "sample":            sample,
        "confidence_intervals": ci_results,
        "anomaly_detection": anomaly,
        "total_reads":       total_reads,
        "human_pct":         round(human_pct, 1),
    }
    json_path = f"{outdir}/{sample}_stats.json"
    with open(json_path, "w") as f:
        json.dump(stats_output, f, indent=2, default=lambda o: int(o) if isinstance(o, bool) else str(o))

    csv_path = f"{outdir}/{sample}_confidence_intervals.csv"
    if ci_results:
        pd.DataFrame(ci_results).to_csv(csv_path, index=False)
        print(f"\nSaved: {sample}_confidence_intervals.csv")
    print(f"Saved: {sample}_stats.json")
    print(f"{'='*70}")

    return stats_output

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 maps_stats.py <sample> <outdir>")
        sys.exit(1)
    sample = sys.argv[1]
    outdir = sys.argv[2]
    specimen_type = sys.argv[3] if len(sys.argv) > 3 else "unknown" 
    run_stats(sample, outdir)
