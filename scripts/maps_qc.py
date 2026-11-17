#!/usr/bin/env python3
"""
MAPS QC Module - maps_qc.py
Validates sample quality after human depletion before running the full pipeline.
Reads Kraken2 human depletion report and fastp JSON to calculate QC metrics.

Hard fail criteria (pipeline stops if any are met):
  - Total reads < 50,000
  - Human DNA > 90%
  - Microbial reads < 1,000
  - Illumina: Q30 rate < 50%  |  Ion Torrent: Q20 rate < 70%
Soft warnings (logged, do not stop the pipeline):
  - Duplication > 80%
  - Insert size < 100bp

Usage:  python3 maps_qc.py SAMPLE OUTDIR [illumina|iontorrent]
Output: OUTDIR/SAMPLE_qc.json — includes "failures" and "warnings" lists
        naming each parameter that failed, its value and the threshold.
"""

import sys
import os
import json


def qc_check(sample, outdir, platform='illumina'):
    human_report    = f"{outdir}/human_depleted/{sample}_human_report.txt"
    fastp_json      = f"{outdir}/fastp/{sample}_fastp.json"
    human_pct       = 0.0
    total_reads     = 0
    unclassified    = 0
    classified      = 0
    duplication_pct = 0.0
    q30_rate        = 0.0
    q20_rate        = 0.0
    insert_size     = 0

    if os.path.exists(human_report):
        with open(human_report) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6 and "Homo sapiens" in parts[5]:
                    human_pct = float(parts[0])
                if len(parts) >= 4 and parts[3].strip() == "U":
                    try: unclassified = int(parts[1])
                    except: pass
                if len(parts) >= 4 and parts[3].strip() == "R":
                    try: classified = int(parts[1])
                    except: pass
        total_reads = unclassified + classified

    if os.path.exists(fastp_json):
        try:
            with open(fastp_json) as f:
                fastp = json.load(f)
            summary = fastp.get("summary", {})
            after   = summary.get("after_filtering", {})
            q30_rate        = after.get("q30_rate", 0.0) * 100
            q20_rate        = after.get("q20_rate", 0.0) * 100
            duplication_pct = fastp.get("duplication", {}).get("rate", 0.0) * 100
            insert_size     = fastp.get("insert_size", {}).get("peak", 0)
        except:
            pass

    microbial_reads = total_reads * (1 - human_pct / 100)

    # Quality metric depends on platform — per APHL guidelines
    # Illumina: hard fail below 50% Q30 | Ion Torrent: hard fail below 70% Q20
    if platform == "iontorrent":
        quality_name, quality_value, quality_min = "Q20 rate", q20_rate, 70.0
    else:
        quality_name, quality_value, quality_min = "Q30 rate", q30_rate, 50.0
    # A rate of 0 means fastp metrics were unavailable — not treated as a failure
    quality_ok = not (0 < quality_value < quality_min)

    checks = {
        "total_reads_ok":     total_reads >= 50000,
        "human_pct_ok":       human_pct <= 90,
        "microbial_reads_ok": microbial_reads >= 1000,
        "quality_ok":         quality_ok,
    }

    failures = []
    if not checks["total_reads_ok"]:
        failures.append({"parameter": "Total reads", "value": f"{total_reads:,}",
                         "threshold": ">= 50,000",
                         "reason": "Insufficient sequencing depth"})
    if not checks["human_pct_ok"]:
        failures.append({"parameter": "Human DNA", "value": f"{human_pct:.1f}%",
                         "threshold": "<= 90%",
                         "reason": "Too much host (human) DNA"})
    if not checks["microbial_reads_ok"]:
        failures.append({"parameter": "Microbial reads", "value": f"{microbial_reads:,.0f}",
                         "threshold": ">= 1,000",
                         "reason": "Too few non-human reads to classify"})
    if not checks["quality_ok"]:
        failures.append({"parameter": quality_name, "value": f"{quality_value:.1f}%",
                         "threshold": f">= {quality_min:.0f}%",
                         "reason": "Poor base-call quality"})

    warnings = []
    if duplication_pct > 80:
        warnings.append({"parameter": "Duplication", "value": f"{duplication_pct:.1f}%",
                         "threshold": "<= 80%",
                         "reason": "High duplication — possible PCR artifact or low complexity library"})
    if 0 < insert_size < 100:
        warnings.append({"parameter": "Insert size", "value": f"{insert_size}bp",
                         "threshold": ">= 100bp",
                         "reason": "Short insert size — possible DNA degradation"})

    passed = not failures
    status = "PASS" if passed else "FAIL"

    print(f"\n{'='*55}")
    print(f"QC REPORT — {sample} ({platform})")
    print(f"{'='*55}")
    print(f"  {'✓' if checks['total_reads_ok'] else '✗'} Total reads: {total_reads:,} (need 50,000+)")
    print(f"  {'✓' if checks['human_pct_ok'] else '✗'} Human DNA: {human_pct:.1f}% (max 90%)")
    print(f"  {'✓' if checks['microbial_reads_ok'] else '✗'} Microbial reads: {microbial_reads:,.0f} (need 1,000+)")
    if quality_value > 0:
        print(f"  {'✓' if quality_ok else '✗'} {quality_name}: {quality_value:.1f}% (need {quality_min:.0f}%+)")
    else:
        print(f"  ? {quality_name}: not available (fastp metrics missing)")
    if platform != "iontorrent" and quality_ok and 0 < q30_rate < 70:
        print(f"  ⚠ Q30 rate low ({q30_rate:.1f}%) — results should be interpreted with caution")
    if duplication_pct > 0:
        print(f"  {'✓' if duplication_pct <= 80 else '⚠'} Duplication: {duplication_pct:.1f}% (warning above 80%)")
    if insert_size > 0:
        print(f"  {'✓' if insert_size >= 100 else '⚠'} Insert size: {insert_size}bp (warning below 100bp)")

    print(f"\n  Final QC status: {status}")
    if failures:
        print("  Failed parameters: " + "; ".join(
            f"{f['parameter']} {f['value']} (required {f['threshold']})" for f in failures))
        print(f"  ⚠ Sample {sample} failed QC — stopping pipeline")

    qc = {
        "sample":           sample,
        "platform":         platform,
        "status":           status,
        "total_reads":      total_reads,
        "human_pct":        round(human_pct, 1),
        "microbial_reads":  int(microbial_reads),
        "q30_rate":         round(q30_rate, 1),
        "q20_rate":         round(q20_rate, 1),
        "quality_metric":   quality_name,
        "duplication_pct":  round(duplication_pct, 1),
        "insert_size":      insert_size,
        "checks":           checks,
        "failures":         failures,
        "warnings":         warnings,
    }
    with open(f"{outdir}/{sample}_qc.json", "w") as f:
        json.dump(qc, f, indent=2)

    return passed

if __name__ == "__main__":
    sample   = sys.argv[1]
    outdir   = sys.argv[2]
    platform = sys.argv[3] if len(sys.argv) > 3 else "illumina"
    qc_check(sample, outdir, platform)
