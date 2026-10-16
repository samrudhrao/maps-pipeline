#!/usr/bin/env python3
"""
MAPS QC Module - maps_qc.py
Validates sample quality after human depletion before running the full pipeline.
Reads Kraken2 human depletion report and fastp JSON to calculate QC metrics.

Hard fail criteria (pipeline stops if any are met):
  - Total reads < 50,000
  - Human DNA > 90%
  - Microbial reads < 1,000
  - Q30 rate < 50%

Usage: python3 maps_qc.py SAMPLE OUTDIR
Output: OUTDIR/SAMPLE_qc.json
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
            duplication_pct = fastp.get("duplication", {}).get("rate", 0.0) * 100
            insert_size     = fastp.get("insert_size", {}).get("peak", 0)
        except:
            pass

    microbial_reads = total_reads * (1 - human_pct / 100)

    # Core pass/fail checks — unchanged
    checks = {
        "total_reads_ok":     (total_reads >= 50000,
                               f"Total reads: {total_reads:,} (need 50,000+)"),
        "human_pct_ok":       (human_pct <= 90,
                               f"Human DNA: {human_pct:.1f}% (max 90%)"),
        "microbial_reads_ok": (microbial_reads >= 1000,
                               f"Microbial reads: {microbial_reads:,.0f} (need 1,000+)"),
    }

    passed = all(v[0] for v in checks.values())
    status = "PASS" if passed else "FAIL"

    print(f"\n{'='*55}")
    print(f"QC REPORT — {sample}: checking...")
    print(f"{'='*55}")
    for check, (ok, msg) in checks.items():
        icon = "✓" if ok else "✗"
        print(f"  {icon} {msg}")

    # Informational metrics — logged but not pass/fail
    if q30_rate > 0:
        q30_flag = "✓" if q30_rate >= 70 else "⚠"
        print(f"  {q30_flag} Q30 rate: {q30_rate:.1f}% (informational)")
    if duplication_pct > 0:
        dup_flag = "✓" if duplication_pct <= 80 else "⚠"
        print(f"  {dup_flag} Duplication: {duplication_pct:.1f}% (informational)")
        if duplication_pct > 80:
            print(f"  ⚠ High duplication rate — possible PCR artifact or low complexity library")
    if insert_size > 0:
        ins_flag = "✓" if insert_size >= 100 else "⚠"
        print(f"  {ins_flag} Insert size: {insert_size}bp (informational)")
        if insert_size < 100:
            print(f"  ⚠ Short insert size — possible DNA degradation")

    # Q30 threshold depends on platform
    # Illumina: hard fail below 50% (platform standard)
    # Ion Torrent: hard fail below 5% (platform naturally lower Q30)
    q30_min = 5.0 if platform == "iontorrent" else 50.0
    if q30_rate > 0 and q30_rate < q30_min:
        passed = False
        status = "FAIL"
        print(f"\n  ✗ Q30 rate critically low ({q30_rate:.1f}%) — sequencing quality insufficient for reliable results")
    elif q30_rate > 0 and q30_rate < 70:
        if platform != "iontorrent":
            print(f"\n  ⚠ Q30 rate low ({q30_rate:.1f}%) — results should be interpreted with caution")

    # Reprint status header with correct pass/fail
    print(f"\n  Final QC status: {status}")

    if not passed:
        print(f"  ⚠ Sample {sample} failed QC — stopping pipeline")

    qc = {
        "sample":           sample,
        "status":           status,
        "total_reads":      total_reads,
        "human_pct":        round(human_pct, 1),
        "microbial_reads":  int(microbial_reads),
        "q30_rate":         round(q30_rate, 1),
        "duplication_pct":  round(duplication_pct, 1),
        "insert_size":      insert_size,
        "checks":           {k: v[0] for k, v in checks.items()},
    }
    with open(f"{outdir}/{sample}_qc.json", "w") as f:
        json.dump(qc, f, indent=2)

    return passed

if __name__ == "__main__":
    sample   = sys.argv[1]
    outdir   = sys.argv[2]
    platform = sys.argv[3] if len(sys.argv) > 3 else "illumina"
    qc_check(sample, outdir, platform)
