import os
import sys
import json

# Clinically significant human viruses to report
CLINICAL_VIRUSES = {
    "cytomegalovirus":          {"abbrev": "CMV",    "significance": "Opportunistic — immunocompromised risk"},
    "epstein-barr":             {"abbrev": "EBV",    "significance": "Lymphoproliferative disease risk"},
    "herpes simplex":           {"abbrev": "HSV",    "significance": "Encephalitis/disseminated infection risk"},
    "human herpesvirus":        {"abbrev": "HHV",    "significance": "Herpesvirus — clinical review needed"},
    "varicella-zoster":         {"abbrev": "VZV",    "significance": "Reactivation risk"},
    "influenza":                {"abbrev": "Flu",    "significance": "Respiratory pathogen"},
    "respiratory syncytial":    {"abbrev": "RSV",    "significance": "Respiratory pathogen"},
    "sars-cov":                 {"abbrev": "COVID",  "significance": "Respiratory pathogen"},
    "adenovirus":               {"abbrev": "AdV",    "significance": "Respiratory/GI pathogen"},
    "parvovirus":               {"abbrev": "B19",    "significance": "Aplastic crisis risk"},
    "hepatitis":                {"abbrev": "HepV",   "significance": "Hepatitis virus"},
    "hiv":                      {"abbrev": "HIV",    "significance": "Immunodeficiency virus"},
    "human papillomavirus":     {"abbrev": "HPV",    "significance": "Oncogenic virus"},
    "polyomavirus":             {"abbrev": "BKV/JCV","significance": "Nephropathy risk in transplant"},
    "norovirus":                {"abbrev": "NoV",    "significance": "GI pathogen"},
    "enterovirus":              {"abbrev": "EV",     "significance": "Meningitis/encephalitis risk"},
    "coronavirus":              {"abbrev": "CoV",    "significance": "Respiratory pathogen"},
    "metapneumovirus":          {"abbrev": "hMPV",   "significance": "Respiratory pathogen"},
    "parainfluenza":            {"abbrev": "PIV",    "significance": "Respiratory pathogen"},
    "bocavirus":                {"abbrev": "HBoV",   "significance": "Respiratory pathogen"},
}

# Known non-human pathogens to exclude
EXCLUDE_VIRAL = [
    "phage", "bacteriophage", "gihfavirus", "siphoviridae",
    "myoviridae", "podoviridae", "caudovirales", "leviviridae",
    "inovirus", "microvirus", "plant virus", "tobacco",
    "pelohabitans", "baculovirus", "nucleopolyhedrovirus",
]

def parse_viral_report(path, min_reads=50, min_pct=0.1):
    results = []
    if not os.path.exists(path):
        return results
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            pct   = float(parts[0])
            reads = int(parts[1])
            rank  = parts[3].strip()
            name  = parts[5].strip().lower()

            if rank != "S" or reads < min_reads or pct < min_pct:
                continue

            # Skip non-human pathogens
            if any(ex in name for ex in EXCLUDE_VIRAL):
                continue

            # Check if clinically significant
            sig = None
            abbrev = None
            for key, info in CLINICAL_VIRUSES.items():
                if key in name:
                    sig    = info["significance"]
                    abbrev = info["abbrev"]
                    break

            if sig:
                results.append({
                    "virus":        parts[5].strip(),
                    "abbrev":       abbrev,
                    "reads":        reads,
                    "pct":          round(pct, 2),
                    "significance": sig,
                })

    return sorted(results, key=lambda x: -x["reads"])

def viral_report(sample, outdir):
    viral_path = f"{outdir}/bact_results/{sample}_viral_report.txt"
    results    = parse_viral_report(viral_path)

    print(f"\n{'='*65}")
    print(f"VIRAL SCREENING REPORT — {sample}")
    print(f"{'='*65}")

    if not results:
        print("  No clinically significant viruses detected")
        print("  (threshold: >=50 reads, >=0.1% abundance)")
    else:
        print(f"\n  {'Virus':<35} {'Abbrev':<8} {'Reads':>6} {'Pct':>6}  Significance")
        print("  " + "─"*90)
        for r in results:
            print(f"  {r['virus']:<35} {r['abbrev']:<8} "
                  f"{r['reads']:>6} {r['pct']:>5.2f}%  {r['significance']}")
        print(f"\n  ⚠ Viral co-infection detected — clinical review recommended")

    print(f"\n  Total clinically significant viruses: {len(results)}")
    print(f"{'='*65}")

    # Update summary json
    json_path = f"{outdir}/{sample}_summary.json"
    if os.path.exists(json_path):
        with open(json_path) as f:
            summary = json.load(f)
        summary["viral_detected"]  = len(results)
        summary["viral_results"]   = results[:3]
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)

    return results

if __name__ == "__main__":
    viral_report(sys.argv[1], sys.argv[2])
