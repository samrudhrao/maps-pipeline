#!/usr/bin/env python3
"""
MAPS AMR Resistance Report - maps_amr_report.py
Parses ABRicate output from CARD and ResFinder databases to identify
antimicrobial resistance genes in assembled contigs.

Detection thresholds:
  - Standard genes: >=80% coverage, >=90% identity
  - Critical genes (mecA, vanA, NDM, KPC, MCR, OXA-48, VIM, IMP, GES, SPM): >=60% coverage

Features:
  - 37 resistance gene families with clinical interpretation
  - 8 high-risk combination alerts (e.g. NDM + OXA-48 = pan-resistant)
  - Contig-level Kraken2 carrier identification
  - Drug class lookup for clinical reporting

To add new resistance genes: update RESISTANCE_MAP dictionary
To add new critical genes: update CRITICAL_GENES list
To add new combination alerts: update HIGH_RISK_COMBOS list

Usage: python3 maps_amr_report.py SAMPLE OUTDIR
Output: printed to stdout, saved to OUTDIR/SAMPLE_amr_report.txt
"""

import pandas as pd
import os
import sys
import json

RESISTANCE_MAP = {
    "meca":   {"drug":"Methicillin/Oxacillin","significance":"MRSA"},
    "vana":   {"drug":"Vancomycin","significance":"VRE"},
    "vanb":   {"drug":"Vancomycin","significance":"VRE"},
    "ctx-m":  {"drug":"Cephalosporins","significance":"ESBL"},
    "shv":    {"drug":"Beta-lactams","significance":"ESBL"},
    "tem":    {"drug":"Penicillins","significance":"Beta-lactamase"},
    "ndm":    {"drug":"Carbapenems","significance":"CRE"},
    "kpc":    {"drug":"Carbapenems","significance":"CRE"},
    "oxa-48": {"drug":"Carbapenems","significance":"CRE"},
    "vim":    {"drug":"Carbapenems","significance":"Metallo-beta-lactamase"},
    "erm":    {"drug":"Macrolides","significance":"Macrolide resistance"},
    "tet":    {"drug":"Tetracyclines","significance":"Tetracycline resistance"},
    "aac":    {"drug":"Aminoglycosides","significance":"Aminoglycoside resistance"},
    "sul":    {"drug":"Sulfonamides","significance":"Sulfonamide resistance"},
    "dhfr":   {"drug":"Trimethoprim","significance":"Trimethoprim resistance"},
    "cfr":    {"drug":"Linezolid","significance":"Linezolid resistance"},
    "mcr":    {"drug":"Colistin","significance":"Colistin resistance"},
    "bla":    {"drug":"Beta-lactams","significance":"Beta-lactamase"},
    "dfra":   {"drug":"Trimethoprim","significance":"Trimethoprim resistance"},
    "dfr":    {"drug":"Trimethoprim","significance":"Trimethoprim resistance"},
    "imp":    {"drug":"Carbapenems","significance":"Metallo-beta-lactamase"},
    "ges":    {"drug":"Carbapenems","significance":"Carbapenemase"},
    "spm":    {"drug":"Carbapenems","significance":"Metallo-beta-lactamase"},
    "van":    {"drug":"Vancomycin","significance":"Glycopeptide resistance"},
    "aph":    {"drug":"Aminoglycosides","significance":"Aminoglycoside resistance"},
    "ant":    {"drug":"Aminoglycosides","significance":"Aminoglycoside resistance"},
    "cat":    {"drug":"Chloramphenicol","significance":"Chloramphenicol resistance"},
    "cml":    {"drug":"Chloramphenicol","significance":"Chloramphenicol resistance"},
    "qnr":    {"drug":"Fluoroquinolones","significance":"Quinolone resistance"},
    "gyra":   {"drug":"Fluoroquinolones","significance":"Quinolone resistance"},
    "pbp":    {"drug":"Beta-lactams","significance":"Penicillin-binding protein alteration"},
    "lnu":    {"drug":"Lincosamides","significance":"Lincosamide resistance"},
    "mph":    {"drug":"Macrolides","significance":"Macrolide resistance"},
    "msr":    {"drug":"Macrolides","significance":"Macrolide resistance"},
    "qnr":    {"drug":"Fluoroquinolones","significance":"Quinolone resistance"},
}

# Gene to possible carrier organisms
# Used to suggest which organism likely carries the resistance gene
# when contig-level Kraken2 classification is unresolved
# Update this when new resistance gene-organism associations are identified
GENE_ORGANISM_MAP = {
    "meca":    ["staphylococcus aureus"],
    "vana":    ["enterococcus faecium","enterococcus faecalis"],
    "vanb":    ["enterococcus faecium","enterococcus faecalis"],
    "ctx-m":   ["escherichia coli","klebsiella pneumoniae"],
    "ndm":     ["klebsiella pneumoniae","escherichia coli"],
    "kpc":     ["klebsiella pneumoniae"],
    "oxa-48":  ["klebsiella pneumoniae","escherichia coli"],
    "vim":     ["pseudomonas aeruginosa","klebsiella pneumoniae"],
    "tem":     ["escherichia coli","klebsiella pneumoniae"],
    "shv":     ["klebsiella pneumoniae"],
    "dfra":    ["escherichia coli"],
    "dfr":     ["escherichia coli"],
    "tet":     ["escherichia coli","klebsiella pneumoniae","staphylococcus aureus"],
    "erm":     ["staphylococcus aureus","streptococcus"],
    "mcr":     ["escherichia coli","klebsiella pneumoniae"],
    "cfr":     ["staphylococcus aureus"],
    "aac":     ["pseudomonas aeruginosa","escherichia coli"],
    "qnr":     ["escherichia coli","klebsiella pneumoniae"],
    "bla":     ["escherichia coli","klebsiella pneumoniae","pseudomonas aeruginosa"],
}

def get_significance(gene_name):
    """
    Look up drug class and clinical significance for a resistance gene.
    Matches gene name against RESISTANCE_MAP using partial string matching.
    Returns dict with 'drug' and 'significance' keys.
    Returns 'Unknown' if gene not in lookup table.
    """
    gene_lower = gene_name.lower()
    for key, info in RESISTANCE_MAP.items():
        if key.lower() in gene_lower:
            return info
    return {"drug":"Unknown","significance":"Unknown resistance gene"}

def get_detected_organisms(outdir, sample):
    """Get list of organisms detected above threshold in this sample.
    Used to cross-reference AMR gene carriers against detected pathogens.
    Reads from abundance CSV generated in Step 9.
    """
    detected = []
    abundance_csv = f"{outdir}/{sample}_abundance.csv"
    if os.path.exists(abundance_csv):
        try:
            df = pd.read_csv(abundance_csv)
            df = df[df["relative_abundance_%"] >= 3.0]
            detected = [o.lower().replace("[bact] ","").replace("[fung] ","")
                       for o in df["organism"].tolist()]
        except:
            pass
    return detected

# Critical genes — lower coverage threshold (60%) due to clinical importance
# These genes confer resistance to last-resort antibiotics
# Even partial detection (60% coverage) is clinically actionable
# To add a new critical gene: add its name (lowercase) to this list
CRITICAL_GENES = ["meca", "vana", "vanb", "ndm", "kpc", "mcr", "oxa-48", "vim", "imp", "ges", "spm"]

# High risk gene combinations that trigger automatic clinical alerts
# Format: ([list of genes that must ALL be present], alert message)
# To add a new combination: append a tuple to this list
# These alerts are printed prominently in the AMR report
HIGH_RISK_COMBOS = [
    (["ndm", "oxa-48"],  "PAN-RESISTANT WARNING — NDM + OXA-48 detected — extremely limited treatment options"),
    (["ndm", "kpc"],     "PAN-RESISTANT WARNING — NDM + KPC detected — carbapenem resistance confirmed"),
    (["meca", "cfr"],    "CRITICAL — MRSA + Linezolid resistance — very few treatment options remaining"),
    (["meca", "vana"],   "CRITICAL — MRSA + VRE genes — multi-drug resistant organism"),
    (["mcr", "ndm"],     "CRITICAL — Colistin + Carbapenem resistance — pan-resistant potential"),
    (["kpc", "vim"],     "CRITICAL — Dual carbapenemase detected — pan-resistant potential"),
    (["ndm", "vim"],     "CRITICAL — Dual carbapenemase detected — pan-resistant potential"),
    (["ndm", "imp"],     "CRITICAL — Dual carbapenemase detected — pan-resistant potential"),
    (["meca", "vanb"],   "CRITICAL — MRSA + VRE genes — multi-drug resistant organism"),
    (["mcr", "kpc"],     "CRITICAL — Colistin + Carbapenem resistance — pan-resistant potential"),
]

def check_high_risk_combos(detected_genes):
    """
    Check if detected resistance genes form any high-risk combinations.
    Returns list of alert messages for any combinations found.
    Called after all genes are collected from CARD and ResFinder results.
    """
    """Check if detected genes form any high risk combinations."""
    detected_lower = [g.lower() for g in detected_genes]
    alerts = []
    for combo, message in HIGH_RISK_COMBOS:
        if all(any(c in g for g in detected_lower) for c in combo):
            alerts.append(message)
    return alerts

def is_critical_gene(gene_name):
    """
    Check if a gene name matches any critical gene in CRITICAL_GENES list.
    Critical genes use a lower 60% coverage threshold instead of 80%.
    """
    gene_lower = gene_name.lower()
    return any(c in gene_lower for c in CRITICAL_GENES)

def classify_contig(contig_name, outdir, sample):
    """
    Run Kraken2 on the specific contig carrying the resistance gene.
    Identifies which organism is likely carrying the AMR gene.
    Returns organism name at species, genus, or class level.
    Falls back to 'Unresolved' if Kraken2 cannot classify the contig.
    """
    import subprocess
    import tempfile

    scaffolds = f"{outdir}/assembly/{sample}/scaffolds.fasta"
    if not os.path.exists(scaffolds):
        return None, None

    # Extract the specific contig sequence from the SPAdes assembly
    # Contig name comes from ABRicate output (e.g. NODE_1_length_50000_cov_10)
    contig_seq = []
    found = False
    with open(scaffolds) as f:
        for line in f:
            if line.startswith(">") and contig_name in line:
                found = True
                contig_seq.append(line)
            elif line.startswith(">") and found:
                break
            elif found:
                contig_seq.append(line)

    if not contig_seq:
        return None, None

    # Write contig to temporary FASTA file and run Kraken2 classification
    # Uses standard bacterial DB to identify the carrier organism
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta',
                                     delete=False) as tmp:
        tmp.write("".join(contig_seq))
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["kraken2", "--db", os.path.expanduser(f"{MAPS_HOME}/kraken2_db"),
             "--output", "/dev/stdout", tmp_path],
            capture_output=True, text=True, timeout=30
        )
        os.unlink(tmp_path)

        # Parse Kraken2 output — C = classified, U = unclassified
        # Returns taxid and contig name for downstream name resolution
        for line in result.stdout.strip().split("\n"):
            if line.startswith("C"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    return parts[2], parts[1]
        return "0", contig_name
    except:
        return None, None

# Taxid to organism name lookup for common resistance gene carriers
# Used when Kraken2 classifies a contig at genus/family/class level
# Add new taxids here if carrier identification fails for common organisms
TAXID_NAMES = {
    "1236":   ("Gammaproteobacteria", "class"),
    "91347":  ("Enterobacterales", "order"),
    "543":    ("Enterobacteriaceae", "family"),
    "1224":   ("Proteobacteria", "phylum"),
    "286":    ("Pseudomonas", "genus"),
    "287":    ("Pseudomonas aeruginosa", "species"),
    "562":    ("Escherichia coli", "species"),
    "544":    ("Citrobacter", "genus"),
    "546":    ("Citrobacter freundii", "species"),
    "573":    ("Klebsiella pneumoniae", "species"),
    "1313":   ("Streptococcus pneumoniae", "species"),
    "1280":   ("Staphylococcus aureus", "species"),
    "2":      ("Bacteria", "kingdom"),
    "0":      ("Unclassified", "unknown"),
}

def get_likely_organism(gene_name, detected_organisms, contig_name=None,
                        outdir=None, sample=None):
    """
    Use contig-level Kraken2 classification to identify carrier organism.
    Falls back to gene-organism lookup if contig classification fails.
    """
    # Try contig-level classification first — most accurate method
    # Extracts the contig carrying the AMR gene and runs Kraken2 on it directly
    if contig_name and outdir and sample:
        taxid, _ = classify_contig(contig_name, outdir, sample)
        if taxid and taxid in TAXID_NAMES:
            name, rank = TAXID_NAMES[taxid]
            if rank == "species":
                return name, "contig-level (species)"
            elif rank == "genus":
                return f"{name} sp.", "contig-level (genus)"
            else:
                return f"Unresolved — {name} ({rank} level)", "contig-level (ambiguous)"

    # Fall back to gene-organism lookup if contig classification fails
    # Uses GENE_ORGANISM_MAP to suggest likely carriers based on known associations
    gene_lower = gene_name.lower()
    possible = []
    for key, organisms in GENE_ORGANISM_MAP.items():
        if key in gene_lower:
            possible = organisms
            break

    if not possible:
        return "Unknown", "inferred"

    for org in possible:
        for detected in detected_organisms:
            if org in detected:
                return org.title(), "inferred (detected in sample)"

    return " / ".join(o.title() for o in possible), "inferred (not detected)"

def parse_amr_report(sample, outdir):
    """
    Parse ABRicate resistance results from CARD and ResFinder databases.
    Applies coverage and identity thresholds, deduplicates genes,
    identifies carrier organisms, and checks for high-risk combinations.
    Returns list of resistance gene dictionaries for clinical reporting.
    """
    results  = []
    detected = get_detected_organisms(outdir, sample)
    res_path = f"{outdir}/virulence/{sample}_resistance.txt"

    if os.path.exists(res_path):
        try:
            seen = set()
            with open(res_path) as f:
                for line in f:
                    if line.startswith("#") or not line.strip(): continue
                    parts = line.strip().split("\t")
                    if len(parts) < 10: continue
                    gene       = parts[5].strip()
                    cov        = float(parts[9])
                    ident      = float(parts[10])
                    resistance = parts[14].strip() if len(parts) > 14 else ""
                    gene_base  = gene.rstrip("_0123456789")
                    # Use lower 60% coverage threshold for critical genes
                    # Standard threshold is 80% coverage for non-critical genes
                    cov_threshold = 60 if is_critical_gene(gene) else 80
                    critical      = is_critical_gene(gene)

                    if cov >= cov_threshold and ident >= 90 and gene_base not in seen:
                        seen.add(gene_base)
                        sig      = get_significance(gene)
                        contig   = parts[1].strip() if len(parts) > 1 else None
                        org, status = get_likely_organism(
                            gene, detected, contig, outdir, sample)
                        results.append({
                            "gene":         gene,
                            "coverage":     round(cov, 1),
                            "identity":     round(ident, 1),
                            "drug":         resistance if resistance else sig["drug"],
                            "significance": sig["significance"],
                            "carrier":      org,
                            "carrier_status": status,
                            "critical":     critical,
                        })
        except Exception as e:
            # Log parse errors but continue — don't crash on malformed ABRicate output
            print(f"  Parse error: {e}")

    print(f"\n{'='*70}")
    print(f"AMR RESISTANCE REPORT — {sample}")
    print(f"{'='*70}")

    if not results:
        print("  No resistance genes detected above threshold")
        print("  (>=80% coverage standard, >=60% for critical genes, >=90% identity)")
    else:
        print(f"\n  {'Gene':<20} {'Drug Class':<25} {'Cov':>5} {'ID':>5}  {'Carrier':<35} {'Status'}")
        print("  " + "─"*100)
        for r in results:
            status_icon = "✓ CONFIRMED" if r["carrier_status"] == "confirmed" else "~ INFERRED"
            critical_flag = " ⚠ CRITICAL GENE" if r.get("critical") else ""
            print(f"  {r['gene']:<20} {r['drug']:<25} "
                  f"{r['coverage']:>4.1f}% {r['identity']:>4.1f}%  "
                  f"{r['carrier']:<35} {status_icon}{critical_flag}")
            print(f"  {'':20} {r['significance']}")
            print()

        # Check for high-risk gene combinations that trigger clinical alerts
        # e.g. NDM + OXA-48 = pan-resistant, MRSA + VRE = multi-drug resistant
        detected_genes = [r['gene'] for r in results]
        alerts = check_high_risk_combos(detected_genes)
        if alerts:
            print(f"  {'='*65}")
            for alert in alerts:
                print(f"  🚨 {alert}")
            print(f"  {'='*65}")

    print(f"  Total resistance genes detected: {len(results)}")
    print(f"{'='*70}")

    # Save AMR results to CSV for downstream use in final clinical report
    if results:
        pd.DataFrame(results).to_csv(f"{outdir}/{sample}_amr_report.csv", index=False)

    # Update summary JSON with AMR findings for inclusion in final report
    json_path = f"{outdir}/{sample}_summary.json"
    if os.path.exists(json_path):
        with open(json_path) as f:
            summary = json.load(f)
        summary["amr_genes_detected"] = len(results)
        summary["amr_results"] = results[:5]
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)

    return results

if __name__ == "__main__":
    parse_amr_report(sys.argv[1], sys.argv[2])
