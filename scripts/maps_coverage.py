#!/usr/bin/env python3
"""
MAPS Coverage Breadth - maps_coverage.py
Step 9b of the MAPS pipeline — genome coverage breadth as informational flag.

IMPORTANT: This script does NOT filter or remove any organisms from results.
It is purely informational — adds breadth percentage alongside abundance.
Low breadth may indicate noise rather than true infection but does NOT
automatically disqualify a detection.

How it works:
  1. For each organism in abundance CSV, looks up RefSeq reference genome
  2. Downloads reference genome if not cached (cached in /maps/ref_genomes/)
  3. Aligns human-depleted reads to reference using Minimap2
  4. Calculates % of reference genome covered by at least 1 read
  5. Reports breadth alongside abundance for physician interpretation

Expected breadth values in clinical metagenomics (short-read, 200-500k reads):
  - Real dominant infection: 0.1-5% breadth
  - Low abundance real infection: 0.05-1% breadth
  - Noise/contamination: <0.05% breadth
  These thresholds are based on BSW validation data — may vary by sample type

To add new organisms: add entry to REFSEQ_MAP with NCBI RefSeq accession
Reference genomes cached in: /maps/ref_genomes/
Requires: minimap2, samtools, NCBI datasets CLI

Usage: python3 maps_coverage.py SAMPLE OUTDIR R1 [R2]
Output: OUTDIR/SAMPLE_coverage.json
"""

import os, sys, subprocess, csv, json, tempfile, shutil

SAMPLE   = sys.argv[1]
OUTDIR   = sys.argv[2]
READS_R1 = sys.argv[3]
READS_R2 = sys.argv[4] if len(sys.argv) > 4 else None
THREADS  = 8
CACHE_DIR = os.path.expanduser("/maps/ref_genomes")
os.makedirs(CACHE_DIR, exist_ok=True)

# RefSeq accessions for common clinical organisms
# Used to download reference genomes for coverage breadth calculation
# To add a new organism: find its RefSeq assembly accession on NCBI
# and add it here as 'organism name (lowercase)': 'GCF_XXXXXXXXX.X'
REFSEQ_MAP = {
    "staphylococcus aureus":       "GCF_000013425.1",
    "staphylococcus epidermidis":  "GCF_000007645.1",
    "staphylococcus cohnii":       "GCF_002811215.1",
    "pseudomonas aeruginosa":      "GCF_000006765.1",
    "escherichia coli":            "GCF_000005845.2",
    "klebsiella pneumoniae":       "GCF_000240185.1",
    "citrobacter freundii":        "GCF_000648755.1",
    "enterococcus faecium":        "GCF_000174395.2",
    "streptococcus pneumoniae":    "GCF_000006885.1",
    "haemophilus influenzae":      "GCF_000027305.1",
    "acinetobacter baumannii":     "GCF_000018445.1",
    "bacillus licheniformis":      "GCF_000009125.1",
    "candida albicans":            "GCF_000182965.3",
    "aspergillus fumigatus":       "GCF_000002655.1",
    "aspergillus flavus":          "GCF_009017415.1",
    "aspergillus terreus":         "GCF_000149645.1",
    "histoplasma capsulatum":      "GCF_000632125.1",
    "coccidioides posadasii":      "GCF_000149335.2",
    "fusarium oxysporum":          "GCF_000149955.1",
    "pneumocystis jirovecii":      "GCF_001477535.1",
    "scedosporium apiospermum":    "GCF_028858705.1",
    "lomentospora prolificans":    "GCF_011037405.1",
    "rhizopus delemar":            "GCF_000402035.1",
    "talaromyces marneffei":       "GCF_000497645.1",
    "mucor circinelloides":        "GCF_000412695.1",
}

def get_genome(org_name):
    """
    Get path to reference genome for an organism.
    Downloads from NCBI if not already cached in /maps/ref_genomes/
    Returns path to gzipped FASTA file or None if organism not in REFSEQ_MAP.
    Genomes are cached after first download — no re-downloading on subsequent runs.
    """
    key = org_name.lower().strip()
    accession = None
    for name, acc in REFSEQ_MAP.items():
        if name in key or all(w in key for w in name.split()[:2]):
            accession = acc
            break
    if not accession:
        return None

    genome_path = os.path.join(CACHE_DIR, f"{accession}.fna.gz")
    if os.path.exists(genome_path):
        return genome_path

    print(f"    Downloading reference: {org_name} ({accession})...")
    try:
        subprocess.run([
            "datasets", "download", "genome", "accession", accession,
            "--filename", f"/tmp/{accession}.zip", "--include", "genome"
        ], capture_output=True, timeout=120, check=True)

        subprocess.run(["unzip", "-o", f"/tmp/{accession}.zip",
                       "-d", f"/tmp/{accession}/"],
                      capture_output=True)

        for root, dirs, files in os.walk(f"/tmp/{accession}/"):
            for fname in files:
                if fname.endswith(".fna"):
                    src = os.path.join(root, fname)
                    with open(src, 'rb') as fin, \
                         open(genome_path, 'wb') as fout:
                        subprocess.run(["gzip", "-c", src], stdout=fout)
                    return genome_path
    except Exception as e:
        print(f"    Warning: Could not download {org_name}: {e}")
    return None

def calc_breadth(genome_path, r1, r2):
    """
    Calculate genome coverage breadth using Minimap2 and samtools.
    Aligns human-depleted reads to reference genome.
    Returns percentage of reference genome covered by at least 1 read.
    Uses temporary directory cleaned up after each calculation.
    Returns None if alignment fails.
    """
    tmp = tempfile.mkdtemp()
    try:
        sam = os.path.join(tmp, "aln.sam")
        bam = os.path.join(tmp, "aln.bam")
        bam_s = bam + ".sorted.bam"

        # Align reads to reference genome using Minimap2
        # -ax sr = short-read preset for Illumina/Ion Torrent reads
        # --secondary=no = exclude secondary alignments to avoid overcounting
        cmd = ["minimap2", "-ax", "sr", "-t", str(THREADS),
               "--secondary=no", genome_path, r1]
        if r2 and r2 != r1 and os.path.exists(r2):
            cmd.append(r2)

        with open(sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f,
                                  stderr=subprocess.PIPE, timeout=300)
        if result.returncode != 0:
            return None

        # Convert SAM to sorted BAM, filter unmapped reads (-F 4)
        subprocess.run(["samtools", "view", "-bS", "-F", "4",
                       sam, "-o", bam], capture_output=True)
        subprocess.run(["samtools", "sort", bam, "-o", bam_s],
                      capture_output=True)
        subprocess.run(["samtools", "index", bam_s], capture_output=True)

        # Calculate per-contig coverage statistics using samtools coverage
        # Column 5 = coverage breadth (% of bases covered)
        res = subprocess.run(["samtools", "coverage", bam_s],
                           capture_output=True, text=True)
        if res.returncode != 0:
            return None

        total = covered = 0
        for line in res.stdout.strip().split('\n'):
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) >= 6:
                try:
                    endpos = int(parts[2])
                    breadth = float(parts[5])
                    total += endpos
                    if breadth > 0:
                        covered += endpos * (breadth / 100)
                except:
                    pass

        return round((covered / total) * 100, 3) if total > 0 else 0.0

    except Exception as e:
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

def main():
    """
    Main function — processes all organisms in abundance CSV.
    Runs coverage breadth calculation for each organism with a reference genome.
    Organisms without reference genomes reported as N/A.
    Results saved to JSON file for physician reference.
    """
    abundance_path = f"{OUTDIR}/{SAMPLE}_abundance.csv"
    if not os.path.exists(abundance_path):
        sys.exit(0)

    with open(abundance_path) as f:
        organisms = list(csv.DictReader(f))

    if not organisms:
        sys.exit(0)

    # Get path to human-depleted reads from Step 3
    # Tries paired-end first, falls back to single-end for Ion Torrent
    r1 = f"{OUTDIR}/human_depleted/{SAMPLE}_clean_1.fastq"
    r2 = f"{OUTDIR}/human_depleted/{SAMPLE}_clean_2.fastq"
    if not os.path.exists(r1):
        r1 = f"{OUTDIR}/human_depleted/{SAMPLE}_clean.fastq"
        r2 = None
    if not os.path.exists(r1):
        print("  Coverage breadth: reads not found — skipping")
        sys.exit(0)

    print(f"\n  Coverage breadth (informational — does not filter results)")
    print(f"  {'Organism':<45} {'Abundance':>10} {'Breadth':>10}")
    print(f"  {'─'*67}")

    coverage_results = {}

    for row in organisms:
        org = row['organism']
        pct = float(row['relative_abundance_%'])
        clean = org.replace('[BACT] ','').replace('[FUNG] ','').lower()

        genome = get_genome(clean)
        if genome is None:
            coverage_results[org] = None
            print(f"  {org:<45} {pct:>9.1f}%  {'N/A':>9} (no reference)")
            continue

        breadth = calc_breadth(genome, r1, r2)
        coverage_results[org] = breadth

        if breadth is None:
            print(f"  {org:<45} {pct:>9.1f}%  {'ERROR':>9}")
        else:
            flag = " ⚠ LOW" if breadth < 1.0 else ""
            print(f"  {org:<45} {pct:>9.1f}%  {breadth:>8.2f}%{flag}")

    # Save coverage breadth results to JSON
    # Includes note that this is informational only
    # Physicians can view breadth alongside abundance in the report
    with open(f"{OUTDIR}/{SAMPLE}_coverage.json", 'w') as f:
        json.dump({
            'sample': SAMPLE,
            'note': 'Informational only — does not affect organism reporting',
            'coverage_breadth': coverage_results
        }, f, indent=2)

    print(f"\n  Coverage saved: {SAMPLE}_coverage.json")

if __name__ == "__main__":
    main()
