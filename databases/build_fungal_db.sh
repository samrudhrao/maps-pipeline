#!/bin/bash
# Build MAPS custom fungal Kraken2 database
# Usage: bash build_fungal_db.sh /path/to/output/db
# Requirements: Kraken2, Python3, Biopython

DB_PATH=${1:-~/maps/kraken2_db_fungal_refseq}
SCRIPT_DIR="$(dirname "$0")/../scripts"

echo "Building MAPS fungal database at: $DB_PATH"

# Step 1 — Download taxonomy
kraken2-build --download-taxonomy --db $DB_PATH

# Step 2 — Download fungal genomes from NCBI using accessions
python3 << PYEOF
from Bio import Entrez
import os, subprocess, gzip

Entrez.email = "your.email@institution.edu"  # Replace with your email

accessions = []
with open("$(dirname "$0")/fungal_db_accessions.txt") as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        parts = line.strip().split()
        if len(parts) >= 2:
            accessions.append((parts[0], parts[1]))

os.makedirs(f"$DB_PATH/genomes", exist_ok=True)

for taxid, acc in accessions:
    print(f"Downloading {acc} (taxid {taxid})...")
    # Download from NCBI FTP
    cmd = f"datasets download genome accession {acc} --filename /tmp/{acc}.zip"
    subprocess.run(cmd, shell=True)
    # Extract and tag with taxid
    subprocess.run(f"unzip -o /tmp/{acc}.zip -d /tmp/{acc}/", shell=True)
    # Find the FASTA file and copy to genomes folder
    for root, dirs, files in os.walk(f"/tmp/{acc}/"):
        for file in files:
            if file.endswith(".fna"):
                src = os.path.join(root, file)
                dst = f"$DB_PATH/genomes/{taxid}_{acc}.fna.gz"
                with open(src, 'rb') as f_in, gzip.open(dst, 'wb') as f_out:
                    f_out.write(f_in.read())
                print(f"  Saved: {dst}")
PYEOF

# Step 3 — Add genomes to database
echo "Adding genomes to database..."
for gz in $DB_PATH/genomes/*.fna.gz; do
    taxid=$(basename $gz | cut -d_ -f1)
    tmp=/tmp/genome_${taxid}.fna
    gunzip -c $gz | awk -v t="$taxid" '
        /^>/{acc=$1; sub(/^>/,"",acc); print ">" acc "|kraken:taxid|" t; next}
        {print}
    ' > $tmp
    kraken2-build --add-to-library $tmp --db $DB_PATH
done

# Step 4 — Build database
echo "Building Kraken2 database..."
kraken2-build --build --db $DB_PATH --threads 8

echo "Done — fungal database built at $DB_PATH"
