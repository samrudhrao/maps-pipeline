#!/bin/bash
# Build MAPS supplemental Kraken2 database
# Contains S. cohnii and B. licheniformis
# Usage: bash build_supplemental_db.sh /path/to/output/db

DB_PATH=${1:-~/maps/kraken2_db_extra}

echo "Building MAPS supplemental database at: $DB_PATH"

# Download taxonomy
kraken2-build --download-taxonomy --db $DB_PATH

# Download genomes
mkdir -p $DB_PATH/genomes

ACCESSIONS=(
    "29382 GCA_060225365.1"
    "29382 GCA_060226205.1"
    "29382 GCF_049943095.1"
    "1402 GCA_055397075.1"
    "1402 GCF_052142905.1"
    "1402 GCF_977007745.1"
)

for entry in "${ACCESSIONS[@]}"; do
    taxid=$(echo $entry | cut -d' ' -f1)
    acc=$(echo $entry | cut -d' ' -f2)
    echo "Downloading $acc (taxid $taxid)..."
    datasets download genome accession $acc --filename /tmp/${acc}.zip
    unzip -o /tmp/${acc}.zip -d /tmp/${acc}/
    find /tmp/${acc}/ -name "*.fna" | while read fna; do
        gzip -c $fna > $DB_PATH/genomes/${taxid}_${acc}.fna.gz
    done
done

# Add to library and build
for gz in $DB_PATH/genomes/*.fna.gz; do
    taxid=$(basename $gz | cut -d_ -f1)
    tmp=/tmp/supp_${taxid}.fna
    gunzip -c $gz | awk -v t="$taxid" '
        /^>/{acc=$1; sub(/^>/,"",acc); print ">" acc "|kraken:taxid|" t; next}
        {print}
    ' > $tmp
    kraken2-build --add-to-library $tmp --db $DB_PATH
done

kraken2-build --build --db $DB_PATH --threads 8
echo "Done — supplemental database built at $DB_PATH"
