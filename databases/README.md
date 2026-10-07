# MAPS Reference Databases

This folder contains accession lists and build scripts to recreate the MAPS reference databases from publicly available NCBI data.

## Databases

### 1. Standard Bacterial Database
Download directly from Kraken2:
```bash
kraken2-build --standard --db kraken2_db --threads 8
```
Size: ~7.5GB

### 2. Custom Fungal Database (26 species)
```bash
bash build_fungal_db.sh /path/to/kraken2_db_fungal_refseq
```
Accessions listed in `fungal_db_accessions.txt`
Size: ~1.5GB

### 3. Supplemental Database (S. cohnii + B. licheniformis)
```bash
bash build_supplemental_db.sh /path/to/kraken2_db_extra
```
Accessions listed in `supplemental_db_accessions.txt`
Size: ~14MB

### 4. Viral Database
```bash
kraken2-build --download-library viral --db kraken2_db_viral
kraken2-build --build --db kraken2_db_viral --threads 8
```
Size: ~619MB

## Requirements
- Kraken2
- NCBI datasets CLI tool
- Python3 with Biopython
