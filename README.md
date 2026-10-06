# MAPS — Metagenomic Abundance and Pathogen Scoring

Clinical metagenomics pipeline for pathogen detection from raw NGS sequencing data.
Developed at Baylor Scott & White Health.

## Overview
MAPS identifies bacterial, fungal, and viral pathogens directly from FASTQ sequencing files using k-mer classification, genome assembly, AI scoring, and resistance gene detection.

## Supported Platforms
- Illumina (paired-end)
- Ion Torrent (single-end)
- PacBio HiFi (BAM auto-converted)

## Pipeline Steps
16-step pipeline including FastQC, fastp, human depletion, Kraken2 classification, SPAdes assembly, MLST, AMR detection, AI ensemble scoring (RF + XGBoost + SHAP), and clinical report generation.

## Usage
```bash
bash full_pipeline.sh SAMPLE R1 R2 OUTDIR [SPECIMEN_TYPE]
```

## Requirements
- Kraken2
- fastp
- SPAdes
- ABRicate
- MLST
- MultiQC
- Python 3.10+
- See scripts/ for Python dependencies

## Web Interface
Streamlit app at http://localhost:8501
```bash
streamlit run scripts/maps_app.py --server.port 8501
```

## Author
Samrudh Rao, Scientist I
Baylor Scott & White Health
