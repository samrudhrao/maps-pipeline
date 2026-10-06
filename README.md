# MAPS — Metagenomic Abundance and Pathogen Scoring

A clinical metagenomics pipeline for pathogen detection from raw NGS sequencing data, developed at Baylor Scott & White Health.

---

## Overview

MAPS identifies bacterial, fungal, and viral pathogens directly from FASTQ or BAM sequencing files using k-mer classification, de novo genome assembly, AI ensemble scoring, and resistance gene detection. It is designed for clinical use in sterile-site infections where culture may be negative or insufficient.

---

## Key Features

- **Multi-platform support** — Illumina paired-end, Ion Torrent single-end
- **Automated QC** — hard fail criteria stop the pipeline on poor quality samples
- **Specimen-type aware scoring** — adjusts pathogen confidence based on biological context
- **AI ensemble scoring** — Random Forest + XGBoost with SHAP explainability
- **AMR detection** — CARD + ResFinder databases with critical gene alerts and combination flags
- **Viral screening** — filters phages, reports only clinically significant viruses
- **Clinical report generation** — ranked pathogen list with confidence tiers and plain-language summary
- **Web interface** — Streamlit app for non-technical users

---

## Pipeline Overview

| Step | Tool | Purpose |
|---|---|---|
| 1 | FastQC | Raw read quality assessment |
| 2 | fastp | Adapter trimming, quality filtering |
| 3 | Kraken2 | Human read depletion |
| 3b | maps_qc.py | QC pass/fail — pipeline stops on fail |
| 4 | Kraken2 | Bacterial classification |
| 4b | Kraken2 | Supplemental DB (rare organisms) |
| 5 | Kraken2 | Fungal classification (26 species) |
| 5b | Kraken2 | Viral classification |
| 5c | maps_viral_report.py | Clinical viral screening |
| 6 | SPAdes | De novo genome assembly |
| 7 | MLST | Strain typing |
| 8 | ABRicate (VFDB) | Virulence gene detection |
| 8b | ABRicate (CARD + ResFinder) | AMR resistance gene detection |
| 9 | Python | Abundance report |
| 10 | run_ai_scoring.py | AI confidence scoring |
| 11 | Python | Co-infection detection |
| 12 | MultiQC | Run summary |
| 13 | maps_final_report.py | Clinical report generation |
| 14 | maps_predict.py | Random Forest prediction |
| 14b | maps_ensemble.py | RF + XGBoost + SHAP ensemble |
| 15 | maps_database.py | Save to SQLite database |
| 16 | maps_stats.py | Confidence intervals + anomaly detection |

---

## Usage

```bash
bash full_pipeline.sh SAMPLE R1 R2 OUTDIR [SPECIMEN_TYPE]
```

**Example:**
```bash
bash full_pipeline.sh 151 \
    /path/to/151_R1.fastq.gz \
    /path/to/151_R2.fastq.gz \
    /path/to/output \
    abscess
```

**Supported specimen types:** abscess, synovial, pleural, peritoneal, plasma, bal, unknown

---

## QC Thresholds

| Metric | Hard Fail |
|---|---|
| Total reads | < 50,000 |
| Human DNA | > 90% |
| Microbial reads | < 1,000 |
| Q30 rate | < 50% |

---

## AI Scoring

Two machine learning models trained on 26 confirmed clinical samples:

- **Random Forest** — 500 trees, 100% LOO cross-validation accuracy
- **XGBoost** — 500 trees, 96.2% LOO cross-validation accuracy
- **Ensemble** — RF (55%) + XGBoost (45%) weighted average
- **SHAP explainability** — top 3 feature drivers reported per organism

Confidence tiers:
- HIGH ≥ 70%
- MODERATE 40-69%
- LOW < 40%

---

## Detection Databases

| Database | Size | Contents |
|---|---|---|
| Kraken2 standard | 7.5GB | Full bacterial RefSeq |
| Custom fungal | 1.5GB | 26 clinical fungal species |
| Supplemental | 14MB | S. cohnii, B. licheniformis |
| Viral | 619MB | Clinical human viruses |
| CARD | — | 6,052 AMR resistance genes |
| ResFinder | — | 3,206 AMR resistance genes |
| VFDB | — | 4,592 virulence factors |

---

## Web Interface

```bash
streamlit run scripts/maps_app.py --server.port 8501
```

Access at http://localhost:8501

---

## Installation

**1. Clone the repository**
```bash
git clone https://github.com/samrudhrao/maps-pipeline.git
cd maps-pipeline
```

**2. Install system dependencies**

The following tools are required. Installation method may vary by system:
- [Kraken2](https://github.com/DerrickWood/kraken2)
- [fastp](https://github.com/OpenGENOME/fastp)
- [SPAdes](https://github.com/ablab/spades)
- [ABRicate](https://github.com/tseemann/abricate)
- [MLST](https://github.com/tseemann/mlst)
- [MultiQC](https://multiqc.info)
- [samtools](http://www.htslib.org)

**3. Install Python packages**
```bash
pip install pandas numpy scikit-learn xgboost shap joblib streamlit
```

**4. Download reference databases**

Kraken2 standard database (~7.5GB):
```bash
kraken2-build --standard --db kraken2_db --threads 8
```

Custom fungal (26 species), viral and supplemental databases are maintained internally. Contact the maintainer for access.

**5. Run the pipeline**
```bash
bash full_pipeline.sh SAMPLE R1.fastq.gz R2.fastq.gz output/ abscess
```

---

## Requirements

**Tools:**
- Kraken2
- fastp
- SPAdes
- ABRicate
- MLST
- MultiQC
- samtools

**Python packages:**
- pandas, numpy, scikit-learn, xgboost, shap, joblib, streamlit

---

## Institution

Baylor Scott & White Health — Temple, TX
Developer: Samrudh Rao, Scientist I

---

## Status

Currently under clinical validation. Not yet FDA cleared. For research use only.
