# MAPS — Metagenomic Abundance and Pathogen Scoring
FROM ubuntu:22.04

LABEL maintainer="Samrudh Rao <Samrudh.Rao@BSWHealth.org>"
LABEL description="MAPS Clinical Metagenomics Pipeline"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# ─── System tools + bioinformatics via apt ────────────────────────────
RUN apt-get update && apt-get install -y \
    wget curl git unzip python3 python3-pip \
    fastp samtools kraken2 multiqc minimap2 \
    ncbi-blast+ default-jre perl \
    libdbi-perl libxml-simple-perl \
    any2fasta \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ─── NCBI Datasets CLI ────────────────────────────────────────────────
RUN wget -q https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets \
    && chmod +x datasets && mv datasets /usr/local/bin/

# ─── SPAdes ───────────────────────────────────────────────────────────
RUN wget -q https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz \
    && tar -xzf SPAdes-3.15.5-Linux.tar.gz \
    && mv SPAdes-3.15.5-Linux /opt/spades \
    && rm SPAdes-3.15.5-Linux.tar.gz
ENV PATH="/opt/spades/bin:$PATH"

# ─── ABRicate via apt ─────────────────────────────────────────────────
RUN apt-get update && apt-get install -y abricate \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ─── MLST ─────────────────────────────────────────────────────────────
RUN git clone https://github.com/tseemann/mlst.git /opt/mlst \
    && chmod +x /opt/mlst/bin/mlst
ENV PATH="/opt/mlst/bin:$PATH"

# ─── Python packages ──────────────────────────────────────────────────
RUN pip3 install --no-cache-dir \
    pandas numpy scikit-learn xgboost shap \
    joblib streamlit biopython

# ─── MAPS pipeline ────────────────────────────────────────────────────
WORKDIR /maps
COPY full_pipeline.sh .
COPY scripts/ ./scripts/
COPY models/ ./models/
COPY filter_list.txt .
RUN chmod +x full_pipeline.sh

# ─── Database mount points ────────────────────────────────────────────
RUN mkdir -p /maps/kraken2_db /maps/kraken2_db_fungal_refseq \
    /maps/kraken2_db_extra /maps/kraken2_db_viral \
    /maps/output /maps/ref_genomes /tmp/maps_runs

EXPOSE 8501
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1
CMD ["streamlit", "run", "scripts/maps_app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
