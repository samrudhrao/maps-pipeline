# MAPS — Metagenomic Abundance and Pathogen Scoring
# Docker image: Pipeline + Models + Streamlit
# Databases mounted at runtime from Azure Blob Storage

FROM ubuntu:22.04

LABEL maintainer="Samrudh Rao <Samrudh.Rao@BSWHealth.org>"
LABEL description="MAPS Clinical Metagenomics Pipeline"

# Prevent interactive prompts during install
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# ─── System dependencies ───────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    wget curl git unzip \
    python3 python3-pip \
    fastp \
    samtools \
    kraken2 \
    multiqc \
    ncbi-blast+ \
    perl \
    default-jre \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ─── SPAdes ────────────────────────────────────────────────────────────
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz \
    && tar -xzf SPAdes-3.15.5-Linux.tar.gz \
    && mv SPAdes-3.15.5-Linux /opt/spades \
    && rm SPAdes-3.15.5-Linux.tar.gz
ENV PATH="/opt/spades/bin:$PATH"

# ─── ABRicate ──────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y cpanminus \
    && cpanm Bio::Perl List::Util \
    && git clone https://github.com/tseemann/abricate.git /opt/abricate \
    && chmod +x /opt/abricate/bin/abricate
ENV PATH="/opt/abricate/bin:$PATH"

# ─── MLST ──────────────────────────────────────────────────────────────
RUN git clone https://github.com/tseemann/mlst.git /opt/mlst \
    && chmod +x /opt/mlst/bin/mlst
ENV PATH="/opt/mlst/bin:$PATH"

# ─── Python packages ───────────────────────────────────────────────────
RUN pip3 install --no-cache-dir \
    pandas numpy scikit-learn xgboost shap \
    joblib streamlit biopython

# ─── MAPS pipeline ─────────────────────────────────────────────────────
WORKDIR /maps
COPY full_pipeline.sh .
COPY scripts/ ./scripts/
COPY models/ ./models/
COPY filter_list.txt .

RUN chmod +x full_pipeline.sh

# ─── Database mount points ─────────────────────────────────────────────
# Databases are mounted at runtime from Azure Blob Storage
# Mount: /maps/kraken2_db
# Mount: /maps/kraken2_db_fungal_refseq
# Mount: /maps/kraken2_db_extra
# Mount: /maps/kraken2_db_viral
RUN mkdir -p \
    /maps/kraken2_db \
    /maps/kraken2_db_fungal_refseq \
    /maps/kraken2_db_extra \
    /maps/kraken2_db_viral \
    /maps/output \
    /tmp/maps_runs

# ─── Streamlit port ────────────────────────────────────────────────────
EXPOSE 8501

# ─── Healthcheck ───────────────────────────────────────────────────────
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# ─── Start Streamlit ───────────────────────────────────────────────────
CMD ["streamlit", "run", "scripts/maps_app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
