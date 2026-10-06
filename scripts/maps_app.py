import streamlit as st
import subprocess
import os
import json
import pandas as pd
import shutil
from datetime import datetime

st.set_page_config(
    page_title="MAPS — Metagenomic Abundance and Pathogen Scoring",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS — BSW blue and white
st.markdown("""
<style>
    .main { background-color: #ffffff; }
    .stApp { background-color: #ffffff; }
    h1, h2, h3 { color: #003087; }
    .stButton > button {
        background-color: #003087;
        color: white;
        border-radius: 8px;
        font-weight: bold;
        padding: 0.5rem 2rem;
    }
    .stButton > button:hover { background-color: #0050b3; }
    .stTabs [data-baseweb="tab"] {
        font-size: 1rem;
        font-weight: 600;
        color: #003087;
    }
    .stTabs [aria-selected="true"] {
        background-color: #003087;
        color: white;
        border-radius: 6px 6px 0 0;
    }
    .metric-card {
        background-color: #f0f4ff;
        border-left: 4px solid #003087;
        padding: 1rem;
        border-radius: 6px;
        margin: 0.5rem 0;
    }
    .stExpander { border-left: 3px solid #003087; }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 style='color:#003087;'>🧬 MAPS</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='color:#003087;'>Metagenomic Abundance and Pathogen Scoring</h3>",
            unsafe_allow_html=True)
st.caption("Baylor Scott & White Health | Clinical Metagenomics Pipeline")
st.divider()

# Tabs
tab1, tab2 = st.tabs(["🏠 Run Pipeline", "📖 About MAPS"])

# ─────────────────────────────────────────────
# TAB 1 — RUN PIPELINE
# ─────────────────────────────────────────────
with tab1:

    # Sidebar
    with st.sidebar:
        st.markdown("<h3 style='color:#003087;'>⚙️ Run Settings</h3>", unsafe_allow_html=True)
        sample_name   = st.text_input("Sample Name", placeholder="e.g. 151, T1, RMID-001")
        specimen_type = st.selectbox("Specimen Type",
                                      ["unknown","abscess","synovial","pleural",
                                       "peritoneal","plasma","bal"])
        platform      = st.radio("Sequencing Platform",
                                  ["Auto-detect","Illumina (paired-end)",
                                   "Ion Torrent (single-end)"])
        st.divider()
        st.markdown("<h4 style='color:#003087;'>Pipeline Steps</h4>", unsafe_allow_html=True)
        for s in ["FastQC + fastp trimming","Human DNA depletion",
                  "Bacterial/Fungal/Viral classification","SPAdes genome assembly",
                  "AMR resistance detection","AI ensemble scoring (RF + XGBoost + SHAP)",
                  "Clinical report generation"]:
            st.markdown(f"✓ {s}")

    # Upload section
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("<h3 style='color:#003087;'>📁 Upload FASTQ Files</h3>",
                    unsafe_allow_html=True)
        r1_file = st.file_uploader(
            "R1 FASTQ (or single-end for Ion Torrent)",
            type=["fastq","fastq.gz","fq","fq.gz","bam"], key="r1"
        )
        r2_file = None
        if platform != "Ion Torrent (single-end)":
            r2_file = st.file_uploader(
                "R2 FASTQ (Illumina paired-end)",
                type=["fastq","fastq.gz","fq","fq.gz","bam"], key="r2"
            )

    with col2:
        st.markdown("<h3 style='color:#003087;'>✅ QC Requirements</h3>",
                    unsafe_allow_html=True)
        st.table(pd.DataFrame({
            "Metric":    ["Minimum reads","Maximum human DNA",
                          "Minimum Q30 rate","Minimum microbial reads"],
            "Threshold": ["50,000","90%","50%","1,000"]
        }))

    st.divider()

    STEPS = [
        "FastQC quality check","fastp trimming","Human read depletion",
        "QC validation","Bacterial classification","Supplemental DB",
        "Fungal classification","Viral classification","Viral screening",
        "Genome assembly","MLST typing","Virulence genes","AMR detection",
        "Abundance report","AI confidence scoring","Co-infection detection",
        "MultiQC","Clinical report","AI prediction","Ensemble prediction",
        "Saving to database","Statistical analysis",
    ]

    run_ready = sample_name and r1_file and \
                (r2_file or platform == "Ion Torrent (single-end)")

    if not run_ready:
        if not sample_name:
            st.warning("⚠️ Enter a sample name to continue")
        elif not r1_file:
            st.warning("⚠️ Upload at least the R1 FASTQ file")

    if st.button("🚀 Run MAPS Pipeline", type="primary",
                 disabled=not run_ready, use_container_width=True):

        run_dir = f"/tmp/maps_{sample_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        outdir  = f"{run_dir}/output"
        os.makedirs(outdir, exist_ok=True)

        r1_path = f"{run_dir}/{r1_file.name}"
        with open(r1_path, "wb") as f:
            f.write(r1_file.getbuffer())

        if r2_file:
            r2_path = f"{run_dir}/{r2_file.name}"
            with open(r2_path, "wb") as f:
                f.write(r2_file.getbuffer())
        else:
            r2_path = r1_path

        # Progress UI
        st.markdown("<h3 style='color:#003087;'>🔄 Running Pipeline</h3>",
                    unsafe_allow_html=True)
        progress_bar = st.progress(0)
        current_step = st.empty()
        step_index   = 0
        failed       = False
        fail_reason  = ""

        cmd = [
            "bash", os.path.expanduser("~/maps/full_pipeline.sh"),
            sample_name, r1_path, r2_path, outdir, specimen_type
        ]

        try:
            process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1
            )
            for line in process.stdout:
                line = line.strip()
                if line.startswith("Step"):
                    step_index = min(step_index + 1, len(STEPS))
                    pct = min(int((step_index / len(STEPS)) * 100), 95)
                    progress_bar.progress(pct)
                    current_step.info(f"⏳ {line}")
                if "QC FAILED" in line or "failed QC" in line:
                    failed = True
                    fail_reason = line
            process.wait()
        except Exception as e:
            failed = True
            fail_reason = str(e)

        if failed:
            progress_bar.progress(100)
            current_step.empty()
            st.error(f"❌ Pipeline stopped at QC: {fail_reason}")
            st.warning("""
            Sample failed quality control. Common reasons:
            - Insufficient reads (need ≥ 50,000)
            - Too much human DNA (> 90%)
            - Low Q30 rate (< 50%) — poor sequencing quality
            - Too few microbial reads (< 1,000)
            """)

        elif os.path.exists(f"{outdir}/{sample_name}_summary.json"):
            progress_bar.progress(100)
            current_step.empty()
            st.success("✅ Pipeline complete!")

            with open(f"{outdir}/{sample_name}_summary.json") as f:
                summary = json.load(f)

            st.divider()
            st.markdown("<h2 style='color:#003087;'>📊 Results</h2>",
                        unsafe_allow_html=True)

            # QC metrics
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Quality",         summary.get("quality","N/A"))
            m2.metric("Human DNA",       f"{summary.get('human_pct',0):.1f}%")
            m3.metric("Microbial Reads", f"{summary.get('microbial_reads',0):,}")
            m4.metric("Infection Type",  summary.get("infection_type","N/A"))

            st.info(f"🩺 **Clinical Summary:** {summary.get('clinical_summary','N/A')}")

            # Pathogen table
            abundance_path = f"{outdir}/{sample_name}_abundance.csv"
            if os.path.exists(abundance_path):
                st.markdown("<h3 style='color:#003087;'>🦠 Detected Pathogens</h3>",
                            unsafe_allow_html=True)
                df = pd.read_csv(abundance_path)
                df = df[df["relative_abundance_%"] >= 3.0].sort_values(
                    "relative_abundance_%", ascending=False).reset_index(drop=True)
                df.index += 1
                df.columns = ["Organism","Abundance (%)"]
                df["Abundance (%)"] = df["Abundance (%)"].round(2)
                st.dataframe(df, use_container_width=True)

            # AI scores
            ai_path = f"{outdir}/{sample_name}_ai_scores.csv"
            if os.path.exists(ai_path):
                ai_df = pd.read_csv(ai_path)
                if len(ai_df) > 0:
                    st.markdown("<h3 style='color:#003087;'>🤖 AI Confidence Scores</h3>",
                                unsafe_allow_html=True)
                    st.dataframe(ai_df, use_container_width=True, hide_index=True)

            # AMR
            amr_path = f"{outdir}/{sample_name}_amr_report.csv"
            if os.path.exists(amr_path):
                amr_df = pd.read_csv(amr_path)
                if len(amr_df) > 0:
                    st.markdown("<h3 style='color:#003087;'>⚠️ AMR Resistance Genes</h3>",
                                unsafe_allow_html=True)
                    st.dataframe(amr_df, use_container_width=True, hide_index=True)
                else:
                    st.success("✅ No AMR resistance genes detected")

            # Viral
            viral = summary.get("viral_results",[])
            if viral:
                st.markdown("<h3 style='color:#003087;'>🦠 Viral Detections</h3>",
                            unsafe_allow_html=True)
                st.dataframe(pd.DataFrame(viral), use_container_width=True, hide_index=True)
            else:
                st.success("✅ No clinically significant viruses detected")

            # Anomaly
            stats_path = f"{outdir}/{sample_name}_stats.json"
            if os.path.exists(stats_path):
                with open(stats_path) as f:
                    stats = json.load(f)
                if stats.get("anomaly_status") == "ANOMALY":
                    st.warning(f"⚠️ Anomaly detected — score: {stats.get('anomaly_score','N/A')}")
                else:
                    st.success("✅ Sample profile normal")

            # Downloads
            st.divider()
            st.markdown("<h3 style='color:#003087;'>📥 Download Results</h3>",
                        unsafe_allow_html=True)
            d1, d2, d3 = st.columns(3)
            with d1:
                with open(f"{outdir}/{sample_name}_summary.json") as f:
                    st.download_button("📄 Summary JSON", f.read(),
                                       file_name=f"{sample_name}_summary.json",
                                       mime="application/json",
                                       use_container_width=True)
            with d2:
                if os.path.exists(abundance_path):
                    with open(abundance_path) as f:
                        st.download_button("📊 Abundance CSV", f.read(),
                                           file_name=f"{sample_name}_abundance.csv",
                                           mime="text/csv",
                                           use_container_width=True)
            with d3:
                if os.path.exists(amr_path):
                    with open(amr_path) as f:
                        st.download_button("🧬 AMR Report", f.read(),
                                           file_name=f"{sample_name}_amr_report.csv",
                                           mime="text/csv",
                                           use_container_width=True)
        else:
            st.error("❌ Pipeline failed — no results generated")

        try:
            shutil.rmtree(run_dir)
        except:
            pass

# ─────────────────────────────────────────────
# TAB 2 — ABOUT MAPS
# ─────────────────────────────────────────────
with tab2:
    st.markdown("<h2 style='color:#003087;'>📖 About MAPS</h2>", unsafe_allow_html=True)
    st.markdown("""
    MAPS (Metagenomic Abundance and Pathogen Scoring) is a clinical metagenomics pipeline 
    developed at Baylor Scott & White Health. It identifies bacterial, fungal, and viral 
    pathogens directly from raw sequencing data using a combination of k-mer classification, 
    genome assembly, AI scoring, and resistance gene detection.
    """)
    st.divider()

    with st.expander("1️⃣ Pipeline Overview"):
        st.markdown("""
        MAPS takes raw FASTQ sequencing files as input and produces a ranked list of 
        detected pathogens with AI confidence scores, resistance genes, and a clinical summary.

        **Input:** Raw FASTQ files (Illumina paired-end or Ion Torrent single-end)

        **Output:**
        - Ranked pathogen list with abundance percentages
        - AI confidence scores (HIGH / MODERATE / LOW)
        - AMR resistance gene report
        - Viral screening report
        - Confidence intervals and anomaly detection
        - Downloadable CSV and JSON reports

        **Pipeline command:**
    bash full_pipeline.sh SAMPLE R1 R2 OUTDIR [SPECIMEN_TYPE]

        The pipeline runs 16 sequential steps taking approximately 20-45 minutes per sample 
        depending on read count and available compute resources.
        """)

    with st.expander("2️⃣ Supported Sample Types"):
        st.markdown("""
        MAPS applies specimen-type aware scoring to adjust pathogen confidence based on 
        the biological context of the sample. Organisms that are implausible for a given 
        specimen type are automatically downscored or excluded.
        """)
        st.table(pd.DataFrame({
            "Specimen Type": ["Abscess","Synovial fluid","Pleural fluid",
                              "Peritoneal fluid","Plasma","BAL"],
            "Typical Human DNA": ["5-40%","70-95%","40-80%","20-60%","95-99%","50-85%"],
            "Common pathogens": [
                "S. aureus, Streptococcus, E. coli",
                "S. aureus, S. epidermidis, N. gonorrhoeae",
                "S. pneumoniae, S. aureus, Klebsiella",
                "E. coli, Candida, Klebsiella",
                "E. coli, S. aureus, Candida",
                "P. aeruginosa, Klebsiella, Aspergillus"
            ],
            "Excluded organisms": [
                "Environmental bacteria",
                "Gut anaerobes, Pneumocystis, oral flora",
                "Gut anaerobes, environmental bacteria",
                "Environmental bacteria",
                "Gut anaerobes, environmental bacteria",
                "Gut anaerobes"
            ]
        }))

    with st.expander("3️⃣ Sequencing Platforms"):
        st.markdown("""
        MAPS supports two sequencing platforms with automatic detection based on input files.

        **Illumina (paired-end)**
        - Two FASTQ files (R1 and R2)
        - Fixed 150bp read length
        - Higher Q30 rates (85-95%)
        - Lower error rate — primarily substitution errors
        - fastp: paired-end mode with adapter detection
        - Kraken2: --paired flag
        - SPAdes: standard --careful mode

        **Ion Torrent (single-end)**
        - Single FASTQ file
        - Variable read length (50-400bp)
        - Lower Q30 rates (70-80% typical)
        - Higher indel error rate (homopolymer errors)
        - fastp: single-end mode
        - Kraken2: single-end mode (no --paired flag)
        - SPAdes: --iontorrent --careful mode

        **Platform detection:** If R1 equals R2 (same file passed twice) → Ion Torrent. 
        If R1 and R2 are different files → Illumina.
        """)

    with st.expander("4️⃣ Quality Control (QC)"):
        st.markdown("""
        Every sample passes through QC at Step 3b. The pipeline stops automatically if 
        any hard-fail criterion is met. Results are never reported from a QC-failed sample.
        """)
        st.table(pd.DataFrame({
            "Metric": ["Total reads","Human DNA %","Microbial reads","Q30 rate"],
            "Hard fail threshold": ["< 50,000","> 90%","< 1,000","< 50%"],
            "What it means": [
                "Insufficient sequencing depth for reliable detection",
                "Too much host DNA — almost no microbial signal left",
                "Not enough microbial reads after human depletion",
                "Sequencing quality too poor — reads have >1% error rate per base"
            ]
        }))
        st.markdown("""
        **Informational warnings (not hard fails):**
        - Q30 between 50-70%: results may be less reliable
        - Duplication rate > 30%: possible PCR over-amplification
        - Insert size < 50bp: possible degraded DNA

        **QC metrics source:** fastp JSON output (Q30, duplication, insert size) 
        and Kraken2 human depletion report (total reads, human %, microbial reads).
        """)

    with st.expander("5️⃣ Pipeline Steps"):
        st.markdown("""
        | Step | Tool | Purpose |
        |---|---|---|
        | 1 | FastQC | Raw read quality assessment |
        | 2 | fastp | Adapter trimming, quality filtering (Q20, length ≥50bp) |
        | 3 | Kraken2 (standard DB) | Human read depletion |
        | 3b | maps_qc.py | QC pass/fail — pipeline stops on fail |
        | 4 | Kraken2 (standard DB) | Bacterial classification |
        | 4b | Kraken2 (supplemental DB) | S. cohnii and other rare organisms |
        | 5 | Kraken2 (fungal DB) | Fungal classification — 26 clinical species |
        | 5b | Kraken2 (viral DB) | Viral classification |
        | 5c | maps_viral_report.py | Clinical viral screening — filters phages |
        | 6 | SPAdes | De novo genome assembly from reads |
        | 7 | MLST | Sequence-based strain typing |
        | 8 | ABRicate (VFDB) | Virulence gene detection |
        | 8b | ABRicate (CARD + ResFinder) | AMR resistance gene detection |
        | 9 | Python inline | Abundance calculation and organism ranking |
        | 10 | run_ai_scoring.py | AI confidence scoring |
        | 11 | Python inline | Co-infection detection |
        | 12 | MultiQC | Sequencing run summary report |
        | 13 | maps_final_report.py | Full clinical report generation |
        | 14 | maps_predict.py | Random Forest AI prediction |
        | 14b | maps_ensemble.py | RF + XGBoost ensemble + SHAP |
        | 15 | maps_database.py | Save results to SQLite database |
        | 16 | maps_stats.py | Wilson confidence intervals + anomaly detection |
        """)

    with st.expander("6️⃣ Detection Databases"):
        st.markdown("""
        **Bacterial Database (Kraken2 standard)**
        - Full NCBI RefSeq standard database
        - 7.5GB — covers all major bacterial species
        - Used for human depletion and bacterial classification

        **Fungal Database (custom, 26 species)**
        """)
        st.table(pd.DataFrame({
            "Category": ["Candida","Candida","Candida","Candida","Candida","Candida",
                         "Aspergillus","Aspergillus","Aspergillus","Aspergillus",
                         "Pneumocystis","Endemic fungi","Endemic fungi","Endemic fungi",
                         "Endemic fungi","Mucorales","Mucorales","Mucorales","Mucorales",
                         "Other","Other","Other","Other","Other","Other","Other"],
            "Species": [
                "Candida albicans","Candida tropicalis","Candida auris",
                "Nakaseomyces glabratus (C. glabrata)","Pichia kudriavzevii (C. krusei)",
                "Lodderomyces parapsilosis (C. parapsilosis)",
                "Aspergillus fumigatus","Aspergillus niger","Aspergillus terreus","Aspergillus flavus",
                "Pneumocystis jirovecii",
                "Histoplasma capsulatum","Coccidioides immitis","Coccidioides posadasii",
                "Talaromyces marneffei",
                "Mucor circinelloides","Rhizopus microsporus","Rhizopus delemar","Lichtheimia corymbifera",
                "Cryptococcus neoformans","Sporothrix schenckii","Fusarium solani","Fusarium oxysporum",
                "Scedosporium apiospermum","Lomentospora prolificans","Alternaria alternata"
            ]
        }))
        st.markdown("""
        **Viral Database**
        - Clinical viruses: CMV, EBV, HSV, VZV, influenza, RSV, adenovirus, and others
        - Phages and non-human pathogens automatically filtered

        **Supplemental Database (rare organisms)**
        - Staphylococcus cohnii — poorly represented in standard DB
        - Bacillus licheniformis
        - Triggers only when ≥500 Staphylococcus reads detected

        **AMR Databases**
        - CARD: 6,052 resistance gene sequences
        - ResFinder: 3,206 resistance gene sequences
        - Detection threshold: ≥80% coverage, ≥90% identity (≥60% for critical genes)

        **Virulence Database**
        - VFDB: 4,592 virulence factor sequences
        """)

    with st.expander("7️⃣ AI Scoring"):
        st.markdown("""
        MAPS uses an ensemble of two machine learning models to score the likelihood 
        that each detected organism is a true pathogen rather than background noise or contamination.

        **Random Forest (RF)**
        - 500 decision trees trained on 26 confirmed clinical samples
        - Each tree votes — majority determines the confidence score
        - 100% leave-one-out cross-validation accuracy on training data
        - Weight in ensemble: 55%

        **XGBoost**
        - Sequential tree boosting — each tree fixes errors of the previous
        - 96.2% leave-one-out cross-validation accuracy
        - Better at edge cases and borderline detections
        - Weight in ensemble: 45%

        **11 Features used by both models:**
        """)
        st.table(pd.DataFrame({
            "Feature": ["abundance","human_pct","microbial_reads","pathogen_score",
                        "above_threshold","high_human","abundance_x_score",
                        "log_reads","low_reads","is_fungal","sample_type_enc"],
            "Description": [
                "Relative abundance % of organism",
                "Human DNA % in sample",
                "Total microbial read count",
                "Known pathogen score (0-1) from lookup table",
                "Whether organism is above 3% threshold",
                "Whether human DNA exceeds 50%",
                "Abundance × pathogen score combined",
                "Log-transformed microbial read count",
                "Whether microbial reads < 5,000",
                "Whether organism is fungal",
                "Encoded specimen type"
            ]
        }))
        st.markdown("""
        **Confidence Tiers:**
        - 🟢 **HIGH** — ensemble score ≥ 70% — likely true pathogen
        - 🟡 **MODERATE** — ensemble score 40-69% — interpret with caution
        - 🔴 **LOW** — ensemble score < 40% — likely background or contaminant

        **Species vs Genus level reporting:**
        - Score ≥ 70%: full species name reported (e.g. *Staphylococcus aureus*)
        - Score < 70%: genus level only (e.g. *Staphylococcus* sp.)

        **SHAP Explainability:**
        Each AI score is accompanied by the top 3 features that drove the prediction.
        For example: `abundance(0.788) | above_threshold(0.225) | abundance_x_score(0.110)`
        means abundance was the dominant factor pushing confidence up by 0.788.

        **Specimen-type scoring:**
        Pathogen scores are adjusted based on specimen type. Organisms that are 
        biologically implausible in a given specimen type receive lower scores.
        """)

    with st.expander("8️⃣ AMR Resistance Report"):
        st.markdown("""
        Resistance genes are detected by running ABRicate against CARD and ResFinder 
        databases on the assembled contigs from SPAdes.

        **Detection thresholds:**
        - Standard genes: ≥80% coverage, ≥90% identity
        - Critical genes (mecA, vanA, vanB, NDM, KPC, MCR, VIM, OXA-48): ≥60% coverage

        **Critical genes** are flagged with ⚠️ CRITICAL GENE and trigger at lower 
        coverage thresholds because even partial detection of these genes is clinically significant.

        **High-risk combination alerts** — the pipeline automatically flags dangerous combinations:
        """)
        st.table(pd.DataFrame({
            "Combination": ["NDM + OXA-48","NDM + KPC","mecA + cfr",
                            "mecA + vanA","MCR + NDM","KPC + VIM"],
            "Clinical significance": [
                "Pan-resistant — extremely limited treatment options",
                "Dual carbapenemase — carbapenem resistance confirmed",
                "MRSA + linezolid resistance — very few options remaining",
                "MRSA + VRE genes — multi-drug resistant",
                "Colistin + carbapenem resistance — pan-resistant potential",
                "Dual carbapenemase — pan-resistant potential"
            ]
        }))
        st.markdown("""
        **Carrier organism identification:**
        The pipeline uses contig-level Kraken2 classification to identify which organism 
        is most likely carrying each resistance gene. Results are reported as:
        - **Confirmed** — organism detected in sample AND contig resolves to species level
        - **Inferred** — gene present but contig only resolves to genus or class level
        - **Unresolved** — contig too short or ambiguous for species-level classification

        Note: Short-read metagenomics cannot definitively prove which organism carries 
        a resistance gene. Long-read sequencing (Nanopore) or culture-based susceptibility 
        testing is required for definitive linkage.
        """)

    with st.expander("9️⃣ Interpreting Results"):
        st.markdown("""
        **Abundance percentage:**
        Relative abundance is calculated as organism reads / total microbial reads × 100%.
        Bacterial and fungal fractions are calculated separately. An organism at 95% 
        bacterial abundance means 95% of all classified bacterial reads mapped to that organism.

        **Detection threshold:**
        Only organisms at ≥3% relative abundance with ≥50 reads are reported. 
        Organisms below this threshold may be present but cannot be reliably distinguished 
        from background noise.

        **Infection types:**
        - **Monomicrobial** — one organism clearly dominant (>80% abundance)
        - **Polymicrobial** — 2-3 organisms above threshold
        - **Complex polymicrobial** — 4+ organisms above threshold
        - **No detection** — no organism above 3% threshold

        **Co-infection detection:**
        When multiple organisms are detected above threshold, MAPS reports them as 
        a co-infection. The clinical summary reflects the dominant organism.

        **Confidence intervals:**
        Wilson score 95% confidence intervals are calculated for each organism's 
        abundance estimate. Wide intervals (Low precision) indicate fewer supporting reads.

        **Anomaly detection:**
        Isolation Forest algorithm compares each sample's profile against the 26 training 
        samples. A MODERATE or SEVERE anomaly means the sample looks statistically unusual 
        and results should be interpreted with additional caution.
        """)

    with st.expander("🔟 Limitations"):
        st.markdown("""
        **DNA vs live organisms:**
        MAPS detects DNA — not live organisms. Patients who received antibiotics before 
        sampling may have culture-negative results but MAPS-positive results because 
        dead bacterial DNA remains detectable. A MAPS-positive result does not confirm 
        active infection in isolation.

        **Read depth requirements:**
        Plasma and blood samples naturally contain 95-99% human DNA. After human depletion, 
        very few microbial reads remain. These sample types require either human DNA 
        depletion before sequencing or higher-output sequencers (MiSeq or NextSeq) 
        to generate sufficient microbial reads.

        **Sequencing platform differences:**
        Ion Torrent data has lower Q30 rates and higher indel error rates than Illumina. 
        At Q30 rates below 5% (as seen in some Ion Torrent runs), misclassification 
        increases significantly and results should be interpreted with caution.

        **Short-read assembly limitations:**
        SPAdes assembly from short reads produces fragmented contigs for complex 
        polymicrobial samples. MLST and AMR carrier identification are less reliable 
        in highly polymicrobial samples.

        **Training data size:**
        The AI models were trained on 26 confirmed clinical samples. This is sufficient 
        for initial validation but models will improve as more confirmed samples are added.

        **Not a replacement for culture:**
        MAPS is a complementary tool. Culture remains the gold standard for organism 
        viability, susceptibility testing, and clinical correlation. MAPS results should 
        always be interpreted alongside clinical presentation and other diagnostic findings.

        **No FDA clearance:**
        MAPS is currently a laboratory-developed tool under validation. Results should 
        not be used as the sole basis for clinical decision-making.
        """)

    with st.expander("1️⃣1️⃣ Glossary"):
        st.markdown("""
        | Term | Definition |
        |---|---|
        | **AMR** | Antimicrobial resistance — the ability of a microorganism to resist the effects of drugs |
        | **Abundance** | The relative proportion of an organism's DNA in a sample |
        | **Assembly** | Reconstructing a genome sequence by piecing together short sequencing reads |
        | **BAL** | Bronchoalveolar lavage — fluid collected from the lungs during bronchoscopy |
        | **CARD** | Comprehensive Antibiotic Resistance Database — a curated database of AMR genes |
        | **Contig** | A contiguous sequence assembled from overlapping reads |
        | **CRE** | Carbapenem-resistant Enterobacteriaceae — resistant to last-resort antibiotics |
        | **ESBL** | Extended-spectrum beta-lactamase — enzyme that breaks down many penicillins and cephalosporins |
        | **FASTQ** | Raw sequencing file format containing reads and quality scores |
        | **Human depletion** | Computationally removing human reads to enrich for microbial signal |
        | **K-mer** | A short subsequence of length k used by Kraken2 for rapid classification |
        | **Kraken2** | K-mer based taxonomic classification tool |
        | **MLST** | Multilocus sequence typing — strain-level identification using housekeeping genes |
        | **MRSA** | Methicillin-resistant Staphylococcus aureus |
        | **Q20** | Base quality score meaning 99% chance the base is called correctly |
        | **Q30** | Base quality score meaning 99.9% chance the base is called correctly |
        | **Relative abundance** | Organism reads / total microbial reads × 100% |
        | **ResFinder** | Database of acquired resistance genes |
        | **SHAP** | SHapley Additive exPlanations — method for explaining AI model predictions |
        | **SPAdes** | De novo genome assembler used to reconstruct organism genomes from reads |
        | **VRE** | Vancomycin-resistant Enterococcus |
        | **VFDB** | Virulence Factor Database — curated database of bacterial virulence genes |
        """)

    st.divider()
    st.caption("MAPS v1.0 | Baylor Scott & White Health | Samrudh Rao, Scientist I")

