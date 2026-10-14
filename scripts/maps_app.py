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
tab1, = st.tabs(["🏠 Run Pipeline"])

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
            "bash", os.path.expanduser("/maps/full_pipeline.sh"),
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

            # QC Results
            qc_path = f"{outdir}/{sample_name}_qc.json"
            if os.path.exists(qc_path):
                with open(qc_path) as f:
                    qc = json.load(f)
                st.markdown("<h3 style='color:#003087;'>🔬 QC Report</h3>",
                            unsafe_allow_html=True)
                qc_status = qc.get("status","UNKNOWN")
                if qc_status == "PASS":
                    st.success(f"✅ QC Status: PASS")
                else:
                    st.error(f"❌ QC Status: {qc_status}")
                q1, q2, q3, q4, q5, q6 = st.columns(6)
                q1.metric("Total Reads",     f"{qc.get('total_reads',0):,}")
                q2.metric("Human DNA",       f"{qc.get('human_pct',0):.1f}%")
                q3.metric("Microbial Reads", f"{qc.get('microbial_reads',0):,}")
                q4.metric("Q30 Rate",        f"{qc.get('q30_rate',0):.1f}%")
                q5.metric("Duplication",     f"{qc.get('duplication_pct',0):.1f}%")
                q6.metric("Insert Size",     f"{qc.get('insert_size',0)}bp")
                checks = qc.get("checks", {})
                if not checks.get("total_reads_ok", True):
                    st.warning("⚠️ Total reads below minimum threshold (50,000)")
                if not checks.get("human_pct_ok", True):
                    st.warning("⚠️ Human DNA above maximum threshold (90%)")
                if not checks.get("microbial_reads_ok", True):
                    st.warning("⚠️ Microbial reads below minimum threshold (1,000)")
                st.divider()

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
