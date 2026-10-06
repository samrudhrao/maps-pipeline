#!/bin/bash
SAMPLE=$1
R1=$2
R2=$3
OUTDIR=$4
SPECIMEN_TYPE=${5:-unknown}
THRESHOLD=3.0

mkdir -p $OUTDIR/{fastqc,fastp,human_depleted,bact_results,fungi_results,assembly,mlst,virulence,multiqc}

# Convert BAM to FASTQ if input is BAM
if [[ "$R1" == *.bam ]]; then
    echo "BAM file detected — converting to FASTQ..."
    BAM_BASE=$OUTDIR/fastp/${SAMPLE}_from_bam
    samtools sort -n $R1 -o ${BAM_BASE}_sorted.bam 2>/dev/null

    # Detect PacBio by checking header
    IS_PACBIO=$(samtools view -H $R1 2>/dev/null | grep -c "PL:PACBIO" || echo 0)

    if [ "$IS_PACBIO" -gt 0 ]; then
        echo "  PacBio HiFi BAM detected — using single-end conversion"
        samtools fastq -0 ${BAM_BASE}_R1.fastq -n ${BAM_BASE}_sorted.bam 2>/dev/null
        R1=${BAM_BASE}_R1.fastq
        R2=${BAM_BASE}_R1.fastq
        echo "  PacBio single-end BAM converted successfully"
    else
        samtools fastq -1 ${BAM_BASE}_R1.fastq.gz -2 ${BAM_BASE}_R2.fastq.gz -0 /dev/null -s /dev/null -n ${BAM_BASE}_sorted.bam 2>/dev/null
        if [ -s ${BAM_BASE}_R2.fastq.gz ]; then
            R1=${BAM_BASE}_R1.fastq.gz
            R2=${BAM_BASE}_R2.fastq.gz
            echo "  Paired-end BAM converted successfully"
        else
            samtools fastq -0 ${BAM_BASE}_R1.fastq.gz -n $R1 2>/dev/null
            R1=${BAM_BASE}_R1.fastq.gz
            R2=${BAM_BASE}_R1.fastq.gz
            echo "  Single-end BAM converted successfully"
        fi
    fi
    echo "  BAM conversion done"
fi

# Detect sequencing platform from input
# If R2 is empty or same as R1 = single-end (Ion Torrent)
# If R2 is different = paired-end (Illumina)
if [ -z "$R2" ] || [ "$R1" = "$R2" ]; then
    PLATFORM="iontorrent"
    echo "Platform detected: Ion Torrent (single-end)"
else
    PLATFORM="illumina"
    echo "Platform detected: Illumina (paired-end)"
fi

echo "================================================"
echo "MAPS Pipeline — Sample: $SAMPLE"
echo "================================================"

# Step 1: FastQC
echo "Step 1: FastQC quality check..."
fastqc $R1 $R2 -o $OUTDIR/fastqc/ --threads 8 -q
echo "  FastQC done"

# Step 2: fastp
echo "Step 2: fastp trimming..."
if [ "$PLATFORM" = "iontorrent" ]; then
    fastp \
        -i $R1 \
        -o $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        --html $OUTDIR/fastp/${SAMPLE}_fastp.html \
        --json $OUTDIR/fastp/${SAMPLE}_fastp.json \
        --thread 8 \
        --qualified_quality_phred 20 \
        --length_required 50
    # Create empty R2 placeholder for downstream compatibility
    touch $OUTDIR/fastp/${SAMPLE}_R2_clean.fastq.gz
else
    fastp \
        -i $R1 -I $R2 \
        -o $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        -O $OUTDIR/fastp/${SAMPLE}_R2_clean.fastq.gz \
        --html $OUTDIR/fastp/${SAMPLE}_fastp.html \
        --json $OUTDIR/fastp/${SAMPLE}_fastp.json \
        --thread 8 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --detect_adapter_for_pe
fi
echo "  fastp done"

# Step 3: Human depletion
echo "Step 3: Human read depletion..."
if [ "$PLATFORM" = "iontorrent" ]; then
    kraken2 \
        --db ~/maps/kraken2_db \
        --gzip-compressed \
        --report $OUTDIR/human_depleted/${SAMPLE}_human_report.txt \
        --output /dev/null \
        --unclassified-out $OUTDIR/human_depleted/${SAMPLE}_clean.fastq \
        --threads 8 \
        $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz
    cp $OUTDIR/human_depleted/${SAMPLE}_clean.fastq \
       $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq
    cp $OUTDIR/human_depleted/${SAMPLE}_clean.fastq \
       $OUTDIR/human_depleted/${SAMPLE}_clean_2.fastq
else
    kraken2 \
        --db ~/maps/kraken2_db \
        --paired --gzip-compressed \
        --report $OUTDIR/human_depleted/${SAMPLE}_human_report.txt \
        --output /dev/null \
        --unclassified-out $OUTDIR/human_depleted/${SAMPLE}_clean#.fastq \
        --threads 8 \
        $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        $OUTDIR/fastp/${SAMPLE}_R2_clean.fastq.gz
fi
echo "  Human depletion done"

# Step 3b: QC check
echo "Step 3b: QC check..."
python3 ~/maps/scripts/maps_qc.py $SAMPLE $OUTDIR
QC_STATUS=$(python3 -c "import json; d=json.load(open('$OUTDIR/${SAMPLE}_qc.json')); print(d['status'])" 2>/dev/null)
if [ "$QC_STATUS" = "FAIL" ]; then
    echo "  ✗ QC FAILED — stopping pipeline"
    echo "  Check $OUTDIR/${SAMPLE}_qc.json for details"
    exit 1
fi
echo "  QC done"

# Step 4: Bacterial classification
echo "Step 4: Bacterial classification..."
if [ "$PLATFORM" = "iontorrent" ]; then
    kraken2 \
        --db ~/maps/kraken2_db \
        --gzip-compressed \
        --report $OUTDIR/bact_results/${SAMPLE}_bact_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz
else
    kraken2 \
        --db ~/maps/kraken2_db \
        --paired --gzip-compressed \
        --report $OUTDIR/bact_results/${SAMPLE}_bact_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        $OUTDIR/fastp/${SAMPLE}_R2_clean.fastq.gz
fi
echo "  Bacterial classification done"

# Step 4b: Supplemental DB (only if Staphylococcus > 5%)
echo "Step 4b: Supplemental database classification..."
STAPH_READS=$(grep -w "Staphylococcus$" $OUTDIR/bact_results/${SAMPLE}_bact_report.txt 2>/dev/null | awk '{print $2}' | head -1)
STAPH_READS=${STAPH_READS:-0}
STAPH_CHECK=$(awk "BEGIN {print ($STAPH_READS > 500) ? 1 : 0}")
if [ "$STAPH_CHECK" = "1" ] && [ -d ~/maps/kraken2_db_extra ]; then
    if [ "$PLATFORM" = "iontorrent" ]; then
        kraken2 \
            --db ~/maps/kraken2_db_extra \
            --report $OUTDIR/bact_results/${SAMPLE}_extra_report.txt \
            --output /dev/null \
            --threads 8 \
            $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq
    else
        kraken2 \
            --db ~/maps/kraken2_db_extra \
            --paired --confidence 0.1 \
            --report $OUTDIR/bact_results/${SAMPLE}_extra_report.txt \
            --output /dev/null \
            --threads 8 \
            $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq \
            $OUTDIR/human_depleted/${SAMPLE}_clean_2.fastq
    fi
    echo "  Supplemental done (Staphylococcus reads: ${STAPH_READS} detected)"
else
    touch $OUTDIR/bact_results/${SAMPLE}_extra_report.txt
    echo "  Staphylococcus not dominant — skipping supplemental DB"
fi

# Step 5: Fungal classification
echo "Step 5: Fungal classification..."
if [ "$PLATFORM" = "iontorrent" ]; then
    kraken2 \
        --db ~/maps/kraken2_db_fungal_refseq \
        --report $OUTDIR/fungi_results/${SAMPLE}_fungi_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq
else
    kraken2 \
        --db ~/maps/kraken2_db_fungal_refseq \
        --paired \
        --report $OUTDIR/fungi_results/${SAMPLE}_fungi_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq \
        $OUTDIR/human_depleted/${SAMPLE}_clean_2.fastq
fi
echo "  Fungal classification done"

# Step 5b: Viral classification
echo "Step 5b: Viral classification..."
if [ "$PLATFORM" = "iontorrent" ]; then
    kraken2 \
        --db ~/maps/kraken2_db_viral \
        --report $OUTDIR/bact_results/${SAMPLE}_viral_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq
else
    kraken2 \
        --db ~/maps/kraken2_db_viral \
        --paired \
        --report $OUTDIR/bact_results/${SAMPLE}_viral_report.txt \
        --output /dev/null \
        --threads 8 \
        $OUTDIR/human_depleted/${SAMPLE}_clean_1.fastq \
        $OUTDIR/human_depleted/${SAMPLE}_clean_2.fastq
fi
echo "  Viral classification done"
echo "  Viruses detected:"
awk '$4=="S" && $1>0.1 {print "   ", $1"%", $6, $7, $8}' \
    $OUTDIR/bact_results/${SAMPLE}_viral_report.txt | head -5

# Step 5c: Viral clinical report
echo "Step 5c: Viral screening report..."
python3 ~/maps/scripts/maps_viral_report.py $SAMPLE $OUTDIR 2>/dev/null
echo "  Viral screening done"

# Step 6: Assembly
echo "Step 6: Genome assembly (SPAdes)..."
if [ "$PLATFORM" = "iontorrent" ]; then
    spades.py \
        -s $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        -o $OUTDIR/assembly/${SAMPLE} \
        --threads 8 \
        --iontorrent \
        --careful
else
    spades.py \
        -1 $OUTDIR/fastp/${SAMPLE}_R1_clean.fastq.gz \
        -2 $OUTDIR/fastp/${SAMPLE}_R2_clean.fastq.gz \
        -o $OUTDIR/assembly/${SAMPLE} \
        --threads 8 \
        --careful
fi
echo "  Assembly done"

# Step 7: MLST
echo "Step 7: MLST strain typing..."
mlst \
    --datadir ~/maps/mlst-2.23.0/db/pubmlst \
    --blastdb ~/maps/mlst-2.23.0/db/blast/mlst \
    $OUTDIR/assembly/${SAMPLE}/scaffolds.fasta \
    > $OUTDIR/mlst/${SAMPLE}_mlst.txt 2>/dev/null
cat $OUTDIR/mlst/${SAMPLE}_mlst.txt
echo "  MLST done"

# Step 8: Virulence genes
echo "Step 8: Virulence gene detection..."
abricate \
    --db vfdb \
    $OUTDIR/assembly/${SAMPLE}/scaffolds.fasta \
    > $OUTDIR/virulence/${SAMPLE}_virulence.txt 2>/dev/null
echo "  Virulence genes: $(grep -v '#' $OUTDIR/virulence/${SAMPLE}_virulence.txt | wc -l)"
echo "  ABRicate done"

# Step 8b: AMR resistance detection
echo "Step 8b: AMR resistance detection..."
abricate \
    --db card \
    $OUTDIR/assembly/${SAMPLE}/scaffolds.fasta \
    > $OUTDIR/virulence/${SAMPLE}_resistance.txt 2>/dev/null
abricate \
    --db resfinder \
    $OUTDIR/assembly/${SAMPLE}/scaffolds.fasta \
    >> $OUTDIR/virulence/${SAMPLE}_resistance.txt 2>/dev/null
python3 ~/maps/scripts/maps_amr_report.py $SAMPLE $OUTDIR
echo "  AMR detection done"

# Step 9: Abundance report
echo "Step 9: Generating abundance report..."
python3 << PYEOF
import pandas as pd
import os

with open(os.path.expanduser("~/maps/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]

THRESHOLD = 3.0

def parse_report(path, min_reads=50):
    rows = {}
    if not os.path.exists(path): return rows
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            pct, reads = float(parts[0]), int(parts[1])
            rank, name = parts[3].strip(), parts[5].strip()
            if rank == "S" and reads >= min_reads:
                if not any(e in name.lower() for e in EXCLUDE):
                    rows[name] = pct
    return rows

bact = parse_report("$OUTDIR/bact_results/${SAMPLE}_bact_report.txt")
fung = parse_report("$OUTDIR/fungi_results/${SAMPLE}_fungi_report.txt")
supp = parse_report("$OUTDIR/bact_results/${SAMPLE}_extra_report.txt")

# Merge supplemental DB — only S. cohnii and B. licheniformis
# Skip S. cohnii if main DB already has >500 reads of any other Staphylococcus
SUPP_KEEP = ["staphylococcus cohnii"]  # Only add S. cohnii from supplemental

# Read raw Staphylococcus read counts from main DB report
main_staph_reads = 0
try:
    with open("$OUTDIR/bact_results/${SAMPLE}_bact_report.txt") as _f:
        for _line in _f:
            _p = _line.strip().split("\t")
            if len(_p) < 6: continue
            if _p[3].strip() == "S" and int(_p[1]) >= 50:
                _name = _p[5].strip().lower()
                # Only count well-known pathogenic Staphylococcus species
                KNOWN_STAPH = ['aureus', 'epidermidis', 'lugdunensis',
                               'saprophyticus', 'haemolyticus', 'capitis',
                               'warneri', 'caprae', 'hominis', 'pasteuri',
                               'pettenkoferi', 'simulans', 'xylosus']
                if 'staphylococcus' in _name and any(s in _name for s in KNOWN_STAPH):
                    main_staph_reads += int(_p[1])
except: pass

for org, pct in supp.items():
    if not any(s in org.lower() for s in SUPP_KEEP):
        continue
    if 'cohnii' in org.lower() and main_staph_reads > 500:
        continue  # Main DB has dominant non-cohnii Staph — skip S. cohnii
    if org not in bact:
        bact[org] = pct

# Apply specimen-type exclusions
specimen_type = "$SPECIMEN_TYPE"
try:
    import sys as _sys
    _sys.path.insert(0, "/home/samrudh/maps")
    from maps_specimen_scores import should_exclude
    bact = {k:v for k,v in bact.items() if not should_exclude(k, specimen_type)}
    fung = {k:v for k,v in fung.items() if not should_exclude(k, specimen_type)}
except:
    pass

bact_df = pd.Series(bact)
fung_df = pd.Series(fung)
if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
bact_df.index = ["[BACT] " + x for x in bact_df.index]
fung_df.index = ["[FUNG] " + x for x in fung_df.index]
combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
above = combined[combined >= THRESHOLD]

print(f"Organisms >= {THRESHOLD}%:")
if len(above) == 0:
    print("  None above threshold")
for org, pct in above.items():
    print(f"  {org:<55} {pct:.2f}%")

combined[combined > 0].reset_index().rename(
    columns={"index":"organism", 0:"relative_abundance_%"}).to_csv(
    "$OUTDIR/${SAMPLE}_abundance.csv", index=False)
print(f"Saved: $OUTDIR/${SAMPLE}_abundance.csv")
PYEOF
echo "  Abundance report done"

echo "DONE: $SAMPLE"

# Step 10: AI scoring
echo "Step 10: AI confidence scoring..."
python3 -W ignore ~/maps/scripts/run_ai_scoring.py $SAMPLE $OUTDIR
echo "  AI scoring done"

# Step 11: Co-infection detection
echo "Step 11: Co-infection detection..."
python3 << PYEOF
import pandas as pd
import os

with open(os.path.expanduser("~/maps/filter_list.txt")) as f:
    EXCLUDE = [line.strip().lower() for line in f if line.strip()]

def parse_report(path, min_reads=50):
    rows = {}
    if not os.path.exists(path): return rows
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            pct, reads = float(parts[0]), int(parts[1])
            rank, name = parts[3].strip(), parts[5].strip()
            if rank == "S" and reads >= min_reads:
                if not any(e in name.lower() for e in EXCLUDE):
                    rows[name] = pct
    return rows

bact = parse_report("$OUTDIR/bact_results/${SAMPLE}_bact_report.txt")
fung = parse_report("$OUTDIR/fungi_results/${SAMPLE}_fungi_report.txt")
supp = parse_report("$OUTDIR/bact_results/${SAMPLE}_extra_report.txt")

# Merge supplemental DB — only S. cohnii and B. licheniformis
# Skip S. cohnii if main DB already has >500 reads of any other Staphylococcus
SUPP_KEEP = ["staphylococcus cohnii"]  # Only add S. cohnii from supplemental

# Read raw Staphylococcus read counts from main DB report
main_staph_reads = 0
try:
    with open("$OUTDIR/bact_results/${SAMPLE}_bact_report.txt") as _f:
        for _line in _f:
            _p = _line.strip().split("\t")
            if len(_p) < 6: continue
            if _p[3].strip() == "S" and int(_p[1]) >= 50:
                _name = _p[5].strip().lower()
                # Only count well-known pathogenic Staphylococcus species
                KNOWN_STAPH = ['aureus', 'epidermidis', 'lugdunensis',
                               'saprophyticus', 'haemolyticus', 'capitis',
                               'warneri', 'caprae', 'hominis', 'pasteuri',
                               'pettenkoferi', 'simulans', 'xylosus']
                if 'staphylococcus' in _name and any(s in _name for s in KNOWN_STAPH):
                    main_staph_reads += int(_p[1])
except: pass

for org, pct in supp.items():
    if not any(s in org.lower() for s in SUPP_KEEP):
        continue
    if 'cohnii' in org.lower() and main_staph_reads > 500:
        continue  # Main DB has dominant non-cohnii Staph — skip S. cohnii
    if org not in bact:
        bact[org] = pct

# Apply specimen-type exclusions
specimen_type = "$SPECIMEN_TYPE"
try:
    import sys as _sys
    _sys.path.insert(0, "/home/samrudh/maps")
    from maps_specimen_scores import should_exclude
    bact = {k:v for k,v in bact.items() if not should_exclude(k, specimen_type)}
    fung = {k:v for k,v in fung.items() if not should_exclude(k, specimen_type)}
except:
    pass
bact_df = pd.Series(bact)
fung_df = pd.Series(fung)
if bact_df.sum() > 0: bact_df = bact_df / bact_df.sum() * 100
if fung_df.sum() > 0: fung_df = fung_df / fung_df.sum() * 100
bact_df.index = ["[BACT] " + x for x in bact_df.index]
fung_df.index = ["[FUNG] " + x for x in fung_df.index]
combined = pd.concat([bact_df, fung_df]).sort_values(ascending=False)
above = combined[combined >= 3.0]

if len(above) >= 2:
    print(f"  *** CO-INFECTION — {len(above)} organisms above 3% ***")
    for org, pct in above.items():
        print(f"      {org:<50} {pct:.2f}%")
elif len(above) == 1:
    print(f"  Single pathogen: {above.index[0]} {above.iloc[0]:.2f}%")
else:
    print(f"  No organisms above 3%")

with open("$OUTDIR/${SAMPLE}_coinfection.txt", "w") as f:
    f.write(f"co_infection={len(above) >= 2}\n")
    f.write(f"organisms_above_threshold={len(above)}\n")
    for org, pct in above.items():
        f.write(f"{org},{pct:.2f}\n")
PYEOF
echo "  Co-infection check done"

# Step 12: MultiQC
echo "Step 12: MultiQC summary..."
multiqc $OUTDIR/ -o $OUTDIR/multiqc/ -q 2>/dev/null
echo "  MultiQC done"

# Step 13: Final report
echo "Step 13: Generating final clinical report..."
python3 -W ignore ~/maps/scripts/maps_final_report.py $SAMPLE $OUTDIR $SPECIMEN_TYPE
echo "  Final report done"

# Step 14: AI prediction
echo "Step 14: AI pathogen prediction..."
python3 -W ignore ~/maps/scripts/maps_predict.py $SAMPLE $OUTDIR 20.0 $SPECIMEN_TYPE
echo "  Prediction done"

# Step 14b: Ensemble prediction
echo "Step 14b: Ensemble prediction (RF + XGBoost + SHAP)..."
python3 -W ignore ~/maps/scripts/maps_ensemble.py $SAMPLE $OUTDIR 20.0 $SPECIMEN_TYPE
echo "  Ensemble done"

# Step 15: Save to database
echo "Step 15: Saving to results database..."
python3 ~/maps/scripts/maps_database.py save $SAMPLE $OUTDIR
echo "  Database updated"


# Step 16: Statistical analysis
echo "Step 16: Statistical analysis (confidence intervals + anomaly detection)..."
python3 -W ignore ~/maps/scripts/maps_stats.py $SAMPLE $OUTDIR $SPECIMEN_TYPE
echo "  Stats done"

echo ""
echo "======================================================="
echo "QC SUMMARY — $SAMPLE"
echo "======================================================="
cat $OUTDIR/${SAMPLE}_qc.json | python3 -c "
import json, sys
d = json.load(sys.stdin)
print(f'  Status:          {d.get(\"status\",\"N/A\")}')
print(f'  Total reads:     {d.get(\"total_reads\",0):,}')
print(f'  Human DNA:       {d.get(\"human_pct\",0):.1f}%')
print(f'  Microbial reads: {d.get(\"microbial_reads\",0):,}')
print(f'  Q30 rate:        {d.get(\"q30_rate\",0):.1f}%')
print(f'  Duplication:     {d.get(\"duplication_pct\",0):.1f}%')
print(f'  Insert size:     {d.get(\"insert_size\",0)}bp')
"
echo "======================================================="
echo ""
echo "================================================"
echo "MAPS COMPLETE — $SAMPLE"
echo "================================================"
