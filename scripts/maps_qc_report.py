#!/usr/bin/env python3
"""
MAPS QC Failure Report - maps_qc_report.py
Creates a PDF explaining why a sample failed QC, from OUTDIR/SAMPLE_qc.json.

Usage:  python3 maps_qc_report.py SAMPLE OUTDIR [SPECIMEN_TYPE]
Output: OUTDIR/SAMPLE_qc_failure_report.pdf
Operator initials and run start come from MAPS_USER / MAPS_RUN_START (set by the app).
"""

import sys
import os
import json
from datetime import datetime
from xml.sax.saxutils import escape

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, HRFlowable)

BSW_BLUE = colors.HexColor("#003087")
FAIL_RED = colors.HexColor("#C00000")
FAIL_BG  = colors.HexColor("#FDE9E7")
WARN_BG  = colors.HexColor("#FFF4CE")
LABEL_BG = colors.HexColor("#DCE6F1")

PLATFORM_NAMES = {"iontorrent": "Ion Torrent (single-end)",
                  "illumina":   "Illumina (paired-end)"}

NEXT_STEPS = {
    "Total reads":     "Total reads: confirm the correct FASTQ was uploaded and the library was loaded "
                       "adequately; re-sequence with greater depth if needed.",
    "Human DNA":       "Human DNA: high host background. Consider host depletion during extraction "
                       "or re-extraction before re-sequencing.",
    "Microbial reads": "Microbial reads: very little non-human DNA. This may be a low-biomass specimen, "
                       "a true negative, or an extraction failure; check the run controls.",
    "Q20 rate":        "Q20 rate: low base quality. Review the Ion Torrent run metrics (chip loading, "
                       "polyclonality) and consider re-sequencing.",
    "Q30 rate":        "Q30 rate: low base quality. Review the Illumina run metrics (cluster density, "
                       "%PF) and consider re-sequencing.",
}


def metric_rows(qc):
    """Every QC parameter as [parameter, result, required, status, type]."""
    checks   = qc.get("checks", {})
    platform = qc.get("platform", "illumina")

    def hard(ok):
        return "PASS" if ok else "FAIL"

    rows = [
        ["Total reads", f"{qc.get('total_reads', 0):,}", ">= 50,000",
         hard(checks.get("total_reads_ok", True)), "Hard fail"],
        ["Human DNA", f"{qc.get('human_pct', 0):.1f}%", "<= 90%",
         hard(checks.get("human_pct_ok", True)), "Hard fail"],
        ["Microbial reads", f"{qc.get('microbial_reads', 0):,}", ">= 1,000",
         hard(checks.get("microbial_reads_ok", True)), "Hard fail"],
    ]
    if platform == "iontorrent":
        q_name, q_val, q_min = "Q20 rate", qc.get("q20_rate", 0), 70
    else:
        q_name, q_val, q_min = "Q30 rate", qc.get("q30_rate", 0), 50
    if q_val:
        q_ok = checks.get("quality_ok", q_val >= q_min)
        rows.append([q_name, f"{q_val:.1f}%", f">= {q_min}%", hard(q_ok), "Hard fail"])
    else:
        rows.append([q_name, "Not available", f">= {q_min}%", "N/A", "Hard fail"])

    dup = qc.get("duplication_pct", 0)
    rows.append(["Duplication", f"{dup:.1f}%" if dup else "Not available", "<= 80%",
                 ("WARN" if dup > 80 else "PASS") if dup else "N/A", "Soft warning"])
    ins = qc.get("insert_size", 0)
    rows.append(["Insert size", f"{ins} bp" if ins else "Not available", ">= 100 bp",
                 ("WARN" if ins < 100 else "PASS") if ins else "N/A", "Soft warning"])
    return rows


def build_report(sample, outdir, specimen_type="unknown", user=None,
                 run_start=None, pdf_path=None):
    with open(os.path.join(outdir, f"{sample}_qc.json")) as f:
        qc = json.load(f)

    pdf_path  = pdf_path or os.path.join(outdir, f"{sample}_qc_failure_report.pdf")
    user      = user or os.environ.get("MAPS_USER") or "N/A"
    run_start = run_start or os.environ.get("MAPS_RUN_START") or "N/A"
    now       = datetime.now().strftime("%Y-%m-%d %H:%M")

    styles = getSampleStyleSheet()
    title  = ParagraphStyle("T", parent=styles["Title"], textColor=BSW_BLUE,
                            fontSize=18, spaceAfter=6)
    h2     = ParagraphStyle("H2", parent=styles["Heading2"], textColor=BSW_BLUE,
                            fontSize=13, spaceBefore=10, spaceAfter=4)
    normal = styles["Normal"]
    cell   = ParagraphStyle("Cell", parent=normal, fontSize=9, leading=11)
    cell_b = ParagraphStyle("CellB", parent=cell, fontName="Helvetica-Bold")
    head   = ParagraphStyle("Head", parent=cell_b, textColor=colors.white)
    banner = ParagraphStyle("Banner", parent=normal, fontSize=11, leading=14,
                            fontName="Helvetica-Bold", textColor=colors.white)
    small  = ParagraphStyle("Small", parent=normal, fontSize=8, textColor=colors.grey)

    def P(text, style=cell):
        return Paragraph(escape(str(text)), style)

    base_style = [
        ("GRID",    (0, 0), (-1, -1), 0.5, colors.grey),
        ("VALIGN",  (0, 0), (-1, -1), "MIDDLE"),
        ("PADDING", (0, 0), (-1, -1), 5),
    ]

    rows     = metric_rows(qc)
    failures = qc.get("failures")
    if failures is None:  # QC JSON written before failures were recorded
        failures = [{"parameter": r[0], "value": r[1], "threshold": r[2], "reason": ""}
                    for r in rows if r[3] == "FAIL"]
    warnings = qc.get("warnings", [])

    story = [
        Paragraph("MAPS QC Failure Report", title),
        Paragraph("Metagenomic Abundance and Pathogen Scoring", normal),
        Paragraph("Baylor Scott &amp; White Health", normal),
        HRFlowable(width="100%", thickness=2, color=BSW_BLUE),
        Spacer(1, 0.12 * inch),
    ]

    n = len(failures)
    msg = (f"QC FAILED — {n} parameter{'s' if n != 1 else ''} failed. The pipeline stopped "
           f"at QC validation (Step 3b); no organism results were generated.")
    t = Table([[P(msg, banner)]], colWidths=[7 * inch])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), FAIL_RED),
                           ("PADDING", (0, 0), (-1, -1), 8)]))
    story.append(t)

    # Run information
    story.append(Paragraph("Run Information", h2))
    platform = qc.get("platform", "")
    info = [
        ["Sample",            sample],
        ["Specimen type",     specimen_type],
        ["Platform",          PLATFORM_NAMES.get(platform, platform or "Not recorded")],
        ["Run started",       run_start],
        ["Report generated",  now],
        ["Operator initials", user],
    ]
    t = Table([[P(a, cell_b), P(b)] for a, b in info], colWidths=[2 * inch, 5 * inch])
    t.setStyle(TableStyle(base_style + [("BACKGROUND", (0, 0), (0, -1), LABEL_BG)]))
    story.append(t)

    # Failed parameters
    story.append(Paragraph("Failed Parameters", h2))
    if failures:
        data = [[P(h, head) for h in ["Parameter", "Result", "Required", "Reason"]]]
        for fl in failures:
            data.append([P(fl["parameter"], cell_b), P(fl["value"]),
                         P(fl["threshold"]), P(fl.get("reason", ""))])
        t = Table(data, colWidths=[1.6 * inch, 1.2 * inch, 1.2 * inch, 3 * inch], repeatRows=1)
        t.setStyle(TableStyle(base_style + [
            ("BACKGROUND", (0, 0), (-1, 0), FAIL_RED),
            ("BACKGROUND", (0, 1), (-1, -1), FAIL_BG),
        ]))
        story.append(t)
    else:
        story.append(Paragraph("No failed parameter was recorded in the QC file. "
                               "Check the pipeline log for the stopping reason.", normal))

    # All QC metrics
    story.append(Paragraph("All QC Metrics", h2))
    data = [[P(h, head) for h in ["Parameter", "Result", "Required", "Status", "Type"]]]
    data += [[P(c) for c in r] for r in rows]
    ts = base_style + [("BACKGROUND", (0, 0), (-1, 0), BSW_BLUE)]
    for i, r in enumerate(rows, start=1):
        if r[3] == "FAIL":
            ts.append(("BACKGROUND", (0, i), (-1, i), FAIL_BG))
        elif r[3] == "WARN":
            ts.append(("BACKGROUND", (0, i), (-1, i), WARN_BG))
    t = Table(data, colWidths=[1.6 * inch, 1.3 * inch, 1.2 * inch, 0.9 * inch, 2 * inch],
              repeatRows=1)
    t.setStyle(TableStyle(ts))
    story.append(t)

    # Soft warnings
    if warnings:
        story.append(Paragraph("Soft Warnings", h2))
        for w in warnings:
            story.append(Paragraph(
                f"• <b>{escape(w['parameter'])}</b>: {escape(w['value'])} — "
                f"{escape(w.get('reason', ''))}", normal))

    # Suggested next steps
    steps = [NEXT_STEPS[fl["parameter"]] for fl in failures if fl["parameter"] in NEXT_STEPS]
    if steps:
        story.append(Paragraph("Suggested Next Steps", h2))
        for s in steps:
            story.append(Paragraph(f"• {escape(s)}", normal))
            story.append(Spacer(1, 0.04 * inch))

    story += [
        Spacer(1, 0.25 * inch),
        HRFlowable(width="100%", thickness=0.5, color=colors.grey),
        Paragraph(escape(f"Generated by MAPS QC module | {now} | Source: {sample}_qc.json"), small),
    ]

    doc = SimpleDocTemplate(pdf_path, pagesize=letter,
                            rightMargin=0.75 * inch, leftMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch)
    doc.build(story)
    return pdf_path


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 maps_qc_report.py SAMPLE OUTDIR [SPECIMEN_TYPE]")
        sys.exit(1)
    path = build_report(sys.argv[1], sys.argv[2],
                        sys.argv[3] if len(sys.argv) > 3 else "unknown")
    print(f"  QC failure report written: {path}")
