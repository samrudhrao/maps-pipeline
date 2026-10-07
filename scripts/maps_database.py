import sqlite3
import json
import os
import sys
from datetime import datetime

DB_PATH = os.path.expanduser("~/maps/maps_results.db")

def init_db():
    conn = sqlite3.connect(DB_PATH)
    c    = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS results (
        id              INTEGER PRIMARY KEY AUTOINCREMENT,
        sample          TEXT NOT NULL,
        run_date        TEXT,
        quality_flag    TEXT,
        human_pct       REAL,
        microbial_reads INTEGER,
        infection_type  TEXT,
        top_organism    TEXT,
        top_ai_score    REAL,
        organisms_above_3pct INTEGER,
        soc_result      TEXT,
        soc_match       TEXT,
        clinical_summary TEXT,
        timestamp       TEXT
    )''')
    conn.commit()
    conn.close()

def save_result(sample, outdir, soc_result=None, soc_match=None):
    init_db()
    json_path = f"{outdir}/{sample}_summary.json"
    if not os.path.exists(json_path):
        print(f"No summary JSON found for {sample}")
        return
    with open(json_path) as f:
        summary = json.load(f)
    conn = sqlite3.connect(DB_PATH)
    c    = conn.cursor()

    # Preserve existing SoC data when rerunning a sample
    c.execute("SELECT soc_result, soc_match FROM results WHERE sample=?", (sample,))
    existing = c.fetchone()
    preserved_soc    = existing[0] if existing and existing[0] not in (None, "Unknown") else soc_result or "Unknown"
    preserved_match  = existing[1] if existing and existing[1] not in (None, "Unknown") else soc_match or "Unknown"

    c.execute("DELETE FROM results WHERE sample=?", (sample,))
    c.execute('''INSERT INTO results
        (sample, run_date, quality_flag, human_pct, microbial_reads,
         infection_type, top_organism, top_ai_score, organisms_above_3pct,
         soc_result, soc_match, clinical_summary, timestamp)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', (
        sample,
        datetime.now().strftime("%Y-%m-%d"),
        summary.get("quality","unknown"),
        summary.get("human_pct",0),
        summary.get("microbial_reads",0),
        summary.get("infection_type","unknown"),
        summary.get("top_organism","None"),
        summary.get("top_ai_score",0),
        summary.get("organisms_above_3pct",0),
        preserved_soc,
        preserved_match,
        summary.get("clinical_summary",""),
        datetime.now().strftime("%Y-%m-%d %H:%M"),
    ))
    conn.commit()
    conn.close()
    print(f"Saved {sample} to database")

def query_all():
    init_db()
    conn = sqlite3.connect(DB_PATH)
    c    = conn.cursor()
    c.execute('''SELECT sample, run_date, quality_flag, top_organism,
                        top_ai_score, soc_match
                 FROM results ORDER BY timestamp DESC''')
    rows = c.fetchall()
    conn.close()
    print(f"\n{'='*80}")
    print(f"MAPS RESULTS DATABASE — {len(rows)} samples")
    print(f"{'='*80}")
    print(f"{'Sample':<14} {'Date':<12} {'Quality':<10} {'Top Organism':<35} {'AI':>7} {'SoC Match':>10}")
    print("─"*80)
    for row in rows:
        sample, date, qual, top_org, top_ai, match = row
        org_short = (top_org or "None")[:33]
        print(f"{sample:<14} {date:<12} {qual:<10} {org_short:<35} "
              f"{(top_ai or 0)*100:>6.1f}% {match or 'Unknown':>10}")

if __name__ == "__main__":
    init_db()
    if len(sys.argv) >= 3 and sys.argv[1] == "save":
        save_result(sys.argv[2], sys.argv[3],
                   sys.argv[4] if len(sys.argv) > 4 else None,
                   sys.argv[5] if len(sys.argv) > 5 else None)
    else:
        query_all()
