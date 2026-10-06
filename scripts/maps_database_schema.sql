-- MAPS Results Database Schema
-- Run this to initialize a fresh database
-- Usage: sqlite3 maps_results.db < maps_database_schema.sql

CREATE TABLE IF NOT EXISTS results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample TEXT NOT NULL,
    run_date TEXT,
    quality_flag TEXT,
    human_pct REAL,
    microbial_reads INTEGER,
    infection_type TEXT,
    top_organism TEXT,
    top_ai_score REAL,
    organisms_above_3pct INTEGER,
    soc_result TEXT DEFAULT 'Unknown',
    soc_match TEXT DEFAULT 'Unknown',
    clinical_summary TEXT,
    timestamp TEXT,
    low_abundance_detections TEXT
);
