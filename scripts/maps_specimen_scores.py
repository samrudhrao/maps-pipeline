#!/usr/bin/env python3
"""
MAPS Specimen Scores - maps_specimen_scores.py
Provides specimen-type specific pathogen score adjustments and exclusions.
Imported by maps_final_report.py, maps_predict.py, maps_ensemble.py, maps_stats.py.

Two main components:

1. SPECIMEN_EXCLUSIONS — organisms excluded entirely for a specimen type
   e.g. Pneumocystis excluded from synovial fluid (not biologically plausible)
   To add new exclusions: add organism name (lowercase) to appropriate list

2. SPECIMEN_ADJUSTMENTS — multipliers applied to base pathogen scores
   Multiplier >1.0 = more likely pathogen in this context
   Multiplier <1.0 = less likely pathogen in this context
   To add new adjustments: add organism-multiplier pair to appropriate dict

Supported specimen types: synovial, pleural, peritoneal, plasma, bal, abscess
To add new specimen types: add entries to both SPECIMEN_EXCLUSIONS and SPECIMEN_ADJUSTMENTS
"""

# Specimen-type specific pathogen score adjustments
# Multipliers applied to base pathogen score based on sample type
# 1.0 = no change, <1.0 = less likely pathogen in this context, >1.0 = more likely

# Organisms to exclude entirely by specimen type
# These organisms are biologically implausible in the listed specimen types
# and should not appear in clinical reports regardless of abundance
SPECIMEN_EXCLUSIONS = {
    "synovial": [
        "pneumocystis",
        "clostridium",
        "anaerocolumna",
        "anaerotignum",
        "fusobacterium periodonticum",
        "prevotella",
        "veillonella",
    ],
    "pleural": [
        "anaerocolumna",
        "anaerotignum",
        "clostridium intestinale",
    ],
    "peritoneal": [],
    "plasma": [
        "anaerocolumna",
        "anaerotignum",
    ],
}

def should_exclude(organism, specimen_type):
    """
    Check if an organism should be excluded from reporting for this specimen type.
    Returns True if the organism is in the exclusion list for this specimen type.
    Called by maps_final_report.py, maps_predict.py, maps_ensemble.py.
    """
    """Return True if organism should be excluded for this specimen type."""
    if not specimen_type or specimen_type == "unknown":
        return False
    org_lower = organism.lower()
    for stype, exclusions in SPECIMEN_EXCLUSIONS.items():
        if stype in specimen_type.lower():
            if any(ex in org_lower for ex in exclusions):
                return True
    return False

# Specimen-type score multipliers
# Applied to base pathogen score from PATHOGEN_SCORES in maps_final_report.py
# Scores are capped at 1.0 after adjustment
# Example: S. aureus in synovial fluid gets 1.2x multiplier (more suspicious)
# Example: S. mitis in synovial fluid gets 0.1x multiplier (likely contamination)
SPECIMEN_ADJUSTMENTS = {
    "synovial": {
        # Oral/skin flora unlikely in synovial fluid
        "streptococcus mitis":        0.1,
        "streptococcus sanguinis":    0.1,
        "streptococcus oralis":       0.1,
        "gemella haemolysans":        0.05,
        "simonsiella muelleri":       0.05,
        "rothia mucilaginosa":        0.1,
        "veillonella parvula":        0.1,
        "prevotella melaninogenica":  0.1,
        "haemophilus parainfluenzae": 0.2,
        "fusobacterium":              0.2,
        # Gut anaerobes — not expected in joint fluid
        "clostridium":                0.05,
        "anaerocolumna":              0.05,
        "anaerotignum":               0.05,
        # Lung pathogens — not expected in joint fluid
        "pneumocystis":               0.05,
        # True joint pathogens
        "staphylococcus aureus":      1.2,
        "streptococcus pyogenes":     1.2,
        "neisseria gonorrhoeae":      1.2,
        "pseudomonas aeruginosa":     1.1,
    },
    "pleural": {
        # Oral flora less likely in pleural fluid
        "streptococcus mitis":        0.2,
        "gemella haemolysans":        0.1,
        "simonsiella muelleri":       0.05,
        # True pleural pathogens
        "streptococcus pneumoniae":   1.2,
        "staphylococcus aureus":      1.1,
        "pseudomonas aeruginosa":     1.1,
        "mycobacterium tuberculosis": 1.3,
        "fusobacterium nucleatum":    1.1,
    },
    "bal": {
        # Oral flora more acceptable in BAL due to contamination during bronchoscopy
        "streptococcus mitis":        0.5,
        "haemophilus parainfluenzae": 0.6,
        # True respiratory pathogens boosted
        "pseudomonas aeruginosa":     1.2,
        "klebsiella pneumoniae":      1.2,
        "aspergillus fumigatus":      1.3,
        "pneumocystis jirovecii":     1.3,
        "staphylococcus aureus":      1.1,
    },
    "plasma": {
        # Blood pathogens
        "candida albicans":           1.2,
        "staphylococcus aureus":      1.2,
        "escherichia coli":           1.1,
        "klebsiella pneumoniae":      1.1,
        "fusobacterium nucleatum":    1.1,
    },
    "peritoneal": {
        # Gut flora more expected but still significant
        "escherichia coli":           1.1,
        "klebsiella pneumoniae":      1.1,
        "candida albicans":           1.2,
        "pseudomonas aeruginosa":     1.1,
        "bacteroides fragilis":       1.1,
    },
    "abscess": {
        "staphylococcus aureus":      1.2,
        "streptococcus pyogenes":     1.2,
        "escherichia coli":           1.1,
        "candida albicans":           1.1,
    },
}

def get_specimen_adjusted_score(base_score, organism, specimen_type):
    """
    Adjust pathogen score based on clinical specimen type context.
    base_score: original pathogen score from PATHOGEN_SCORES (0.0-1.0)
    organism: organism name string (with or without [BACT]/[FUNG] prefix)
    specimen_type: clinical specimen type string
    Returns adjusted score capped at 1.0
    """
    """Adjust pathogen score based on specimen type context."""
    if not specimen_type:
        return base_score
    
    specimen_type = specimen_type.lower().strip()
    organism_lower = organism.lower()
    
    # Find matching specimen type
    adjustments = None
    for stype in SPECIMEN_ADJUSTMENTS:
        if stype in specimen_type:
            adjustments = SPECIMEN_ADJUSTMENTS[stype]
            break
    
    if not adjustments:
        return base_score
    
    # Find matching organism
    multiplier = 1.0
    for org_key, mult in adjustments.items():
        if org_key in organism_lower:
            multiplier = mult
            break
    
    return min(1.0, base_score * multiplier)

if __name__ == "__main__":
    # Quick test to verify specimen adjustments are working correctly
    # Run: python3 maps_specimen_scores.py
    # Expected: S. mitis score reduced in synovial, S. aureus score increased
    tests = [
        ("Streptococcus mitis", "synovial", 0.7),
        ("Staphylococcus aureus", "synovial", 0.7),
        ("Aspergillus fumigatus", "bal", 0.8),
        ("Fusobacterium nucleatum", "plasma", 0.6),
    ]
    print("Specimen-type score adjustments:")
    for org, stype, base in tests:
        adjusted = get_specimen_adjusted_score(base, org, stype)
        print(f"  {org:<35} {stype:<12} base:{base:.2f} → adjusted:{adjusted:.2f}")
