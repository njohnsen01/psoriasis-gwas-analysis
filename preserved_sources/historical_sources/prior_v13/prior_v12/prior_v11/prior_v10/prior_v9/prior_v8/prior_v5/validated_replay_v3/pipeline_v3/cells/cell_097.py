# ============================================================
# Step 4: Chunked classification of biomarker lab names
# Lipid/metabolic: LDL, HDL, TG
# Inflammatory: platelet count, ANC, ALC → SII/inflammatory index later
# ============================================================

import pandas as pd
import numpy as np

LAB_FILE = "psoriasis_biomarker_labs_raw.tsv"
OUT_FILE = "psoriasis_biomarker_labs_classified.tsv"
CHUNKSIZE = 250_000

KEEP_COLS = [
    "person_id",
    "measurement_concept_id",
    "lab_name",
    "lab_code",
    "measurement_datetime",
    "value_as_number",
    "unit_name",
    "unit_source_value",
]

# Concept IDs based on your QC output
LDL_IDS = {"3028288", "3028437", "3009966", "3053341", "3035899"}
HDL_IDS = {"3007070", "3053286"}
TG_IDS  = {"3022192", "3022038", "3027997", "3025839"}

PLT_IDS = {"3024929", "3007461", "1497147"}  # platelet count only; avoids MPV/INR/aPTT
ANC_IDS = {"3013650", "3017732", "3017501", "3015586"}
ALC_IDS = {"3004327", "3019198", "3003215"}

ABS_UNITS = [
    "thousand per microliter",
    "10*3/ul",
    "10*3/mm3",
    "k/ul",
    "10*9/l",
    "billion per liter",
    "cells per microliter",
    "cells/ul",
    "/ul",
    "/mm3",
    "x10(3)/mcl",
]

def classify_lab_chunk(df):
    for col in KEEP_COLS:
        if col not in df.columns:
            df[col] = np.nan

    df = df[KEEP_COLS].copy()

    df["measurement_concept_id"] = df["measurement_concept_id"].astype(str)
    df["lab_name_lower"] = df["lab_name"].astype(str).str.lower()
    df["lab_code_lower"] = df["lab_code"].astype(str).str.lower()
    df["unit_text"] = (
        df["unit_name"].astype(str).str.lower()
        + " | "
        + df["unit_source_value"].astype(str).str.lower()
    )

    df["biomarker"] = np.nan

    # Lipids — keep mass/volume lipid biomarkers, avoid ratios/non-HDL/VLDL
    df.loc[df["measurement_concept_id"].isin(LDL_IDS), "biomarker"] = "LDL"
    df.loc[df["measurement_concept_id"].isin(HDL_IDS), "biomarker"] = "HDL"
    df.loc[df["measurement_concept_id"].isin(TG_IDS), "biomarker"] = "TG"

    # CBC inflammatory components
    abs_unit_mask = df["unit_text"].apply(lambda x: any(u in x for u in ABS_UNITS))

    df.loc[
        df["measurement_concept_id"].isin(PLT_IDS) & abs_unit_mask,
        "biomarker"
    ] = "platelet_count"

    df.loc[
        df["measurement_concept_id"].isin(ANC_IDS) & abs_unit_mask,
        "biomarker"
    ] = "ANC"

    df.loc[
        df["measurement_concept_id"].isin(ALC_IDS) & abs_unit_mask,
        "biomarker"
    ] = "ALC"

    df = df.dropna(subset=["biomarker"]).copy()

    return df[
        [
            "person_id",
            "measurement_concept_id",
            "lab_name",
            "lab_code",
            "measurement_datetime",
            "value_as_number",
            "unit_name",
            "unit_source_value",
            "biomarker",
        ]
    ]

first_write = True
total_in = 0
total_kept = 0
kept_counts = {}

for i, chunk in enumerate(
    pd.read_csv(
        LAB_FILE,
        sep="\t",
        usecols=lambda c: c in KEEP_COLS,
        dtype=str,
        chunksize=CHUNKSIZE,
        low_memory=True,
    )
):
    total_in += len(chunk)

    classified = classify_lab_chunk(chunk)
    total_kept += len(classified)

    counts = classified["biomarker"].value_counts().to_dict()
    for k, v in counts.items():
        kept_counts[k] = kept_counts.get(k, 0) + int(v)

    classified.to_csv(
        OUT_FILE,
        sep="\t",
        index=False,
        mode="w" if first_write else "a",
        header=first_write,
    )

    first_write = False

    print(
        f"[INFO] Chunk {i+1}: input={len(chunk):,}, "
        f"kept={len(classified):,}, total_kept={total_kept:,}"
    )

print("\n[INFO] Finished chunked Step 4")
print("[INFO] Total input rows:", total_in)
print("[INFO] Total classified rows:", total_kept)
print("[INFO] Biomarker counts:", kept_counts)
print("[INFO] Saved:", OUT_FILE)