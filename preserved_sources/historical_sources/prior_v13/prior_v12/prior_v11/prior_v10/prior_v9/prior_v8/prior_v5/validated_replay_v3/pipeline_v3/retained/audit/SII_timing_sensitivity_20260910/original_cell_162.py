
# ============================================================
# Block 7 — Classify laboratory records
# ============================================================

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

LDL_IDS = {"3028288", "3028437", "3009966", "3053341", "3035899"}
HDL_IDS = {"3007070", "3053286"}
TG_IDS = {"3022192", "3022038", "3027997", "3025839"}
PLT_IDS = {"3024929", "3007461", "1497147"}
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
    for column in KEEP_COLS:
        if column not in df.columns:
            df[column] = np.nan

    df = df[KEEP_COLS].copy()
    df["measurement_concept_id"] = df["measurement_concept_id"].astype(str)
    df["unit_text"] = (
        df["unit_name"].astype(str).str.lower()
        + " | "
        + df["unit_source_value"].astype(str).str.lower()
    )
    df["biomarker"] = np.nan

    df.loc[df["measurement_concept_id"].isin(LDL_IDS), "biomarker"] = "LDL"
    df.loc[df["measurement_concept_id"].isin(HDL_IDS), "biomarker"] = "HDL"
    df.loc[df["measurement_concept_id"].isin(TG_IDS), "biomarker"] = "TG"

    absolute_unit = df["unit_text"].apply(
        lambda text: any(unit in text for unit in ABS_UNITS)
    )

    df.loc[
        df["measurement_concept_id"].isin(PLT_IDS) & absolute_unit,
        "biomarker",
    ] = "platelet_count"
    df.loc[
        df["measurement_concept_id"].isin(ANC_IDS) & absolute_unit,
        "biomarker",
    ] = "ANC"
    df.loc[
        df["measurement_concept_id"].isin(ALC_IDS) & absolute_unit,
        "biomarker",
    ] = "ALC"

    return df.dropna(subset=["biomarker"])[
        KEEP_COLS + ["biomarker"]
    ]

first_write = True
total_input = 0
total_classified = 0
classified_counts = {}

for chunk_number, chunk in enumerate(
    pd.read_csv(
        RAW_LABS,
        sep="\t",
        dtype=str,
        usecols=lambda column: column in KEEP_COLS,
        chunksize=250_000,
        low_memory=True,
    ),
    start=1,
):
    total_input += len(chunk)
    classified = classify_lab_chunk(chunk)
    total_classified += len(classified)

    for name, count in classified["biomarker"].value_counts().items():
        classified_counts[name] = (
            classified_counts.get(name, 0) + int(count)
        )

    classified.to_csv(
        CLASSIFIED_LABS,
        sep="\t",
        index=False,
        mode="w" if first_write else "a",
        header=first_write,
    )
    first_write = False

    print(
        f"Chunk {chunk_number}: input={len(chunk):,}, "
        f"kept={len(classified):,}"
    )

classification_qc = {
    "raw_rows": total_input,
    "classified_rows": total_classified,
    "classified_counts": classified_counts,
}
(PROVENANCE / "biomarker_classification_qc.json").write_text(
    json.dumps(classification_qc, indent=2),
    encoding="utf-8",
)
print(json.dumps(classification_qc, indent=2))
