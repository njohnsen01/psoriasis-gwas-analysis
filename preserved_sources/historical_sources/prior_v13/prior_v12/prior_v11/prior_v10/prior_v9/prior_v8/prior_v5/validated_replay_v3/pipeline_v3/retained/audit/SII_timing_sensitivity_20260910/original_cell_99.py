# ============================================================
# Step 6: Create derived biomarkers
# ============================================================

biomarker_df = labs_person.copy()

# Ensure expected columns exist
for col in ["LDL", "HDL", "TG", "platelet_count", "ANC", "ALC"]:
    if col not in biomarker_df.columns:
        biomarker_df[col] = np.nan

# SII-like inflammatory index = platelets × ANC / ALC
biomarker_df["inflammatory_index"] = np.where(
    (biomarker_df["platelet_count"] > 0) &
    (biomarker_df["ANC"] > 0) &
    (biomarker_df["ALC"] > 0),
    biomarker_df["platelet_count"] * biomarker_df["ANC"] / biomarker_df["ALC"],
    np.nan
)

# Optional broad sanity cap for derived SII-like index
biomarker_df.loc[
    ~np.isfinite(biomarker_df["inflammatory_index"]) |
    (biomarker_df["inflammatory_index"] <= 0),
    "inflammatory_index"
] = np.nan

# Log transforms
biomarker_df["log_TG"] = np.where(
    biomarker_df["TG"] > 0,
    np.log(biomarker_df["TG"]),
    np.nan
)

biomarker_df["log_inflammatory_index"] = np.where(
    biomarker_df["inflammatory_index"] > 0,
    np.log(biomarker_df["inflammatory_index"]),
    np.nan
)

# Z-score outcomes
outcome_cols = ["LDL", "HDL", "log_TG", "log_inflammatory_index"]

for col in outcome_cols:
    mean = biomarker_df[col].mean(skipna=True)
    sd = biomarker_df[col].std(skipna=True)
    biomarker_df[col + "_z"] = (
        (biomarker_df[col] - mean) / sd
        if pd.notna(sd) and sd > 0
        else np.nan
    )

print("[INFO] Biomarker outcomes created.")

print("\n[INFO] Non-missing outcome counts:")
print(biomarker_df[[c + "_z" for c in outcome_cols]].notna().sum())

display(
    biomarker_df[
        [
            "person_id", "LDL", "HDL", "TG", "log_TG",
            "platelet_count", "ANC", "ALC",
            "inflammatory_index", "log_inflammatory_index"
        ]
    ].head()
)

biomarker_df.to_csv(
    "psoriasis_biomarker_person_level.tsv",
    sep="\t",
    index=False
)