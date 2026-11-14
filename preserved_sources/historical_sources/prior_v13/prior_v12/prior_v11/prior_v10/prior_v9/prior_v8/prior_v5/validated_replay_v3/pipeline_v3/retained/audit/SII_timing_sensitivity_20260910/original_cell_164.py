# ============================================================
# Block 9 — Derive and standardize the four outcomes
# ============================================================

biomarker_df = labs_person.copy()

for column in ["LDL", "HDL", "TG", "platelet_count", "ANC", "ALC"]:
    if column not in biomarker_df.columns:
        biomarker_df[column] = np.nan

biomarker_df["inflammatory_index"] = np.where(
    (biomarker_df["platelet_count"] > 0)
    & (biomarker_df["ANC"] > 0)
    & (biomarker_df["ALC"] > 0),
    biomarker_df["platelet_count"]
    * biomarker_df["ANC"]
    / biomarker_df["ALC"],
    np.nan,
)

biomarker_df.loc[
    ~np.isfinite(biomarker_df["inflammatory_index"])
    | (biomarker_df["inflammatory_index"] <= 0),
    "inflammatory_index",
] = np.nan

biomarker_df["log_TG"] = np.where(
    biomarker_df["TG"] > 0,
    np.log(biomarker_df["TG"]),
    np.nan,
)
biomarker_df["log_inflammatory_index"] = np.where(
    biomarker_df["inflammatory_index"] > 0,
    np.log(biomarker_df["inflammatory_index"]),
    np.nan,
)

for source in ["LDL", "HDL", "log_TG", "log_inflammatory_index"]:
    mean = biomarker_df[source].mean(skipna=True)
    sd = biomarker_df[source].std(skipna=True, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise RuntimeError(f"Invalid standard deviation for {source}: {sd}")
    biomarker_df[f"{source}_z"] = (biomarker_df[source] - mean) / sd

biomarker_df.to_csv(
    PERSON_BIOMARKERS,
    sep="\t",
    index=False,
)

OUTCOMES = [
    "HDL_z",
    "LDL_z",
    "log_TG_z",
    "log_inflammatory_index_z",
]

print("Nonmissing standardized outcomes:")
print(biomarker_df[OUTCOMES].notna().sum().to_dict())
