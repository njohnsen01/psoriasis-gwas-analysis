# ============================================================
# Step 7: Merge phenotype/covariates + biomarkers
# ============================================================

import pandas as pd
import numpy as np

# Reload if kernel restarted
if "analysis_base" not in globals():
    pheno = pd.read_csv("psoriasis_pheno.tsv", sep="\t", dtype={"FID": str, "IID": str})
    covar = pd.read_csv("psoriasis_covariates.tsv", sep="\t", dtype={"FID": str, "IID": str})
    analysis_base = pheno.merge(covar, on=["FID", "IID"], how="inner")
    analysis_base["person_id"] = analysis_base["IID"].astype(str)

if "biomarker_df" not in globals():
    biomarker_df = pd.read_csv(
        "psoriasis_biomarker_person_level.tsv",
        sep="\t",
        dtype={"person_id": str},
        low_memory=False
    )

# Ensure 1 row/person before merge
analysis_base["person_id"] = analysis_base["person_id"].astype(str)
biomarker_df["person_id"] = biomarker_df["person_id"].astype(str)

analysis_base = analysis_base.drop_duplicates(subset=["person_id"]).copy()
biomarker_df = biomarker_df.drop_duplicates(subset=["person_id"]).copy()

outcomes = ["LDL_z", "HDL_z", "log_TG_z", "log_inflammatory_index_z"]

analysis_biomarker_base = analysis_base.merge(
    biomarker_df,
    on="person_id",
    how="left",
    validate="one_to_one"
)

print("[INFO] analysis_biomarker_base shape:", analysis_biomarker_base.shape)
print("[INFO] Unique participants:", analysis_biomarker_base["person_id"].nunique())

print("\n[INFO] Psoriasis case counts:")
print(analysis_biomarker_base["psoriasis_case"].value_counts(dropna=False))

print("\n[INFO] Non-missing biomarker counts:")
print(analysis_biomarker_base[outcomes].notna().sum())

print("\n[INFO] Non-missing among psoriasis cases:")
print(
    analysis_biomarker_base.loc[
        analysis_biomarker_base["psoriasis_case"] == 1,
        outcomes
    ].notna().sum()
)

analysis_biomarker_base.to_csv(
    "psoriasis_biomarker_analysis_base_no_genotypes.tsv",
    sep="\t",
    index=False
)

display(analysis_biomarker_base.head())