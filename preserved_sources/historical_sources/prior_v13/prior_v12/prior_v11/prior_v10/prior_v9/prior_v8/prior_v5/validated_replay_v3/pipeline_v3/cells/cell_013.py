## STEP 1 — Merge phenotype + demographics + PCs
import pandas as pd

pheno_df = pheno_pc_df.copy()
demo_df  = pso_combined_cohort_demographics_df.copy()

# Convert IDs to string for REGENIE compatibility
pheno_df["person_id"] = pheno_df["person_id"].astype(str)
demo_df["person_id"]  = demo_df["person_id"].astype(str)

# Merge phenotype + demographics (PCs already inside pheno_pc_df)
merged_df = (
    pheno_df
    .merge(demo_df[["person_id","date_of_birth","sex_at_birth"]],
           on="person_id", how="left")
)

merged_df.head()