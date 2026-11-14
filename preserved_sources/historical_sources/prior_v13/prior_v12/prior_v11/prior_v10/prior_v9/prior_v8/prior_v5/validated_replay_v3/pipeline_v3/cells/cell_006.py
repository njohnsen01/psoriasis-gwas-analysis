#run the phenotype construction code, where we define who has and doesnt have Pso 
import pandas as pd

# Make safe copies
demo_df = pso_combined_cohort_demographics_df.copy()
psoriasis_cond_df = pso_combined_cohort_psoriasis_conditions_df.copy()

# Ensure datetime type
psoriasis_cond_df["condition_start_datetime"] = pd.to_datetime(
    psoriasis_cond_df["condition_start_datetime"], errors="coerce"
)

# Use only the date component
psoriasis_cond_df["condition_start_date"] = (
    psoriasis_cond_df["condition_start_datetime"].dt.date
)

# Count DISTINCT psoriasis diagnosis dates per person
psoriasis_date_counts = (
    psoriasis_cond_df
    .groupby("person_id")["condition_start_date"]
    .nunique()
)

# CASES = ≥2 psoriatic diagnosis dates
case_ids = psoriasis_date_counts[psoriasis_date_counts >= 2].index

print("Number of psoriasis cases (≥2 dates):", len(case_ids))

# Build phenotype dataframe (case/control)
pso_pheno_df = (
    demo_df[["person_id"]]
    .drop_duplicates()
    .assign(
        psoriasis_case=lambda df: df["person_id"].isin(case_ids).astype(int)
    )
)

print(pso_pheno_df["psoriasis_case"].value_counts())
pso_pheno_df.head()