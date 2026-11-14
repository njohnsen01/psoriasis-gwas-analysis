## STEP 2 — Build REGENIE covariate dataset

covar_df = merged_df.copy()

# Convert date_of_birth safely
covar_df["date_of_birth"] = pd.to_datetime(covar_df["date_of_birth"], errors="coerce", utc=True)

# Define CDR v8 cutoff date
cdr_cutoff = pd.Timestamp("2023-10-01", tz="UTC")

# Compute age in years
covar_df["age"] = (cdr_cutoff - covar_df["date_of_birth"]).dt.days / 365.25

# Encode sex
covar_df["sex"] = covar_df["sex_at_birth"].map({"Male":1, "Female":0})

# Add age²
covar_df["age2"] = covar_df["age"] ** 2

# Add age × sex
covar_df["age_sex"] = covar_df["age"] * covar_df["sex"]

# Keep only REGENIE-ready columns
covar_df = covar_df[[
    "person_id",
    "sex", "age", "age2", "age_sex",
    "PC1","PC2","PC3","PC4","PC5",
    "PC6","PC7","PC8","PC9","PC10"
]]

covar_df.head()