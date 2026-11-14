## STEP 4 — Build covariate output

covar_out = covar_df.copy()

covar_out["FID"] = covar_out["person_id"]
covar_out["IID"] = covar_out["person_id"]

covar_out = covar_out[
    ["FID","IID",
     "sex","age","age2","age_sex",
     "PC1","PC2","PC3","PC4","PC5",
     "PC6","PC7","PC8","PC9","PC10"]
]

covar_out.head()