## STEP 3 — Build phenotype file

phenotype_out = merged_df[["person_id", "psoriasis_case"]].copy()

phenotype_out["FID"] = phenotype_out["person_id"]
phenotype_out["IID"] = phenotype_out["person_id"]

phenotype_out = phenotype_out[["FID","IID","psoriasis_case"]]

phenotype_out.head()