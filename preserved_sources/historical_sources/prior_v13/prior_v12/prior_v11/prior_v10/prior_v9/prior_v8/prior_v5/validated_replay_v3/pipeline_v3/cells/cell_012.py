## Merge PCs onto your phenotype data
pso_pheno_df = pso_pheno_df.copy()
demo_df = pso_combined_cohort_demographics_df.copy()
pc_df = pc_df.copy()

pso_pheno_df["person_id"] = pso_pheno_df["person_id"].astype(int)
demo_df["person_id"] = demo_df["person_id"].astype(int)
pc_df["research_id"] = pc_df["research_id"].astype(int)

pheno_pc_df = (
    pso_pheno_df
    .merge(pc_df, left_on="person_id", right_on="research_id", how="inner")
)

pheno_pc_df[["person_id", "psoriasis_case", "ancestry_pred", "PC1","PC2"]].head()