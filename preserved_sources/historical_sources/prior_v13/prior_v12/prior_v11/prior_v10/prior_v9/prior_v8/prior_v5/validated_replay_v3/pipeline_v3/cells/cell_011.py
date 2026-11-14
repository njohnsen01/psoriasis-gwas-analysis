import ast

# Start from your ancestry_df
pc_work = ancestry_df.copy()

# Parse the stringified list in pca_features into real Python lists
pc_work["pca_features"] = pc_work["pca_features"].apply(ast.literal_eval)

# Turn the first 10 elements of pca_features into separate PC columns
pc_cols = (
    pc_work["pca_features"]
    .apply(lambda x: pd.Series(x[:10]))
    .rename(columns={i: f"PC{i+1}" for i in range(10)})
)

pc_df = pd.concat(
    [
        pc_work[["research_id", "ancestry_pred", "ancestry_pred_other"]],
        pc_cols,
    ],
    axis=1,
)

pc_df.head()