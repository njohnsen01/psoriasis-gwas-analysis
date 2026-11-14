# ============================================================
# Memory-safe Step 3 QC: Inspect raw biomarker lab names
# ============================================================

import pandas as pd
from collections import Counter

LAB_FILE = "psoriasis_biomarker_labs_raw.tsv"
CHUNKSIZE = 250_000

cols = [
    "measurement_concept_id",
    "lab_name",
    "lab_code",
    "unit_name",
    "unit_source_value"
]

counter = Counter()
total_rows = 0

for i, chunk in enumerate(pd.read_csv(
    LAB_FILE,
    sep="\t",
    usecols=lambda c: c in cols,
    dtype=str,
    chunksize=CHUNKSIZE,
    low_memory=True
)):
    total_rows += len(chunk)
    chunk = chunk.fillna("")

    counts = chunk.groupby(cols, dropna=False).size()

    for key, val in counts.items():
        counter[key] += int(val)

    print(f"[INFO] Finished chunk {i + 1}: rows so far = {total_rows}")

lab_counts = (
    pd.DataFrame(
        [(k[0], k[1], k[2], k[3], k[4], v) for k, v in counter.items()],
        columns=cols + ["count"]
    )
    .sort_values("count", ascending=False)
    .reset_index(drop=True)
)

print("[INFO] Total raw lab rows:", total_rows)
print("[INFO] Unique lab names:", lab_counts["lab_name"].nunique())

display(lab_counts.head(100))

lab_counts.to_csv("AoU_biomarker_lab_name_counts.tsv", sep="\t", index=False)

search_terms = ["ldl", "hdl", "trig", "platelet", "neut", "lymph", "wbc", "leukocyte"]

for term in search_terms:
    print("\n" + "=" * 70)
    print(f"[INFO] Matching lab names for: {term}")
    print("=" * 70)

    tmp = lab_counts[
        lab_counts["lab_name"].str.lower().str.contains(term, na=False) |
        lab_counts["lab_code"].str.lower().str.contains(term, na=False)
    ]

    display(tmp.head(50))

    tmp.to_csv(f"review_{term}_concepts.tsv", sep="\t", index=False)

print("[INFO] Saved AoU_biomarker_lab_name_counts.tsv and review_*_concepts.tsv files.")