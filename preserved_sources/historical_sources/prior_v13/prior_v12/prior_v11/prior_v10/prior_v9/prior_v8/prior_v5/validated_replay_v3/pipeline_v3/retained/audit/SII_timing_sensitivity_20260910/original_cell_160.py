# ============================================================
# Block 5 — Locate retained biomarker checkpoints
# ============================================================

checkpoint_names = [
    "psoriasis_biomarker_labs_raw.tsv",
    "psoriasis_biomarker_labs_classified.tsv",
    "psoriasis_biomarker_person_level.tsv",
]

checkpoint_paths = {}

for filename in checkpoint_names:
    path = copy_first_bucket_match(
        filename,
        SECURE / filename,
        required=False,
    )
    checkpoint_paths[filename] = path

for filename, path in checkpoint_paths.items():
    print(filename, "=>", path)

RAW_LABS = SECURE / "psoriasis_biomarker_labs_raw.tsv"
CLASSIFIED_LABS = SECURE / "psoriasis_biomarker_labs_classified.tsv"
PERSON_BIOMARKERS = SECURE / "psoriasis_biomarker_person_level.tsv"
