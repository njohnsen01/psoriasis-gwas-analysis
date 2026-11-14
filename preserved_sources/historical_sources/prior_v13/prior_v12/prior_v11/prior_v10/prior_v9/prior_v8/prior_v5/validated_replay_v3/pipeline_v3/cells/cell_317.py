# ============================================================
# Block 2A-R — Finish assembly and replace stale analysis_base
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd

QKI_VARIANT = "chr6:163803171:CTTATT:C"

CHR_X_TARGETS = [
    "chrX:22917619:C:T",
    "chrX:5998724:A:G",
    "chrX:67842269:G:T",
]

CHR_X_COLUMN_MAP = {
    "chrX:22917619:C:T": "effect_dosage_23",
    "chrX:5998724:A:G": "effect_dosage_22",
    "chrX:67842269:G:T": "effect_dosage_25",
}

required_covariates = [
    "age",
    "age2",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
]

# Recreate the validated sex variables.
sex_numeric = pd.to_numeric(
    covariates["sex"],
    errors="coerce",
)

sex_counts = sex_numeric.value_counts().to_dict()

if (
    sex_counts.get(0, 0) != 237917
    or sex_counts.get(1, 0) != 154238
):
    raise RuntimeError(
        f"Unexpected frozen sex counts: {sex_counts}"
    )

covariates["sex_group"] = sex_numeric.map(
    {0: "female", 1: "male"}
)
covariates["sex_indicator"] = sex_numeric.map(
    {0: 1.0, 1: 0.0}
)

for column in required_covariates:
    covariates[column] = pd.to_numeric(
        covariates[column],
        errors="coerce",
    )

covariates["age_sex"] = (
    covariates["age"]
    * covariates["sex_indicator"]
)

if "participant_key" not in covariates.columns:
    covariate_id = select_id_column(covariates)
    covariates["participant_key"] = normalize_id(
        covariates[covariate_id]
    )

covariates_for_merge = covariates[
    [
        "participant_key",
        "sex_group",
        "sex_indicator",
        "age",
        "age2",
        "age_sex",
        "PC1",
        "PC2",
        "PC3",
        "PC4",
        "PC5",
        "PC6",
        "PC7",
        "PC8",
        "PC9",
        "PC10",
    ]
].copy()

if covariates_for_merge["participant_key"].duplicated().any():
    raise RuntimeError(
        "Duplicate participant identifiers in covariates."
    )


# ------------------------------------------------------------
# Rebuild phenotype with the standardized merge key
# ------------------------------------------------------------

phenotype_source = pd.read_csv(
    PHENOTYPE_PATH,
    sep="\t",
)
phenotype_id = select_id_column(phenotype_source)

phenotype_for_merge = phenotype_source[
    [phenotype_id, "psoriasis_case"]
].copy()

phenotype_for_merge["participant_key"] = normalize_id(
    phenotype_for_merge[phenotype_id]
)
phenotype_for_merge["psoriasis_case"] = pd.to_numeric(
    phenotype_for_merge["psoriasis_case"],
    errors="coerce",
)
phenotype_for_merge = phenotype_for_merge[
    ["participant_key", "psoriasis_case"]
]

if phenotype_for_merge["participant_key"].duplicated().any():
    raise RuntimeError(
        "Duplicate participant identifiers in phenotype file."
    )


# ------------------------------------------------------------
# Ancestry
# ------------------------------------------------------------

ancestry_source = pd.read_csv(
    ANCESTRY_PATH,
    sep="\t",
)
ancestry_id = select_id_column(ancestry_source)

if "ancestry_pred" not in ancestry_source.columns:
    raise RuntimeError(
        "ancestry_pred is missing from the ancestry file."
    )

ancestry_for_merge = ancestry_source[
    [ancestry_id, "ancestry_pred"]
].copy()

ancestry_for_merge["participant_key"] = normalize_id(
    ancestry_for_merge[ancestry_id]
)
ancestry_for_merge["ancestry_pred"] = (
    ancestry_for_merge["ancestry_pred"]
    .astype("string")
    .str.strip()
    .str.lower()
)

ancestry_for_merge = ancestry_for_merge[
    ["participant_key", "ancestry_pred"]
].drop_duplicates("participant_key")


# ------------------------------------------------------------
# Corrected chromosome X effect dosages
# ------------------------------------------------------------

chrX_dosage = pd.read_csv(
    DOSAGE25_PATH,
    sep="\t",
    usecols=[
        "research_id",
        *CHR_X_COLUMN_MAP.values(),
    ],
)

chrX_dosage["participant_key"] = normalize_id(
    chrX_dosage["research_id"]
)

for variant_id, source_column in CHR_X_COLUMN_MAP.items():
    chrX_dosage[f"dosage::{variant_id}"] = pd.to_numeric(
        chrX_dosage[source_column],
        errors="coerce",
    )

chrX_dosage = chrX_dosage[
    ["participant_key"]
    + [
        f"dosage::{variant_id}"
        for variant_id in CHR_X_TARGETS
    ]
]

if chrX_dosage["participant_key"].duplicated().any():
    raise RuntimeError(
        "Duplicate participant identifiers in chrX dosage file."
    )


# ------------------------------------------------------------
# Locate and load the retained QKI dosage
# ------------------------------------------------------------

qki_search_roots = [
    Path(
        "/home/jupyter/psoriasis_phase2_reproducibility/"
        "secure_participant_level"
    ),
    Path("/home/jupyter/lead_variant_dosage_extraction"),
]

qki_candidates = []

for search_root in qki_search_roots:
    if not search_root.exists():
        continue

    for candidate_path in search_root.rglob("*"):
        if not candidate_path.is_file():
            continue

        candidate_name = candidate_path.name.lower()

        if not any(
            candidate_name.endswith(extension)
            for extension in [
                ".tsv",
                ".tsv.gz",
                ".csv",
                ".csv.gz",
                ".raw",
            ]
        ):
            continue

        try:
            candidate_header = pd.read_csv(
                candidate_path,
                sep=None,
                engine="python",
                nrows=0,
                compression="infer",
            )
        except Exception:
            continue

        qki_matches = [
            column
            for column in candidate_header.columns
            if (
                str(column) == QKI_VARIANT
                or str(column).startswith(
                    QKI_VARIANT + "_"
                )
            )
        ]

        if qki_matches:
            qki_candidates.append({
                "path": candidate_path,
                "columns": list(candidate_header.columns),
                "qki_matches": qki_matches,
            })

if not qki_candidates:
    raise RuntimeError(
        "No retained QKI dosage source was found. "
        "Do not download or recreate participant-level data."
    )

qki_candidates.sort(
    key=lambda item: (
        len(item["columns"]),
        str(item["path"]),
    )
)
qki_source = qki_candidates[0]

qki_id_candidates = [
    column
    for column in [
        "research_id",
        "person_id",
        "IID",
        "iid",
        "#IID",
    ]
    if column in qki_source["columns"]
]

if not qki_id_candidates:
    raise RuntimeError(
        "No recognized ID column in QKI dosage source."
    )

qki_id = qki_id_candidates[0]

effect_columns = [
    column
    for column in qki_source["qki_matches"]
    if str(column).endswith("_C")
]

reference_columns = [
    column
    for column in qki_source["qki_matches"]
    if str(column).endswith("_CTTATT")
]

if effect_columns:
    qki_column = effect_columns[0]
    qki_requires_flip = False
elif reference_columns:
    qki_column = reference_columns[0]
    qki_requires_flip = True
else:
    raise RuntimeError(
        "Could not determine the allele counted by the QKI "
        f"dosage column: {qki_source['qki_matches']}"
    )

qki_dosage = pd.read_csv(
    qki_source["path"],
    sep=None,
    engine="python",
    usecols=[qki_id, qki_column],
    compression="infer",
)

qki_dosage["participant_key"] = normalize_id(
    qki_dosage[qki_id]
)

qki_values = pd.to_numeric(
    qki_dosage[qki_column],
    errors="coerce",
)

if qki_requires_flip:
    qki_values = 2.0 - qki_values

qki_dosage[f"dosage::{QKI_VARIANT}"] = qki_values

qki_dosage = qki_dosage[
    [
        "participant_key",
        f"dosage::{QKI_VARIANT}",
    ]
]

if qki_dosage["participant_key"].duplicated().any():
    raise RuntimeError(
        "Duplicate participant identifiers in QKI dosage source."
    )


# ------------------------------------------------------------
# Explicitly replace the stale analysis_base
# ------------------------------------------------------------

analysis_base = (
    phenotype_for_merge
    .merge(
        covariates_for_merge,
        on="participant_key",
        how="inner",
        validate="one_to_one",
    )
    .merge(
        ancestry_for_merge,
        on="participant_key",
        how="left",
        validate="one_to_one",
    )
    .merge(
        chrX_dosage,
        on="participant_key",
        how="left",
        validate="one_to_one",
    )
    .merge(
        qki_dosage,
        on="participant_key",
        how="left",
        validate="one_to_one",
    )
)

required_analysis_columns = [
    "participant_key",
    "psoriasis_case",
    "sex_group",
    "sex_indicator",
    "age",
    "age2",
    "age_sex",
    "ancestry_pred",
    f"dosage::{QKI_VARIANT}",
    *[
        f"dosage::{variant_id}"
        for variant_id in CHR_X_TARGETS
    ],
]

missing_analysis_columns = [
    column
    for column in required_analysis_columns
    if column not in analysis_base.columns
]

if missing_analysis_columns:
    raise RuntimeError(
        f"Assembly is incomplete: {missing_analysis_columns}"
    )

n_total = len(analysis_base)
n_cases = int(analysis_base["psoriasis_case"].sum())
n_controls = n_total - n_cases

if (
    n_total != 392155
    or n_cases != 5805
    or n_controls != 386350
):
    raise RuntimeError(
        "Unexpected merged cohort totals: "
        f"{n_total:,} participants, {n_cases:,} cases, "
        f"{n_controls:,} controls."
    )

print("=" * 72)
print("CORRECTED SENSITIVITY DATASET VALIDATED")
print("=" * 72)
print("Participants:", f"{n_total:,}")
print("Cases:", f"{n_cases:,}")
print("Controls:", f"{n_controls:,}")
print(
    "Sex counts:",
    analysis_base["sex_group"].value_counts().to_dict(),
)
print(
    "Analyzable ancestry counts:",
    analysis_base.loc[
        analysis_base["ancestry_pred"].isin(
            ["eur", "amr", "afr"]
        ),
        "ancestry_pred",
    ].value_counts().to_dict(),
)
print("QKI dosage source:", qki_source["path"])
print(
    "QKI dosage conversion:",
    "CTTATT dosage converted to C dosage"
    if qki_requires_flip
    else "C dosage used directly",
)
print(
    "Participant-level rows remain in memory only and "
    "were not written to aggregate outputs."
)