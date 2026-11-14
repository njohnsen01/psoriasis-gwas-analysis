# ============================================================
# Block A1 — Build ancestry-analysis matrix for 25 leads + QKI
# ============================================================

import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
from IPython.display import display


PRIMARY_ROOT = Path(
    "/home/jupyter/psoriasis_regenie_loco_step2/"
    f"attempt_{os.environ['PSO_PRIMARY_RUN_ID']}"
)

# This earlier file contains 14 broader groups. It is retained only
# as a source check and will not be modified.
LEAD_GROUPS_PATH = (
    PRIMARY_ROOT
    / "final_aggregate/"
    "psoriasis_corrected_non_mhc_lead_groups.tsv"
)

# Private participant-level dosages. These remain inside Workbench.
LEAD_DOSAGE_PATH = Path(
    "/home/jupyter/psoriasis_corrected_downstream_2026-07-31/"
    "private_participant_level/"
    "corrected_25_lead_dosages_effect_allele.tsv.gz"
)

# Aggregate-only file containing the locked 25-variant crosswalk.
LEAD_BIOMARKER_AGGREGATE_PATH = Path(
    "/home/jupyter/psoriasis_corrected_downstream_2026-07-31/"
    "aggregate_outputs/"
    "corrected_25_lead_biomarker_results_full_cohort.tsv"
)

ANCESTRY_OUTPUT_ROOT = Path(
    "/home/jupyter/psoriasis_corrected_downstream_2026-08-01/"
    "ancestry_all_leads_aggregate_outputs"
)
ANCESTRY_OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

QKI_ID = "chr6:163803171:CTTATT:C"
QKI_DOSAGE = f"dosage::{QKI_ID}"
QKI_EFFECT_ALLELE = "C"
QKI_OTHER_ALLELE = "CTTATT"


def normalize_identifier(values):
    return (
        values.astype("string")
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


# ------------------------------------------------------------
# Validate the previously prepared analysis_base object
# ------------------------------------------------------------

if "analysis_base" not in globals():
    raise RuntimeError(
        "analysis_base is not present. Rerun the previously "
        "validated corrected sensitivity-dataset block first."
    )

required_base_columns = [
    "participant_key",
    "psoriasis_case",
    "sex_indicator",
    "ancestry_pred",
    "age",
    "age2",
    *[f"PC{i}" for i in range(1, 11)],
    QKI_DOSAGE,
]

missing_base_columns = [
    column
    for column in required_base_columns
    if column not in analysis_base.columns
]

if missing_base_columns:
    raise RuntimeError(
        "analysis_base is missing required columns: "
        f"{missing_base_columns}"
    )


# ------------------------------------------------------------
# Validate the required dosage and aggregate-only files
# ------------------------------------------------------------

for required_path in [
    LEAD_DOSAGE_PATH,
    LEAD_BIOMARKER_AGGREGATE_PATH,
]:
    if (
        not required_path.exists()
        or required_path.stat().st_size == 0
    ):
        raise FileNotFoundError(required_path)


# ------------------------------------------------------------
# Read the earlier lead-group file, if available
# ------------------------------------------------------------

if (
    LEAD_GROUPS_PATH.exists()
    and LEAD_GROUPS_PATH.stat().st_size > 0
):
    leads = pd.read_csv(
        LEAD_GROUPS_PATH,
        sep="\t",
        dtype={"CHROM": str},
    )
else:
    leads = pd.DataFrame()

lead_source = str(LEAD_GROUPS_PATH)

print(
    "Initial lead-group rows:",
    len(leads),
)


# ------------------------------------------------------------
# Recover the locked 25 lead-centered variants
# ------------------------------------------------------------
# The earlier lead-group TSV contains 14 broader physical groups.
# The frozen 25-column dosage matrix follows the later locked
# 25 lead-centered crosswalk. Recover that exact crosswalk from
# the aggregate biomarker-results file.

if (
    len(leads) != 25
    or "ID" not in leads.columns
    or leads["ID"].nunique() != 25
):
    print(
        "The earlier file does not contain the locked "
        "25-row crosswalk."
    )
    print(
        "Recovering the 25 corrected lead variants from:"
    )
    print(LEAD_BIOMARKER_AGGREGATE_PATH)

    biomarker_leads = pd.read_csv(
        LEAD_BIOMARKER_AGGREGATE_PATH,
        sep="\t",
        dtype={"chromosome": str},
    )

    required_biomarker_columns = {
        "physical_group",
        "chromosome",
        "position",
        "variant_id",
        "reference_allele",
        "effect_allele",
        "corrected_gwas_beta",
        "corrected_gwas_or",
        "corrected_gwas_p",
    }

    missing_biomarker_columns = (
        required_biomarker_columns
        - set(biomarker_leads.columns)
    )

    if missing_biomarker_columns:
        raise RuntimeError(
            "The 25-lead aggregate biomarker file is missing: "
            f"{sorted(missing_biomarker_columns)}"
        )

    biomarker_leads = (
        biomarker_leads[
            [
                "physical_group",
                "chromosome",
                "position",
                "variant_id",
                "reference_allele",
                "effect_allele",
                "corrected_gwas_beta",
                "corrected_gwas_or",
                "corrected_gwas_p",
            ]
        ]
        .drop_duplicates("variant_id")
        .reset_index(drop=True)
    )

    if len(biomarker_leads) != 25:
        raise RuntimeError(
            "Expected 25 unique corrected lead variants in "
            "the aggregate biomarker file; observed "
            f"{len(biomarker_leads)}."
        )

    leads = biomarker_leads.rename(
        columns={
            "chromosome": "CHROM",
            "position": "POS",
            "variant_id": "ID",
            "reference_allele": "REF",
            "effect_allele": "ALT",
            "corrected_gwas_beta": "BETA",
            "corrected_gwas_or": "OR",
            "corrected_gwas_p": "P",
        }
    )

    # These values are not required for the ancestry models.
    leads["SE"] = np.nan
    leads["A1FREQ"] = np.nan

    lead_source = str(
        LEAD_BIOMARKER_AGGREGATE_PATH
    )


# ------------------------------------------------------------
# Validate the recovered 25-row lead table
# ------------------------------------------------------------

required_lead_columns = {
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "BETA",
    "SE",
    "OR",
    "P",
    "A1FREQ",
    "physical_group",
}

missing_lead_columns = (
    required_lead_columns
    - set(leads.columns)
)

if missing_lead_columns:
    raise RuntimeError(
        "The recovered lead table is missing columns: "
        f"{sorted(missing_lead_columns)}"
    )

if (
    len(leads) != 25
    or leads["ID"].nunique() != 25
):
    raise RuntimeError(
        "Expected 25 unique corrected lead variants; "
        f"observed {len(leads)} rows and "
        f"{leads['ID'].nunique()} unique IDs."
    )


# ------------------------------------------------------------
# Put the lead variants into frozen genomic order
# ------------------------------------------------------------

chromosome_order = {
    str(value): value
    for value in range(1, 23)
}
chromosome_order["X"] = 23

leads["CHROM"] = (
    leads["CHROM"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
)

leads["_chromosome_order"] = (
    leads["CHROM"].map(chromosome_order)
)

if leads["_chromosome_order"].isna().any():
    unexpected = (
        leads.loc[
            leads["_chromosome_order"].isna(),
            "CHROM",
        ]
        .unique()
        .tolist()
    )
    raise RuntimeError(
        "Unexpected chromosome labels in the corrected "
        f"lead table: {unexpected}"
    )

leads = (
    leads.sort_values(
        [
            "_chromosome_order",
            "POS",
            "ID",
        ]
    )
    .drop(columns="_chromosome_order")
    .reset_index(drop=True)
)


# ------------------------------------------------------------
# Validate dosage-file structure
# ------------------------------------------------------------

dosage_header = pd.read_csv(
    LEAD_DOSAGE_PATH,
    sep="\t",
    nrows=0,
).columns.tolist()

if "research_id" not in dosage_header:
    raise RuntimeError(
        "Expected research_id in the 25-lead dosage file."
    )

dosage_columns = sorted(
    [
        column
        for column in dosage_header
        if re.fullmatch(
            r"effect_dosage_\d{2}",
            column,
        )
    ],
    key=lambda value: int(
        value.rsplit("_", 1)[1]
    ),
)

expected_dosage_columns = [
    f"effect_dosage_{index:02d}"
    for index in range(1, 26)
]

if dosage_columns != expected_dosage_columns:
    raise RuntimeError(
        "Unexpected dosage columns: "
        f"{dosage_columns}"
    )


# ------------------------------------------------------------
# Construct the frozen dosage-to-variant crosswalk
# ------------------------------------------------------------

lead_dosage_map = leads[
    [
        "physical_group",
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "BETA",
        "SE",
        "OR",
        "P",
        "A1FREQ",
    ]
].copy()

lead_dosage_map["dosage_source_column"] = (
    dosage_columns
)

lead_dosage_map["analysis_dosage_column"] = (
    "dosage::"
    + lead_dosage_map["ID"].astype(str)
)


# ------------------------------------------------------------
# Validate known chromosome X mapping
# ------------------------------------------------------------

known_mapping = {
    "effect_dosage_22": "chrX:5998724:A:G",
    "effect_dosage_23": "chrX:22917619:C:T",
    "effect_dosage_24": "chrX:67736589:T:C",
    "effect_dosage_25": "chrX:67842269:G:T",
}

observed_mapping = dict(
    zip(
        lead_dosage_map[
            "dosage_source_column"
        ],
        lead_dosage_map["ID"],
    )
)

mapping_mismatches = {
    dosage_column: (
        expected_variant,
        observed_mapping.get(dosage_column),
    )
    for dosage_column, expected_variant
    in known_mapping.items()
    if (
        observed_mapping.get(dosage_column)
        != expected_variant
    )
}

if mapping_mismatches:
    raise RuntimeError(
        "The frozen dosage-column mapping did not "
        f"validate: {mapping_mismatches}"
    )


# ------------------------------------------------------------
# Read and rename the 25 participant-level dosages
# ------------------------------------------------------------
# These rows remain inside Workbench and are not exported.

lead_dosages = pd.read_csv(
    LEAD_DOSAGE_PATH,
    sep="\t",
    usecols=[
        "research_id",
        *dosage_columns,
    ],
    dtype={"research_id": "string"},
)

lead_dosages["participant_key"] = (
    normalize_identifier(
        lead_dosages.pop("research_id")
    )
)

if lead_dosages[
    "participant_key"
].duplicated().any():
    raise RuntimeError(
        "Duplicate participant keys in the "
        "25-lead dosage file."
    )

lead_dosages = lead_dosages.rename(
    columns=dict(
        zip(
            lead_dosage_map[
                "dosage_source_column"
            ],
            lead_dosage_map[
                "analysis_dosage_column"
            ],
        )
    )
)


# ------------------------------------------------------------
# Merge dosages with the validated analysis dataset
# ------------------------------------------------------------

ancestry_analysis = analysis_base[
    required_base_columns
].copy()

ancestry_analysis["participant_key"] = (
    normalize_identifier(
        ancestry_analysis["participant_key"]
    )
)

ancestry_analysis["ancestry_pred"] = (
    ancestry_analysis["ancestry_pred"]
    .astype("string")
    .str.strip()
    .str.lower()
)

if ancestry_analysis[
    "participant_key"
].duplicated().any():
    raise RuntimeError(
        "Duplicate participant keys in analysis_base."
    )

ancestry_analysis = ancestry_analysis.merge(
    lead_dosages,
    on="participant_key",
    how="left",
    validate="one_to_one",
)

if len(ancestry_analysis) != 392155:
    raise RuntimeError(
        "Expected 392,155 participants after merging; "
        f"observed {len(ancestry_analysis):,}."
    )

observed_case_count = int(
    ancestry_analysis["psoriasis_case"].sum()
)

if observed_case_count != 5805:
    raise RuntimeError(
        "Psoriasis case count changed during the "
        f"dosage merge: {observed_case_count:,}."
    )


# ------------------------------------------------------------
# Create the 25-lead target manifest
# ------------------------------------------------------------

lead_target_manifest = (
    lead_dosage_map.rename(
        columns={
            "ID": "variant_id",
            "REF": "other_allele",
            "ALT": "effect_allele",
            "BETA": "primary_beta",
            "SE": "primary_se",
            "OR": "primary_or",
            "P": "primary_p",
            "A1FREQ": (
                "primary_effect_allele_frequency"
            ),
        }
    )[
        [
            "variant_id",
            "physical_group",
            "CHROM",
            "POS",
            "effect_allele",
            "other_allele",
            "analysis_dosage_column",
            "dosage_source_column",
            "primary_beta",
            "primary_se",
            "primary_or",
            "primary_p",
            "primary_effect_allele_frequency",
        ]
    ]
)


# ------------------------------------------------------------
# Add QKI as the 26th targeted follow-up variant
# ------------------------------------------------------------

qki_manifest = pd.DataFrame(
    [
        {
            "variant_id": QKI_ID,
            "physical_group": "QKI follow-up",
            "CHROM": "6",
            "POS": 163803171,
            "effect_allele": (
                QKI_EFFECT_ALLELE
            ),
            "other_allele": (
                QKI_OTHER_ALLELE
            ),
            "analysis_dosage_column": (
                QKI_DOSAGE
            ),
            "dosage_source_column": (
                "previously validated QKI dosage"
            ),
            "primary_beta": 0.186252,
            "primary_se": 0.034079,
            "primary_or": float(
                np.exp(0.186252)
            ),
            "primary_p": 9.931389e-8,
            "primary_effect_allele_frequency": (
                0.072297
            ),
        }
    ]
)

ancestry_target_manifest = pd.concat(
    [
        lead_target_manifest,
        qki_manifest,
    ],
    ignore_index=True,
)

if (
    len(ancestry_target_manifest) != 26
    or ancestry_target_manifest[
        "variant_id"
    ].nunique() != 26
):
    raise RuntimeError(
        "Expected 25 corrected lead variants plus "
        "one distinct QKI variant."
    )


# ------------------------------------------------------------
# Confirm every dosage column is present
# ------------------------------------------------------------

analysis_dosage_columns = (
    ancestry_target_manifest[
        "analysis_dosage_column"
    ].tolist()
)

missing_dosage_columns = [
    column
    for column in analysis_dosage_columns
    if column not in ancestry_analysis.columns
]

if missing_dosage_columns:
    raise RuntimeError(
        "Merged analysis matrix is missing dosage "
        f"columns: {missing_dosage_columns}"
    )


# ------------------------------------------------------------
# Summarize eligible ancestry strata
# ------------------------------------------------------------

ancestry_counts = (
    ancestry_analysis.loc[
        ancestry_analysis[
            "ancestry_pred"
        ].isin(
            [
                "eur",
                "amr",
                "afr",
            ]
        )
    ]
    .groupby(
        "ancestry_pred",
        observed=True,
    )["psoriasis_case"]
    .agg(
        n_total="size",
        n_cases="sum",
    )
)

ancestry_counts["n_controls"] = (
    ancestry_counts["n_total"]
    - ancestry_counts["n_cases"]
)


# ------------------------------------------------------------
# Final validation report
# ------------------------------------------------------------

print("=" * 72)
print("ANCESTRY FOLLOW-UP MATRIX VALIDATED")
print("=" * 72)

print(
    "Targets: 25 corrected non-MHC "
    "lead variants plus QKI"
)
print(
    "Locked 25-lead crosswalk source:",
    lead_source,
)
print(
    "Participants:",
    f"{len(ancestry_analysis):,}",
)
print(
    "Cases:",
    f"{observed_case_count:,}",
)

print("\nEligible ancestry counts:")
display(ancestry_counts.reset_index())

print("\nValidated chromosome X dosage mapping:")
for dosage_column, variant_id in known_mapping.items():
    print(
        f"  {dosage_column} -> {variant_id}"
    )

print(
    "\nParticipant-level rows remain in memory only "
    "and will not be written to aggregate outputs."
)