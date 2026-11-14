# ============================================================
# Block A4 — Save and upload aggregate-only ancestry results
# ============================================================

import os
import hashlib
import json
import subprocess
from datetime import date

from openpyxl import load_workbook
from openpyxl.styles import (
    Alignment,
    Font,
    PatternFill,
)


# ------------------------------------------------------------
# File helpers
# ------------------------------------------------------------

def sha256_file(path):
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(
                8 * 1024**2
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def prepare_for_export(dataframe):
    result = dataframe.copy()

    result = result.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    return result


# ------------------------------------------------------------
# Construct the external-replication template
# ------------------------------------------------------------

replication_template = (
    ancestry_stratified_results.loc[
        estimated_mask,
        [
            "variant_id",
            "physical_group",
            "chromosome",
            "position_grch38",
            "effect_allele",
            "other_allele",
            "ancestry",
            "n_total",
            "n_cases",
            "n_controls",
            "effect_allele_frequency",
            "beta",
            "se",
            "odds_ratio",
            "ci_lower_or",
            "ci_upper_or",
            "p_value",
        ],
    ]
    .copy()
    .reset_index(drop=True)
)

replication_template[
    "required_external_action"
] = (
    "Test the same effect allele in an independent "
    "cohort of the same ancestry using the "
    "prespecified covariate model; harmonize alleles "
    "before meta-analysis."
)

replication_template[
    "within_aou_replication_status"
] = "not independent replication"


# ------------------------------------------------------------
# Construct a concise QKI summary
# ------------------------------------------------------------

qki_summary = qki_stratified_check.merge(
    qki_meta_check[
        [
            "variant_id",
            "q_statistic",
            "q_df",
            "q_p",
            "i2_percent",
            "fixed_beta",
            "fixed_se",
            "fixed_or",
            "fixed_ci_lower_or",
            "fixed_ci_upper_or",
            "fixed_p",
            "random_beta",
            "random_se",
            "random_or",
            "random_ci_lower_or",
            "random_ci_upper_or",
            "random_p",
            "heterogeneity_bonferroni_adjusted_p",
            "heterogeneity_fdr_q_value",
            "heterogeneity_after_bonferroni",
            "heterogeneity_after_fdr",
        ]
    ],
    on="variant_id",
    how="left",
    validate="many_to_one",
)


# ------------------------------------------------------------
# Define aggregate-only output tables
# ------------------------------------------------------------

output_tables = {
    "target_manifest": (
        ancestry_target_manifest
    ),
    "ancestry_stratified_results": (
        ancestry_stratified_results
    ),
    "meta_heterogeneity_results": (
        ancestry_meta_results
    ),
    "leave_one_out_results": (
        ancestry_leave_one_out_results
    ),
    "ancestry_gws_results": (
        ancestry_gws_results
    ),
    "corrected_heterogeneity_results": (
        corrected_heterogeneity_results
    ),
    "qki_summary": (
        qki_summary
    ),
    "same_ancestry_replication_template": (
        replication_template
    ),
}


# ------------------------------------------------------------
# Confirm no participant identifiers can be exported
# ------------------------------------------------------------

for forbidden_column in [
    "participant_key",
    "research_id",
    "person_id",
    "FID",
    "IID",
]:
    for table_name, dataframe in (
        output_tables.items()
    ):
        if (
            forbidden_column
            in dataframe.columns
        ):
            raise RuntimeError(
                "Participant identifier column "
                f"{forbidden_column} was found in "
                f"{table_name}."
            )


# ------------------------------------------------------------
# Write aggregate TSV outputs
# ------------------------------------------------------------

ANCESTRY_OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

output_paths = {}

for table_name, dataframe in (
    output_tables.items()
):
    output_path = (
        ANCESTRY_OUTPUT_ROOT
        / f"{table_name}.tsv"
    )

    prepare_for_export(
        dataframe
    ).to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    output_paths[
        table_name
    ] = output_path


# ------------------------------------------------------------
# Build the aggregate analysis summary
# ------------------------------------------------------------

qki_meta_row = (
    qki_meta_check.iloc[0]
)

summary = {
    "analysis_date": str(
        date.today()
    ),
    "analysis_type": (
        "targeted within-study ancestry-stratified "
        "sensitivity analysis"
    ),
    "n_target_variants": int(
        len(
            ancestry_target_manifest
        )
    ),
    "target_definition": (
        "25 corrected non-MHC lead variants plus QKI"
    ),
    "ancestries": [
        ancestry.upper()
        for ancestry in ANCESTRIES
    ],
    "n_planned_stratified_tests": int(
        len(tasks)
    ),
    "n_estimated_stratified_tests": int(
        n_estimated_tests
    ),
    "n_qc_not_estimated_tests": int(
        (
            ancestry_stratified_results[
                "model_status"
            ]
            == "qc_not_estimated"
        ).sum()
    ),
    "n_valid_heterogeneity_tests": int(
        n_heterogeneity_tests
    ),
    "targeted_bonferroni_threshold": float(
        0.05
        / n_estimated_tests
    ),
    "genome_wide_threshold": float(
        GENOME_WIDE_THRESHOLD
    ),
    "qki_ancestry_strata_crossing_genome_wide_threshold": (
        qki_gws_strata
    ),
    "qki_raw_heterogeneity_p": float(
        qki_meta_row["q_p"]
    ),
    "qki_i2_percent": float(
        qki_meta_row[
            "i2_percent"
        ]
    ),
    "qki_expanded_family_bonferroni_adjusted_p": float(
        qki_meta_row[
            "heterogeneity_bonferroni_adjusted_p"
        ]
    ),
    "qki_expanded_family_fdr_q_value": float(
        qki_meta_row[
            "heterogeneity_fdr_q_value"
        ]
    ),
    "qki_heterogeneity_survives_expanded_family_correction": bool(
        qki_meta_row[
            "heterogeneity_after_bonferroni"
        ]
        or qki_meta_row[
            "heterogeneity_after_fdr"
        ]
    ),
    "n_within_study_ancestry_estimates_crossing_genome_wide_threshold": int(
        len(
            ancestry_gws_results
        )
    ),
    "within_study_ancestry_gws_variants": (
        ancestry_gws_results[
            [
                "variant_id",
                "ancestry",
            ]
        ]
        .to_dict("records")
    ),
    "n_variants_with_heterogeneity_surviving_expanded_family_correction": int(
        len(
            corrected_heterogeneity_results
        )
    ),
    "independent_replication_completed": False,
    "interpretation_guardrail": (
        "Within-study ancestry-stratified results "
        "do not constitute independent replication. "
        "Population-specific discovery requires "
        "genome-wide significance together with "
        "independent same-ancestry replication."
    ),
}

summary_path = (
    ANCESTRY_OUTPUT_ROOT
    / "ancestry_followup_summary.json"
)

summary_path.write_text(
    json.dumps(
        summary,
        indent=2,
    )
    + "\n"
)

output_paths[
    "analysis_summary"
] = summary_path


# ------------------------------------------------------------
# Create the aggregate Excel workbook
# ------------------------------------------------------------

workbook_path = (
    ANCESTRY_OUTPUT_ROOT
    / (
        "psoriasis_26_variant_ancestry_"
        "followup_aggregate_results.xlsx"
    )
)

sheet_map = {
    "Target_manifest": (
        ancestry_target_manifest
    ),
    "Stratified_results": (
        ancestry_stratified_results
    ),
    "Meta_heterogeneity": (
        ancestry_meta_results
    ),
    "Leave_one_out": (
        ancestry_leave_one_out_results
    ),
    "Ancestry_GWS": (
        ancestry_gws_results
    ),
    "Corrected_hetero": (
        corrected_heterogeneity_results
    ),
    "QKI_summary": (
        qki_summary
    ),
    "Replication_template": (
        replication_template
    ),
}

with pd.ExcelWriter(
    workbook_path,
    engine="openpyxl",
) as writer:
    for sheet_name, dataframe in (
        sheet_map.items()
    ):
        prepare_for_export(
            dataframe
        ).to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )


# ------------------------------------------------------------
# Apply readable workbook formatting
# ------------------------------------------------------------

workbook = load_workbook(
    workbook_path
)

header_fill = PatternFill(
    fill_type="solid",
    fgColor="E7E6E6",
)

header_font = Font(
    bold=True,
)

header_alignment = Alignment(
    wrap_text=True,
    vertical="top",
)

body_alignment = Alignment(
    vertical="top",
)

for worksheet in workbook.worksheets:
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = (
        worksheet.dimensions
    )

    for cell in worksheet[1]:
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = (
            header_alignment
        )

    for row in worksheet.iter_rows(
        min_row=2
    ):
        for cell in row:
            cell.alignment = (
                body_alignment
            )

    for column_cells in (
        worksheet.columns
    ):
        inspected_cells = (
            column_cells[:200]
        )

        maximum_length = max(
            (
                len(str(cell.value))
                if cell.value is not None
                else 0
            )
            for cell in inspected_cells
        )

        column_letter = (
            column_cells[0]
            .column_letter
        )

        worksheet.column_dimensions[
            column_letter
        ].width = min(
            max(
                maximum_length + 2,
                10,
            ),
            38,
        )

workbook.save(
    workbook_path
)

output_paths[
    "workbook"
] = workbook_path


# ------------------------------------------------------------
# Validate all aggregate output files
# ------------------------------------------------------------

for logical_name, path in (
    output_paths.items()
):
    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Output {logical_name} "
            f"was not created correctly: {path}"
        )


# ------------------------------------------------------------
# Build the aggregate-only manifest
# ------------------------------------------------------------

manifest_rows = []

for logical_name, path in (
    output_paths.items()
):
    manifest_rows.append(
        {
            "logical_name": (
                logical_name
            ),
            "filename": (
                path.name
            ),
            "bytes": int(
                path.stat().st_size
            ),
            "sha256": (
                sha256_file(path)
            ),
            "content_level": (
                "aggregate_variant_or_analysis_level"
            ),
            "contains_participant_identifiers": False,
            "contains_participant_rows": False,
            "safe_for_manuscript_review": True,
        }
    )

manifest = pd.DataFrame(
    manifest_rows
)

manifest_path = (
    ANCESTRY_OUTPUT_ROOT
    / "ancestry_followup_manifest.tsv"
)

manifest.to_csv(
    manifest_path,
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Upload only the aggregate output directory
# ------------------------------------------------------------

workspace_bucket = os.environ["WORKSPACE_BUCKET"].removeprefix("gs://").strip("/")
remote_destination = (
    f"gs://{workspace_bucket}/"
    "psoriasis_corrected_downstream/"
    "2026-08-01/"
    "ancestry_all_leads"
)

upload_command = [
    "gcloud",
    "storage",
    "cp",
    "--recursive",
    "--billing-project="
    f"{os.environ['GOOGLE_PROJECT']}",
    str(
        ANCESTRY_OUTPUT_ROOT
    ),
    remote_destination,
]

print(
    "Uploading aggregate-only "
    "ancestry package..."
)

subprocess.run(
    upload_command,
    check=True,
)


# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

print("\n" + "=" * 72)
print(
    "TARGETED ANCESTRY FOLLOW-UP "
    "PACKAGE COMPLETE"
)
print("=" * 72)

print(
    "Workbook:",
    workbook_path,
)

print(
    "Manifest:",
    manifest_path,
)

print(
    "Remote destination:",
    remote_destination,
)

print(
    "\nWithin-study ancestry estimates "
    "crossing P < 5e-8:"
)

if len(ancestry_gws_results):
    display(
        ancestry_gws_results
    )
else:
    print("None")

print(
    "\nQKI summary:"
)

display(
    qki_summary
)

print(
    "\nExpanded-family corrected "
    "heterogeneity results:"
)

if len(
    corrected_heterogeneity_results
):
    display(
        corrected_heterogeneity_results
    )
else:
    print("None")

print(
    "\nDownload only:"
)

print(
    "1.",
    workbook_path,
)

print(
    "2.",
    manifest_path,
)

print(
    "\nDo not download anything from "
    "private_participant_level."
)