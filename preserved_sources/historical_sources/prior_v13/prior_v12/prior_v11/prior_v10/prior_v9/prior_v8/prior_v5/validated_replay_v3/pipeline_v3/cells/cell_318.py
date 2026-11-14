# ============================================================
# Block 2B — Corrected subgroup models and heterogeneity tests
# ============================================================

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, norm


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

DOWNSTREAM_ROOT = Path(
    "/home/jupyter/psoriasis_corrected_downstream_2026-07-31"
)
AGG = DOWNSTREAM_ROOT / "aggregate_outputs"
AGG.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "chr6:163803171:CTTATT:C": {
        "label": "QKI region",
        "effect_allele": "C",
    },
    "chrX:22917619:C:T": {
        "label": "DDX53/PTCHD1-AS region",
        "effect_allele": "T",
    },
    "chrX:5998724:A:G": {
        "label": "NLGN4X region",
        "effect_allele": "G",
    },
    "chrX:67842269:G:T": {
        "label": "YIPF6 region",
        "effect_allele": "T",
    },
}

CHR_X_TARGETS = [
    "chrX:22917619:C:T",
    "chrX:5998724:A:G",
    "chrX:67842269:G:T",
]

SEX_MODEL_COVARIATES = [
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

ANCESTRY_MODEL_COVARIATES = [
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


# ------------------------------------------------------------
# Final input validation
# ------------------------------------------------------------

required_columns = [
    "psoriasis_case",
    "sex_group",
    "ancestry_pred",
    *SEX_MODEL_COVARIATES,
    "sex_indicator",
    "age_sex",
    *[
        f"dosage::{variant_id}"
        for variant_id in TARGETS
    ],
]

missing_columns = [
    column
    for column in required_columns
    if column not in analysis_base.columns
]

if missing_columns:
    raise RuntimeError(
        f"analysis_base is missing: {missing_columns}"
    )

if len(analysis_base) != 392155:
    raise RuntimeError(
        f"Unexpected analysis-base size: {len(analysis_base):,}"
    )


# ------------------------------------------------------------
# Model function
# ------------------------------------------------------------

def fit_additive_logistic(
    dataframe,
    variant_id,
    covariate_columns,
    stratum_type,
    stratum,
):
    dosage_column = f"dosage::{variant_id}"

    required = [
        "psoriasis_case",
        dosage_column,
        *covariate_columns,
    ]

    model_data = dataframe[required].copy()

    for column in required:
        model_data[column] = pd.to_numeric(
            model_data[column],
            errors="coerce",
        )

    model_data = (
        model_data
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )

    n_total = int(len(model_data))
    n_cases = int(model_data["psoriasis_case"].sum())
    n_controls = int(n_total - n_cases)

    base_result = {
        "variant_id": variant_id,
        "region_label": TARGETS[variant_id]["label"],
        "effect_allele": TARGETS[variant_id][
            "effect_allele"
        ],
        "stratum_type": stratum_type,
        "stratum": stratum,
        "n_total": n_total,
        "n_cases": n_cases,
        "n_controls": n_controls,
    }

    if n_cases < 100 or n_controls < 1000:
        return {
            **base_result,
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "ci_lower_or": np.nan,
            "ci_upper_or": np.nan,
            "p_value": np.nan,
            "model_status": "not_estimated_underpowered",
        }

    if (
        model_data[dosage_column]
        .nunique(dropna=True)
        <= 1
    ):
        return {
            **base_result,
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "ci_lower_or": np.nan,
            "ci_upper_or": np.nan,
            "p_value": np.nan,
            "model_status": "not_estimated_no_dosage_variation",
        }

    y = model_data["psoriasis_case"].astype(float)

    X = model_data[
        [dosage_column, *covariate_columns]
    ].astype(float)

    # Remove covariates that are constant within the subgroup.
    variable_columns = [
        column
        for column in X.columns
        if X[column].nunique(dropna=True) > 1
    ]

    if dosage_column not in variable_columns:
        return {
            **base_result,
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "ci_lower_or": np.nan,
            "ci_upper_or": np.nan,
            "p_value": np.nan,
            "model_status": "not_estimated_no_dosage_variation",
        }

    X = sm.add_constant(
        X[variable_columns],
        has_constant="add",
    )

    try:
        fitted = sm.GLM(
            y,
            X,
            family=sm.families.Binomial(),
        ).fit(
            maxiter=100,
            tol=1e-8,
            disp=0,
        )

        beta = float(fitted.params[dosage_column])
        se = float(fitted.bse[dosage_column])
        p_value = float(
            fitted.pvalues[dosage_column]
        )

        return {
            **base_result,
            "beta": beta,
            "se": se,
            "odds_ratio": float(np.exp(beta)),
            "ci_lower_or": float(
                np.exp(beta - 1.96 * se)
            ),
            "ci_upper_or": float(
                np.exp(beta + 1.96 * se)
            ),
            "p_value": p_value,
            "model_status": "estimated",
        }

    except Exception as error:
        return {
            **base_result,
            "beta": np.nan,
            "se": np.nan,
            "odds_ratio": np.nan,
            "ci_lower_or": np.nan,
            "ci_upper_or": np.nan,
            "p_value": np.nan,
            "model_status": (
                "failed_" + type(error).__name__
            ),
        }


# ------------------------------------------------------------
# Sex-specific chromosome X analyses
# ------------------------------------------------------------

print("=" * 72, flush=True)
print("RUNNING SEX-SPECIFIC CHROMOSOME X MODELS", flush=True)
print("=" * 72, flush=True)

sex_rows = []

for sex_group in ["female", "male"]:
    subgroup = analysis_base.loc[
        analysis_base["sex_group"] == sex_group
    ].copy()

    subgroup_cases = int(
        subgroup["psoriasis_case"].sum()
    )
    subgroup_controls = int(
        len(subgroup) - subgroup_cases
    )

    print(
        f"\n{sex_group}: "
        f"{len(subgroup):,} participants, "
        f"{subgroup_cases:,} cases, "
        f"{subgroup_controls:,} controls",
        flush=True,
    )

    for variant_id in CHR_X_TARGETS:
        print(
            f"  Fitting {variant_id}...",
            flush=True,
        )

        result = fit_additive_logistic(
            dataframe=subgroup,
            variant_id=variant_id,
            covariate_columns=SEX_MODEL_COVARIATES,
            stratum_type="sex",
            stratum=sex_group,
        )

        sex_rows.append(result)

        print(
            "   ",
            result["model_status"],
            "OR =",
            (
                f"{result['odds_ratio']:.4f}"
                if pd.notna(result["odds_ratio"])
                else "NA"
            ),
            "P =",
            (
                f"{result['p_value']:.4g}"
                if pd.notna(result["p_value"])
                else "NA"
            ),
            flush=True,
        )

sex_results = pd.DataFrame(sex_rows)


# ------------------------------------------------------------
# Sex-heterogeneity tests
# ------------------------------------------------------------

sex_heterogeneity_rows = []

for variant_id in CHR_X_TARGETS:
    rows = sex_results.loc[
        (sex_results["variant_id"] == variant_id)
        & (sex_results["model_status"] == "estimated")
    ].set_index("stratum")

    if {"female", "male"}.issubset(rows.index):
        beta_female = float(
            rows.loc["female", "beta"]
        )
        beta_male = float(
            rows.loc["male", "beta"]
        )
        se_female = float(
            rows.loc["female", "se"]
        )
        se_male = float(
            rows.loc["male", "se"]
        )

        beta_difference = beta_female - beta_male
        se_difference = float(
            np.sqrt(se_female**2 + se_male**2)
        )
        z_heterogeneity = float(
            beta_difference / se_difference
        )
        p_heterogeneity = float(
            2.0 * norm.sf(abs(z_heterogeneity))
        )
        status = "estimated"

    else:
        beta_female = np.nan
        beta_male = np.nan
        se_female = np.nan
        se_male = np.nan
        beta_difference = np.nan
        se_difference = np.nan
        z_heterogeneity = np.nan
        p_heterogeneity = np.nan
        status = "not_estimated"

    sex_heterogeneity_rows.append({
        "variant_id": variant_id,
        "region_label": TARGETS[variant_id]["label"],
        "effect_allele": TARGETS[variant_id][
            "effect_allele"
        ],
        "beta_female": beta_female,
        "se_female": se_female,
        "beta_male": beta_male,
        "se_male": se_male,
        "beta_difference_female_minus_male": (
            beta_difference
        ),
        "se_difference": se_difference,
        "z_heterogeneity": z_heterogeneity,
        "p_heterogeneity": p_heterogeneity,
        "heterogeneity_status": status,
    })

sex_heterogeneity = pd.DataFrame(
    sex_heterogeneity_rows
)


# ------------------------------------------------------------
# Ancestry-stratified analyses
# ------------------------------------------------------------

print("\n" + "=" * 72, flush=True)
print("RUNNING ANCESTRY-STRATIFIED MODELS", flush=True)
print("=" * 72, flush=True)

ancestry_rows = []

for ancestry_group in ["eur", "amr", "afr"]:
    subgroup = analysis_base.loc[
        analysis_base["ancestry_pred"] == ancestry_group
    ].copy()

    subgroup_cases = int(
        subgroup["psoriasis_case"].sum()
    )
    subgroup_controls = int(
        len(subgroup) - subgroup_cases
    )

    print(
        f"\n{ancestry_group.upper()}: "
        f"{len(subgroup):,} participants, "
        f"{subgroup_cases:,} cases, "
        f"{subgroup_controls:,} controls",
        flush=True,
    )

    for variant_id in TARGETS:
        print(
            f"  Fitting {variant_id}...",
            flush=True,
        )

        result = fit_additive_logistic(
            dataframe=subgroup,
            variant_id=variant_id,
            covariate_columns=(
                ANCESTRY_MODEL_COVARIATES
            ),
            stratum_type="ancestry",
            stratum=ancestry_group,
        )

        ancestry_rows.append(result)

        print(
            "   ",
            result["model_status"],
            "OR =",
            (
                f"{result['odds_ratio']:.4f}"
                if pd.notna(result["odds_ratio"])
                else "NA"
            ),
            "P =",
            (
                f"{result['p_value']:.4g}"
                if pd.notna(result["p_value"])
                else "NA"
            ),
            flush=True,
        )

ancestry_results = pd.DataFrame(ancestry_rows)


# ------------------------------------------------------------
# Cochran Q ancestry-heterogeneity tests
# ------------------------------------------------------------

ancestry_heterogeneity_rows = []

for variant_id in TARGETS:
    rows = ancestry_results.loc[
        (ancestry_results["variant_id"] == variant_id)
        & (
            ancestry_results["model_status"]
            == "estimated"
        )
    ].copy()

    rows = rows.loc[
        rows["beta"].notna()
        & rows["se"].notna()
        & (rows["se"] > 0)
    ]

    k = int(len(rows))

    if k >= 2:
        weights = (
            1.0
            / np.square(rows["se"].astype(float))
        )

        fixed_effect_beta = float(
            np.sum(weights * rows["beta"])
            / np.sum(weights)
        )

        q_statistic = float(
            np.sum(
                weights
                * np.square(
                    rows["beta"]
                    - fixed_effect_beta
                )
            )
        )

        q_df = int(k - 1)
        p_heterogeneity = float(
            chi2.sf(q_statistic, q_df)
        )
        status = "estimated"

    else:
        fixed_effect_beta = np.nan
        q_statistic = np.nan
        q_df = np.nan
        p_heterogeneity = np.nan
        status = "not_estimated"

    ancestry_heterogeneity_rows.append({
        "variant_id": variant_id,
        "region_label": TARGETS[variant_id]["label"],
        "effect_allele": TARGETS[variant_id][
            "effect_allele"
        ],
        "n_ancestry_groups_estimated": k,
        "fixed_effect_beta": fixed_effect_beta,
        "q_statistic": q_statistic,
        "q_df": q_df,
        "p_heterogeneity": p_heterogeneity,
        "heterogeneity_status": status,
    })

ancestry_heterogeneity = pd.DataFrame(
    ancestry_heterogeneity_rows
)


# ------------------------------------------------------------
# Output validation
# ------------------------------------------------------------

expected_dimensions = {
    "sex_results": (len(sex_results), 6),
    "sex_heterogeneity": (
        len(sex_heterogeneity),
        3,
    ),
    "ancestry_results": (
        len(ancestry_results),
        12,
    ),
    "ancestry_heterogeneity": (
        len(ancestry_heterogeneity),
        4,
    ),
}

incorrect_dimensions = {
    name: observed
    for name, (observed, expected)
    in expected_dimensions.items()
    if observed != expected
}

if incorrect_dimensions:
    raise RuntimeError(
        "Unexpected result dimensions: "
        f"{incorrect_dimensions}"
    )

sex_models_estimated = int(
    (sex_results["model_status"] == "estimated").sum()
)
ancestry_models_estimated = int(
    (
        ancestry_results["model_status"]
        == "estimated"
    ).sum()
)

if sex_models_estimated != 6:
    print(
        "WARNING: Not all sex-specific models were estimated.",
        flush=True,
    )

if ancestry_models_estimated != 12:
    print(
        "WARNING: Not all ancestry-specific models were "
        "estimated.",
        flush=True,
    )


# ------------------------------------------------------------
# Save aggregate-only outputs
# ------------------------------------------------------------

SEX_RESULTS_PATH = (
    AGG
    / "corrected_priority_sex_stratified_results.tsv"
)
SEX_HETEROGENEITY_PATH = (
    AGG
    / "corrected_priority_sex_heterogeneity.tsv"
)
ANCESTRY_RESULTS_PATH = (
    AGG
    / "corrected_priority_ancestry_stratified_results.tsv"
)
ANCESTRY_HETEROGENEITY_PATH = (
    AGG
    / "corrected_priority_ancestry_heterogeneity.tsv"
)
SENSITIVITY_SUMMARY_PATH = (
    AGG
    / "corrected_priority_sensitivity_summary.json"
)

sex_results.to_csv(
    SEX_RESULTS_PATH,
    sep="\t",
    index=False,
)
sex_heterogeneity.to_csv(
    SEX_HETEROGENEITY_PATH,
    sep="\t",
    index=False,
)
ancestry_results.to_csv(
    ANCESTRY_RESULTS_PATH,
    sep="\t",
    index=False,
)
ancestry_heterogeneity.to_csv(
    ANCESTRY_HETEROGENEITY_PATH,
    sep="\t",
    index=False,
)

summary = {
    "status": "complete",
    "created_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "analysis": (
        "corrected priority sex- and "
        "ancestry-stratified sensitivity analysis"
    ),
    "participants": int(len(analysis_base)),
    "cases": int(
        analysis_base["psoriasis_case"].sum()
    ),
    "controls": int(
        len(analysis_base)
        - analysis_base["psoriasis_case"].sum()
    ),
    "sex_specific_models_expected": 6,
    "sex_specific_models_estimated": (
        sex_models_estimated
    ),
    "sex_heterogeneity_tests_estimated": int(
        (
            sex_heterogeneity[
                "heterogeneity_status"
            ]
            == "estimated"
        ).sum()
    ),
    "ancestry_specific_models_expected": 12,
    "ancestry_specific_models_estimated": (
        ancestry_models_estimated
    ),
    "ancestry_heterogeneity_tests_estimated": int(
        (
            ancestry_heterogeneity[
                "heterogeneity_status"
            ]
            == "estimated"
        ).sum()
    ),
    "ancestry_groups": ["eur", "amr", "afr"],
    "sex_encoding": "0=female; 1=male",
    "participant_rows_written_to_aggregate_outputs": (
        False
    ),
}

SENSITIVITY_SUMMARY_PATH.write_text(
    json.dumps(summary, indent=2)
)


# ------------------------------------------------------------
# Display aggregate results only
# ------------------------------------------------------------

print("\n" + "=" * 72)
print("CORRECTED PRIORITY SENSITIVITY MODELS COMPLETE")
print("=" * 72)

print("\nSex-specific estimates:")
print(
    sex_results[
        [
            "region_label",
            "variant_id",
            "stratum",
            "n_cases",
            "n_controls",
            "odds_ratio",
            "ci_lower_or",
            "ci_upper_or",
            "p_value",
            "model_status",
        ]
    ].to_string(index=False)
)

print("\nSex heterogeneity:")
print(
    sex_heterogeneity[
        [
            "region_label",
            "variant_id",
            "z_heterogeneity",
            "p_heterogeneity",
            "heterogeneity_status",
        ]
    ].to_string(index=False)
)

print("\nAncestry-specific estimates:")
print(
    ancestry_results[
        [
            "region_label",
            "variant_id",
            "stratum",
            "n_cases",
            "n_controls",
            "odds_ratio",
            "ci_lower_or",
            "ci_upper_or",
            "p_value",
            "model_status",
        ]
    ].to_string(index=False)
)

print("\nAncestry heterogeneity:")
print(
    ancestry_heterogeneity[
        [
            "region_label",
            "variant_id",
            "n_ancestry_groups_estimated",
            "q_statistic",
            "q_df",
            "p_heterogeneity",
            "heterogeneity_status",
        ]
    ].to_string(index=False)
)

print("\nSaved aggregate outputs:")
for output_path in [
    SEX_RESULTS_PATH,
    SEX_HETEROGENEITY_PATH,
    ANCESTRY_RESULTS_PATH,
    ANCESTRY_HETEROGENEITY_PATH,
    SENSITIVITY_SUMMARY_PATH,
]:
    print(output_path)

print(
    "\nNo participant identifiers or participant-level rows "
    "were written to aggregate outputs."
)