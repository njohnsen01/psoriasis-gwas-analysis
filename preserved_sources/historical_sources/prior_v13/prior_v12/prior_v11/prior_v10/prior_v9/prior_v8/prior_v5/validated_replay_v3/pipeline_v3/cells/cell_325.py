# ============================================================
# Block A2 — Run ancestry-stratified additive logistic models
# ============================================================

import math
import warnings
from concurrent.futures import ThreadPoolExecutor

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


# ------------------------------------------------------------
# Prespecified analysis settings
# ------------------------------------------------------------

ANCESTRIES = [
    "eur",
    "amr",
    "afr",
]

BASE_COVARIATES = [
    "sex_indicator",
    "age",
    "age2",
    *[f"PC{i}" for i in range(1, 11)],
]

MIN_CALL_RATE = 0.98
MIN_CASES = 200
MIN_CONTROLS = 500

MIN_MINOR_ALLELE_COUNT_CASES = 20.0
MIN_MINOR_ALLELE_COUNT_CONTROLS = 20.0
MIN_MINOR_ALLELE_COUNT_TOTAL = 100.0

GENOME_WIDE_THRESHOLD = 5e-8

# Four concurrent models should be efficient without placing
# excessive memory pressure on the Workbench.
MAX_WORKERS = 4


def safe_zscore(values):
    values = np.asarray(
        values,
        dtype=float,
    )

    standard_deviation = np.nanstd(
        values,
        ddof=0,
    )

    if (
        not np.isfinite(standard_deviation)
        or standard_deviation == 0
    ):
        return values - np.nanmean(values)

    return (
        values - np.nanmean(values)
    ) / standard_deviation


def fit_one_ancestry_variant(task):
    ancestry, target = task

    variant_id = target["variant_id"]
    dosage_column = target[
        "analysis_dosage_column"
    ]

    columns = [
        "psoriasis_case",
        *BASE_COVARIATES,
        dosage_column,
    ]

    frame = ancestry_analysis.loc[
        ancestry_analysis[
            "ancestry_pred"
        ].eq(ancestry),
        columns,
    ].copy()

    frame = frame.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    n_stratum = len(frame)

    n_nonmissing_dosage = int(
        frame[dosage_column]
        .notna()
        .sum()
    )

    call_rate = (
        n_nonmissing_dosage / n_stratum
        if n_stratum
        else np.nan
    )

    # Complete-case filtering is performed separately
    # for each ancestry and variant.
    frame = frame.dropna(
        subset=columns
    )

    n_total = len(frame)

    n_cases = int(
        frame["psoriasis_case"].sum()
    )

    n_controls = (
        n_total - n_cases
    )

    dosage = (
        frame[dosage_column]
        .astype(float)
        .to_numpy()
    )

    phenotype = (
        frame["psoriasis_case"]
        .astype(float)
        .to_numpy()
    )

    effect_allele_count_cases = float(
        dosage[phenotype == 1].sum()
    )

    effect_allele_count_controls = float(
        dosage[phenotype == 0].sum()
    )

    effect_allele_count_total = float(
        dosage.sum()
    )

    minor_allele_count_cases = min(
        effect_allele_count_cases,
        (
            2.0 * n_cases
            - effect_allele_count_cases
        ),
    )

    minor_allele_count_controls = min(
        effect_allele_count_controls,
        (
            2.0 * n_controls
            - effect_allele_count_controls
        ),
    )

    minor_allele_count_total = min(
        effect_allele_count_total,
        (
            2.0 * n_total
            - effect_allele_count_total
        ),
    )

    effect_allele_frequency = (
        effect_allele_count_total
        / (2.0 * n_total)
        if n_total
        else np.nan
    )

    base_result = {
        "variant_id": variant_id,
        "physical_group": target[
            "physical_group"
        ],
        "chromosome": target["CHROM"],
        "position_grch38": int(
            target["POS"]
        ),
        "effect_allele": target[
            "effect_allele"
        ],
        "other_allele": target[
            "other_allele"
        ],
        "ancestry": ancestry.upper(),
        "n_stratum_before_complete_case_filtering": (
            n_stratum
        ),
        "n_total": n_total,
        "n_cases": n_cases,
        "n_controls": n_controls,
        "call_rate": call_rate,
        "effect_allele_frequency": (
            effect_allele_frequency
        ),
        "minor_allele_count_cases": (
            minor_allele_count_cases
        ),
        "minor_allele_count_controls": (
            minor_allele_count_controls
        ),
        "minor_allele_count_total": (
            minor_allele_count_total
        ),
        "beta": np.nan,
        "se": np.nan,
        "odds_ratio": np.nan,
        "ci_lower_or": np.nan,
        "ci_upper_or": np.nan,
        "p_value": np.nan,
        "model_status": "not_estimated",
        "model_note": "",
    }

    # --------------------------------------------------------
    # Variant- and ancestry-specific QC
    # --------------------------------------------------------

    qc_failures = []

    if call_rate < MIN_CALL_RATE:
        qc_failures.append(
            f"call_rate<{MIN_CALL_RATE}"
        )

    if n_cases < MIN_CASES:
        qc_failures.append(
            f"cases<{MIN_CASES}"
        )

    if n_controls < MIN_CONTROLS:
        qc_failures.append(
            f"controls<{MIN_CONTROLS}"
        )

    if (
        minor_allele_count_cases
        < MIN_MINOR_ALLELE_COUNT_CASES
    ):
        qc_failures.append(
            "case_MAC<"
            f"{MIN_MINOR_ALLELE_COUNT_CASES:g}"
        )

    if (
        minor_allele_count_controls
        < MIN_MINOR_ALLELE_COUNT_CONTROLS
    ):
        qc_failures.append(
            "control_MAC<"
            f"{MIN_MINOR_ALLELE_COUNT_CONTROLS:g}"
        )

    if (
        minor_allele_count_total
        < MIN_MINOR_ALLELE_COUNT_TOTAL
    ):
        qc_failures.append(
            "total_MAC<"
            f"{MIN_MINOR_ALLELE_COUNT_TOTAL:g}"
        )

    dosage_standard_deviation = (
        np.nanstd(dosage)
    )

    if (
        not np.isfinite(
            dosage_standard_deviation
        )
        or dosage_standard_deviation == 0
    ):
        qc_failures.append(
            "zero_dosage_variance"
        )

    if qc_failures:
        base_result[
            "model_status"
        ] = "qc_not_estimated"

        base_result[
            "model_note"
        ] = ";".join(qc_failures)

        return base_result

    # --------------------------------------------------------
    # Construct the adjusted logistic-regression design
    # --------------------------------------------------------

    design = pd.DataFrame(
        {
            "effect_dosage": dosage,
            "sex_indicator": (
                frame["sex_indicator"]
                .astype(float)
                .to_numpy()
            ),
            "age": safe_zscore(
                frame["age"]
            ),
            "age2": safe_zscore(
                frame["age2"]
            ),
        }
    )

    # Centering age before constructing the interaction
    # preserves the model while improving conditioning.
    design["age_sex"] = (
        design["age"]
        * design["sex_indicator"]
    )

    for pc in [
        f"PC{i}"
        for i in range(1, 11)
    ]:
        design[pc] = safe_zscore(
            frame[pc]
        )

    design = sm.add_constant(
        design,
        has_constant="add",
    )

    # --------------------------------------------------------
    # Fit the ancestry-stratified additive model
    # --------------------------------------------------------

    try:
        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore"
            )

            fitted = sm.GLM(
                phenotype,
                design,
                family=(
                    sm.families.Binomial()
                ),
            ).fit(
                maxiter=100,
                disp=0,
            )

        beta = float(
            fitted.params[
                "effect_dosage"
            ]
        )

        se = float(
            fitted.bse[
                "effect_dosage"
            ]
        )

        p_value = float(
            fitted.pvalues[
                "effect_dosage"
            ]
        )

        base_result.update(
            {
                "beta": beta,
                "se": se,
                "odds_ratio": (
                    math.exp(beta)
                ),
                "ci_lower_or": (
                    math.exp(
                        beta - 1.96 * se
                    )
                ),
                "ci_upper_or": (
                    math.exp(
                        beta + 1.96 * se
                    )
                ),
                "p_value": p_value,
                "model_status": (
                    "estimated"
                    if bool(
                        getattr(
                            fitted,
                            "converged",
                            True,
                        )
                    )
                    else (
                        "estimated_"
                        "not_converged"
                    )
                ),
                "model_note": "",
            }
        )

    except Exception as error:
        base_result[
            "model_status"
        ] = "model_failed"

        base_result[
            "model_note"
        ] = (
            f"{type(error).__name__}: "
            f"{error}"
        )

    return base_result


# ------------------------------------------------------------
# Construct all 78 prespecified model tasks
# ------------------------------------------------------------

tasks = [
    (
        ancestry,
        target._asdict(),
    )
    for ancestry in ANCESTRIES
    for target in (
        ancestry_target_manifest
        .itertuples(index=False)
    )
]

expected_model_count = (
    len(ANCESTRIES)
    * len(ancestry_target_manifest)
)

if len(tasks) != 78:
    raise RuntimeError(
        "Expected 78 ancestry-stratified models "
        f"but constructed {len(tasks)}."
    )

if len(tasks) != expected_model_count:
    raise RuntimeError(
        "Model-task count does not match the "
        "target manifest."
    )


print("=" * 72)
print(
    "RUNNING TARGETED "
    "ANCESTRY-STRATIFIED MODELS"
)
print("=" * 72)

print(
    f"Models planned: {len(tasks)} "
    "(26 variants × 3 ancestry strata)"
)

print(
    f"Parallel workers: {MAX_WORKERS}"
)


# ------------------------------------------------------------
# Fit the models
# ------------------------------------------------------------

stratified_rows = []

with ThreadPoolExecutor(
    max_workers=MAX_WORKERS
) as executor:
    for completed, result in enumerate(
        executor.map(
            fit_one_ancestry_variant,
            tasks,
        ),
        start=1,
    ):
        stratified_rows.append(
            result
        )

        if (
            result["variant_id"] == QKI_ID
            or completed % 10 == 0
            or completed == len(tasks)
        ):
            if result[
                "model_status"
            ].startswith("estimated"):
                estimate_text = (
                    f"OR="
                    f"{result['odds_ratio']:.4f}, "
                    f"P="
                    f"{result['p_value']:.4g}"
                )
            else:
                estimate_text = (
                    result["model_status"]
                )

            print(
                f"[{completed:02d}/"
                f"{len(tasks)}] "
                f"{result['ancestry']} "
                f"{result['variant_id']}: "
                f"{estimate_text}"
            )


# ------------------------------------------------------------
# Assemble and correct the stratified results
# ------------------------------------------------------------

ancestry_stratified_results = (
    pd.DataFrame(
        stratified_rows
    )
)

estimated_mask = (
    ancestry_stratified_results[
        "model_status"
    ].str.startswith("estimated")
    & ancestry_stratified_results[
        "p_value"
    ].notna()
)

n_estimated_tests = int(
    estimated_mask.sum()
)

if n_estimated_tests == 0:
    raise RuntimeError(
        "No ancestry-stratified model was "
        "successfully estimated."
    )

ancestry_stratified_results[
    "targeted_bonferroni_adjusted_p"
] = np.nan

ancestry_stratified_results[
    "targeted_fdr_q_value"
] = np.nan

ancestry_stratified_results.loc[
    estimated_mask,
    "targeted_bonferroni_adjusted_p",
] = np.minimum(
    (
        ancestry_stratified_results.loc[
            estimated_mask,
            "p_value",
        ]
        * n_estimated_tests
    ),
    1.0,
)

ancestry_stratified_results.loc[
    estimated_mask,
    "targeted_fdr_q_value",
] = multipletests(
    ancestry_stratified_results.loc[
        estimated_mask,
        "p_value",
    ],
    method="fdr_bh",
)[1]

ancestry_stratified_results[
    "crosses_genome_wide_threshold"
] = (
    ancestry_stratified_results[
        "p_value"
    ]
    < GENOME_WIDE_THRESHOLD
)

ancestry_stratified_results[
    "targeted_bonferroni_significant"
] = (
    ancestry_stratified_results[
        "targeted_bonferroni_adjusted_p"
    ]
    < 0.05
)

ancestry_stratified_results[
    "targeted_fdr_significant"
] = (
    ancestry_stratified_results[
        "targeted_fdr_q_value"
    ]
    < 0.05
)


# ------------------------------------------------------------
# Validate reproduction of the previous QKI estimates
# ------------------------------------------------------------

qki_stratified_check = (
    ancestry_stratified_results.loc[
        ancestry_stratified_results[
            "variant_id"
        ].eq(QKI_ID),
        [
            "variant_id",
            "physical_group",
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
            "model_status",
        ],
    ]
    .sort_values("ancestry")
    .reset_index(drop=True)
)

expected_qki_betas = {
    "EUR": 0.196642,
    "AMR": -0.229068,
    "AFR": 0.394949,
}

for row in qki_stratified_check.itertuples(
    index=False
):
    if row.model_status.startswith(
        "estimated"
    ):
        difference = abs(
            row.beta
            - expected_qki_betas[
                row.ancestry
            ]
        )

        if difference > 0.02:
            raise RuntimeError(
                f"QKI {row.ancestry} beta "
                "did not reproduce the "
                "previously validated estimate. "
                f"Observed {row.beta}; expected "
                f"{expected_qki_betas[row.ancestry]}."
            )


# ------------------------------------------------------------
# Final results
# ------------------------------------------------------------

print("\n" + "=" * 72)
print(
    "TARGETED ANCESTRY-STRATIFIED "
    "MODELS COMPLETE"
)
print("=" * 72)

print(
    "Estimated tests:",
    n_estimated_tests,
)

print(
    "Targeted Bonferroni threshold:",
    f"{0.05 / n_estimated_tests:.6g}",
)

print(
    "Genome-wide threshold retained:",
    GENOME_WIDE_THRESHOLD,
)

print("\nModel-status summary:")
display(
    ancestry_stratified_results[
        "model_status"
    ]
    .value_counts(dropna=False)
    .rename_axis("model_status")
    .reset_index(name="n_models")
)

print("\nQKI ancestry-stratified results:")
display(
    qki_stratified_check
)

print(
    "\nThese are within-study ancestry "
    "sensitivity analyses. They are not "
    "independent replication."
)