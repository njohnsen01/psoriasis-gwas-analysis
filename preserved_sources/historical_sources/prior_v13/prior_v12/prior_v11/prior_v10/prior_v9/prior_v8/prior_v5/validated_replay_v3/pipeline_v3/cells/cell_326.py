# ============================================================
# Block A3 — Meta-analysis, heterogeneity, and pattern audit
# ============================================================

from scipy.stats import chi2, norm


meta_rows = []
leave_one_out_rows = []


# ------------------------------------------------------------
# Analyze each of the 25 corrected leads plus QKI
# ------------------------------------------------------------

for target in ancestry_target_manifest.to_dict(
    "records"
):
    variant_id = target["variant_id"]

    valid = ancestry_stratified_results.loc[
        ancestry_stratified_results[
            "variant_id"
        ].eq(variant_id)
        & ancestry_stratified_results[
            "model_status"
        ].str.startswith("estimated")
        & ancestry_stratified_results[
            "beta"
        ].notna()
        & ancestry_stratified_results[
            "se"
        ].gt(0)
    ].copy()

    # --------------------------------------------------------
    # Record variants with fewer than two valid strata
    # --------------------------------------------------------

    if len(valid) < 2:
        meta_rows.append(
            {
                "variant_id": variant_id,
                "physical_group": target[
                    "physical_group"
                ],
                "effect_allele": target[
                    "effect_allele"
                ],
                "other_allele": target[
                    "other_allele"
                ],
                "n_ancestries": len(valid),
                "ancestries_included": ",".join(
                    valid["ancestry"]
                ),
                "fixed_beta": np.nan,
                "fixed_se": np.nan,
                "fixed_or": np.nan,
                "fixed_ci_lower_or": np.nan,
                "fixed_ci_upper_or": np.nan,
                "fixed_p": np.nan,
                "q_statistic": np.nan,
                "q_df": np.nan,
                "q_p": np.nan,
                "i2_percent": np.nan,
                "tau2": np.nan,
                "random_beta": np.nan,
                "random_se": np.nan,
                "random_or": np.nan,
                "random_ci_lower_or": np.nan,
                "random_ci_upper_or": np.nan,
                "random_p": np.nan,
                "ancestries_crossing_genome_wide_threshold": "",
                "within_study_interpretation": (
                    "Insufficient valid ancestry estimates"
                ),
            }
        )
        continue

    # --------------------------------------------------------
    # Fixed-effects meta-analysis
    # --------------------------------------------------------

    beta = valid["beta"].to_numpy(
        dtype=float
    )

    se = valid["se"].to_numpy(
        dtype=float
    )

    weights = 1.0 / np.square(se)

    fixed_beta = float(
        np.sum(weights * beta)
        / np.sum(weights)
    )

    fixed_se = float(
        np.sqrt(
            1.0 / np.sum(weights)
        )
    )

    fixed_z = (
        fixed_beta / fixed_se
    )

    fixed_p = float(
        2.0
        * norm.sf(
            abs(fixed_z)
        )
    )

    # --------------------------------------------------------
    # Cochran Q and I-squared
    # --------------------------------------------------------

    q_statistic = float(
        np.sum(
            weights
            * np.square(
                beta - fixed_beta
            )
        )
    )

    q_df = len(valid) - 1

    q_p = float(
        chi2.sf(
            q_statistic,
            q_df,
        )
    )

    i2_percent = (
        max(
            0.0,
            (
                q_statistic - q_df
            )
            / q_statistic,
        )
        * 100.0
        if q_statistic > 0
        else 0.0
    )

    # --------------------------------------------------------
    # DerSimonian-Laird random-effects meta-analysis
    # --------------------------------------------------------

    c_value = float(
        np.sum(weights)
        - (
            np.sum(
                np.square(weights)
            )
            / np.sum(weights)
        )
    )

    tau2 = (
        max(
            0.0,
            (
                q_statistic - q_df
            )
            / c_value,
        )
        if c_value > 0
        else 0.0
    )

    random_weights = (
        1.0
        / (
            np.square(se)
            + tau2
        )
    )

    random_beta = float(
        np.sum(
            random_weights * beta
        )
        / np.sum(random_weights)
    )

    random_se = float(
        np.sqrt(
            1.0
            / np.sum(
                random_weights
            )
        )
    )

    random_z = (
        random_beta / random_se
    )

    random_p = float(
        2.0
        * norm.sf(
            abs(random_z)
        )
    )

    # --------------------------------------------------------
    # Identify within-study ancestry estimates crossing GWS
    # --------------------------------------------------------

    gws_ancestries = (
        valid.loc[
            valid["p_value"]
            < GENOME_WIDE_THRESHOLD,
            "ancestry",
        ]
        .tolist()
    )

    if gws_ancestries:
        interpretation = (
            "Within-study ancestry-stratum "
            "genome-wide significance in "
            + ",".join(
                gws_ancestries
            )
            + "; independent same-ancestry "
            "replication required"
        )
    else:
        interpretation = (
            "No ancestry stratum reached the "
            "genome-wide threshold; within-study "
            "sensitivity evidence only"
        )

    meta_rows.append(
        {
            "variant_id": variant_id,
            "physical_group": target[
                "physical_group"
            ],
            "effect_allele": target[
                "effect_allele"
            ],
            "other_allele": target[
                "other_allele"
            ],
            "n_ancestries": len(valid),
            "ancestries_included": ",".join(
                valid["ancestry"]
            ),
            "fixed_beta": fixed_beta,
            "fixed_se": fixed_se,
            "fixed_or": math.exp(
                fixed_beta
            ),
            "fixed_ci_lower_or": math.exp(
                fixed_beta
                - 1.96 * fixed_se
            ),
            "fixed_ci_upper_or": math.exp(
                fixed_beta
                + 1.96 * fixed_se
            ),
            "fixed_p": fixed_p,
            "q_statistic": q_statistic,
            "q_df": q_df,
            "q_p": q_p,
            "i2_percent": i2_percent,
            "tau2": tau2,
            "random_beta": random_beta,
            "random_se": random_se,
            "random_or": math.exp(
                random_beta
            ),
            "random_ci_lower_or": math.exp(
                random_beta
                - 1.96 * random_se
            ),
            "random_ci_upper_or": math.exp(
                random_beta
                + 1.96 * random_se
            ),
            "random_p": random_p,
            "ancestries_crossing_genome_wide_threshold": (
                ",".join(
                    gws_ancestries
                )
            ),
            "within_study_interpretation": (
                interpretation
            ),
        }
    )

    # --------------------------------------------------------
    # Leave-one-ancestry-out fixed-effects analyses
    # --------------------------------------------------------

    if len(valid) >= 3:
        for ancestry_left_out in valid[
            "ancestry"
        ]:
            retained = valid.loc[
                valid["ancestry"].ne(
                    ancestry_left_out
                )
            ]

            retained_beta = (
                retained["beta"]
                .to_numpy(dtype=float)
            )

            retained_se = (
                retained["se"]
                .to_numpy(dtype=float)
            )

            retained_weights = (
                1.0
                / np.square(
                    retained_se
                )
            )

            pooled_beta = float(
                np.sum(
                    retained_weights
                    * retained_beta
                )
                / np.sum(
                    retained_weights
                )
            )

            pooled_se = float(
                np.sqrt(
                    1.0
                    / np.sum(
                        retained_weights
                    )
                )
            )

            pooled_p = float(
                2.0
                * norm.sf(
                    abs(
                        pooled_beta
                        / pooled_se
                    )
                )
            )

            leave_one_out_rows.append(
                {
                    "variant_id": variant_id,
                    "effect_allele": target[
                        "effect_allele"
                    ],
                    "ancestry_left_out": (
                        ancestry_left_out
                    ),
                    "ancestries_included": (
                        ",".join(
                            retained[
                                "ancestry"
                            ]
                        )
                    ),
                    "fixed_beta": pooled_beta,
                    "fixed_se": pooled_se,
                    "fixed_or": math.exp(
                        pooled_beta
                    ),
                    "fixed_ci_lower_or": (
                        math.exp(
                            pooled_beta
                            - 1.96
                            * pooled_se
                        )
                    ),
                    "fixed_ci_upper_or": (
                        math.exp(
                            pooled_beta
                            + 1.96
                            * pooled_se
                        )
                    ),
                    "fixed_p": pooled_p,
                }
            )


# ------------------------------------------------------------
# Assemble results
# ------------------------------------------------------------

ancestry_meta_results = pd.DataFrame(
    meta_rows
)

ancestry_leave_one_out_results = (
    pd.DataFrame(
        leave_one_out_rows
    )
)


# ------------------------------------------------------------
# Correct the heterogeneity tests
# ------------------------------------------------------------

valid_heterogeneity = (
    ancestry_meta_results[
        "q_p"
    ].notna()
)

n_heterogeneity_tests = int(
    valid_heterogeneity.sum()
)

ancestry_meta_results[
    "heterogeneity_bonferroni_adjusted_p"
] = np.nan

ancestry_meta_results[
    "heterogeneity_fdr_q_value"
] = np.nan

if n_heterogeneity_tests:
    ancestry_meta_results.loc[
        valid_heterogeneity,
        "heterogeneity_bonferroni_adjusted_p",
    ] = np.minimum(
        (
            ancestry_meta_results.loc[
                valid_heterogeneity,
                "q_p",
            ]
            * n_heterogeneity_tests
        ),
        1.0,
    )

    ancestry_meta_results.loc[
        valid_heterogeneity,
        "heterogeneity_fdr_q_value",
    ] = multipletests(
        ancestry_meta_results.loc[
            valid_heterogeneity,
            "q_p",
        ],
        method="fdr_bh",
    )[1]

ancestry_meta_results[
    "heterogeneity_after_bonferroni"
] = (
    ancestry_meta_results[
        "heterogeneity_bonferroni_adjusted_p"
    ]
    < 0.05
)

ancestry_meta_results[
    "heterogeneity_after_fdr"
] = (
    ancestry_meta_results[
        "heterogeneity_fdr_q_value"
    ]
    < 0.05
)


# ------------------------------------------------------------
# Extract QKI results
# ------------------------------------------------------------

qki_meta_check = (
    ancestry_meta_results.loc[
        ancestry_meta_results[
            "variant_id"
        ].eq(QKI_ID)
    ]
    .copy()
)

qki_leave_one_out_check = (
    ancestry_leave_one_out_results.loc[
        ancestry_leave_one_out_results[
            "variant_id"
        ].eq(QKI_ID)
    ]
    .copy()
)

if len(qki_meta_check) != 1:
    raise RuntimeError(
        "QKI meta-analysis row was not "
        "uniquely recovered."
    )

qki_gws_strata = (
    qki_stratified_check.loc[
        qki_stratified_check[
            "p_value"
        ]
        < GENOME_WIDE_THRESHOLD,
        "ancestry",
    ]
    .tolist()
)


# ------------------------------------------------------------
# Identify all within-study GWS ancestry results
# ------------------------------------------------------------

ancestry_gws_results = (
    ancestry_stratified_results.loc[
        ancestry_stratified_results[
            "crosses_genome_wide_threshold"
        ],
        [
            "variant_id",
            "physical_group",
            "ancestry",
            "effect_allele",
            "effect_allele_frequency",
            "odds_ratio",
            "ci_lower_or",
            "ci_upper_or",
            "p_value",
            "n_cases",
            "n_controls",
        ],
    ]
    .sort_values(
        [
            "p_value",
            "variant_id",
            "ancestry",
        ]
    )
    .reset_index(drop=True)
)


# ------------------------------------------------------------
# Identify corrected ancestry heterogeneity
# ------------------------------------------------------------

corrected_heterogeneity_results = (
    ancestry_meta_results.loc[
        ancestry_meta_results[
            "heterogeneity_after_fdr"
        ]
        | ancestry_meta_results[
            "heterogeneity_after_bonferroni"
        ],
        [
            "variant_id",
            "physical_group",
            "n_ancestries",
            "q_p",
            "i2_percent",
            "heterogeneity_bonferroni_adjusted_p",
            "heterogeneity_fdr_q_value",
            "heterogeneity_after_bonferroni",
            "heterogeneity_after_fdr",
        ],
    ]
    .sort_values(
        [
            "q_p",
            "variant_id",
        ]
    )
    .reset_index(drop=True)
)


# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

print("=" * 72)
print(
    "ANCESTRY META-ANALYSIS AND "
    "HETEROGENEITY AUDIT COMPLETE"
)
print("=" * 72)

print(
    "Variants with valid heterogeneity tests:",
    n_heterogeneity_tests,
)

print("\nQKI meta-analysis:")
display(
    qki_meta_check
)

print(
    "\nQKI leave-one-ancestry-out analyses:"
)
display(
    qki_leave_one_out_check
)

print(
    "\nQKI ancestry strata reaching "
    "P < 5e-8:",
    (
        qki_gws_strata
        if qki_gws_strata
        else "none"
    ),
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
    "\nVariants with heterogeneity surviving "
    "either Bonferroni or FDR correction:"
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
    "\nThese are within-study ancestry "
    "sensitivity analyses, not independent "
    "replication."
)