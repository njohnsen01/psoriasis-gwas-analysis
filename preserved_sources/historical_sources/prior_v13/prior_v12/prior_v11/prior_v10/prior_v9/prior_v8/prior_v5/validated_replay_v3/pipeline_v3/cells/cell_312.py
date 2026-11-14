# ============================================================
# Corrected 25-lead dosage harmonization and biomarker analyses
# Primary GWAS is NOT rerun by this block.
# ============================================================

import os
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
from scipy import stats

ATTEMPT = f"{os.environ['PSO_PRIMARY_RUN_ID']}"
GWAS_ROOT = Path(
    "/home/jupyter/psoriasis_regenie_loco_step2"
) / f"attempt_{ATTEMPT}" / "final_aggregate"
DOWNSTREAM_ROOT = Path(
    "/home/jupyter/psoriasis_corrected_downstream_2026-07-31"
)
PRIVATE = DOWNSTREAM_ROOT / "private_participant_level"
AGG = DOWNSTREAM_ROOT / "aggregate_outputs"
PRIVATE.mkdir(parents=True, exist_ok=True)
AGG.mkdir(parents=True, exist_ok=True)

GROUPS_PATH = GWAS_ROOT / "psoriasis_corrected_lead_centered_100kb_groups.tsv"
BIOMARKER_BASE_PATH = Path(
    f"{os.environ['PSO_WORKSPACE_MOUNT'].rstrip('/')}/"
    "psoriasis_biomarker_analysis_base_with_lead_dosages.tsv"
)
NEW8_PATH = (
    PRIVATE / "corrected_8_missing_lead_dosages_effect_allele.tsv.gz"
)

for required_path in [GROUPS_PATH, BIOMARKER_BASE_PATH, NEW8_PATH]:
    if not required_path.exists() or required_path.stat().st_size == 0:
        raise FileNotFoundError(f"Required input missing: {required_path}")

print("Loading corrected lead groups...", flush=True)
groups = pd.read_csv(GROUPS_PATH, sep="\t")
groups.columns = [str(c).lstrip("#") for c in groups.columns]
required_group_columns = {"ID", "REF", "ALT"}
missing_group_columns = required_group_columns.difference(groups.columns)
if missing_group_columns:
    raise RuntimeError(
        f"Corrected group file lacks columns: {sorted(missing_group_columns)}"
    )
if len(groups) != 25 or groups["ID"].nunique() != 25:
    raise RuntimeError(
        f"Expected 25 unique corrected leads; found {len(groups)} rows "
        f"and {groups['ID'].nunique()} unique IDs."
    )

print("Loading frozen biomarker analysis base...", flush=True)
base = pd.read_csv(BIOMARKER_BASE_PATH, sep="\t", low_memory=False)
print(f"Biomarker base: {len(base):,} rows, {len(base.columns):,} columns", flush=True)

print("Loading eight newly extracted effect-allele dosages...", flush=True)
new8 = pd.read_csv(NEW8_PATH, sep="\t", compression="gzip", low_memory=False)
print(f"New dosage matrix: {len(new8):,} rows, {len(new8.columns):,} columns", flush=True)

def choose_id_column(frame, label):
    preferred = [
        "research_id", "person_id", "IID", "iid",
        "participant_id", "sample_id"
    ]
    for candidate in preferred:
        if candidate in frame.columns:
            return candidate
    for candidate in frame.columns:
        low = str(candidate).lower()
        if low in {"research_id", "person_id", "iid", "participant_id"}:
            return candidate
    raise RuntimeError(f"Could not identify participant ID column in {label}.")

base_id_col = choose_id_column(base, "biomarker base")
new8_id_col = choose_id_column(new8, "new eight-dosage file")
base_ids = base[base_id_col].astype(str)
new8_ids = new8[new8_id_col].astype(str)

if base_ids.duplicated().any():
    raise RuntimeError("Biomarker base participant IDs are not unique.")
if new8_ids.duplicated().any():
    raise RuntimeError("New dosage participant IDs are not unique.")
if len(base) != 392155:
    raise RuntimeError(f"Expected 392,155 biomarker-base rows; found {len(base):,}.")
if len(new8) != 392155:
    raise RuntimeError(f"Expected 392,155 new-dosage rows; found {len(new8):,}.")
if set(base_ids) != set(new8_ids):
    raise RuntimeError("Participant sets differ between biomarker base and new dosages.")

new8_ordered = new8.assign(_join_id=new8_ids).set_index("_join_id").loc[base_ids]

def dosage_candidates(frame, variant_id):
    exact = []
    prefixed = []
    for column in frame.columns:
        name = str(column)
        if name == variant_id:
            exact.append(column)
        elif name.startswith(variant_id + "_"):
            prefixed.append(column)
    return exact + prefixed

def allele_from_column(column_name, variant_id):
    name = str(column_name)
    if name == variant_id:
        return None
    suffix = name[len(variant_id) + 1:]
    if suffix and ":" not in suffix and "_" not in suffix:
        return suffix
    return None

new8_variant_ids = set()
for variant_id in groups["ID"].astype(str):
    if dosage_candidates(new8, variant_id):
        new8_variant_ids.add(variant_id)

if len(new8_variant_ids) != 8:
    raise RuntimeError(
        "Expected eight corrected leads in the new dosage file; found "
        f"{len(new8_variant_ids)}: {sorted(new8_variant_ids)}"
    )

aligned = pd.DataFrame({base_id_col: base_ids.to_numpy()})
qc_rows = []
variant_to_aligned_column = {}

for lead_index, row in groups.reset_index(drop=True).iterrows():
    variant_id = str(row["ID"])
    ref = str(row["REF"])
    alt = str(row["ALT"])
    source_name = None
    source_frame = None
    source_series = None

    new_candidates = dosage_candidates(new8, variant_id)
    base_candidates = dosage_candidates(base, variant_id)

    if variant_id in new8_variant_ids:
        if len(new_candidates) != 1:
            raise RuntimeError(
                f"Expected one new dosage column for {variant_id}; "
                f"found {new_candidates}"
            )
        source_name = "new_targeted_extraction"
        source_frame = "new8"
        source_column = new_candidates[0]
        source_series = new8_ordered[source_column]
    else:
        if len(base_candidates) != 1:
            raise RuntimeError(
                f"Expected one existing biomarker-base dosage column for "
                f"{variant_id}; found {base_candidates}"
            )
        source_name = "existing_biomarker_base"
        source_frame = "base"
        source_column = base_candidates[0]
        source_series = base[source_column]

    raw = pd.to_numeric(source_series, errors="coerce").to_numpy(dtype=float)
    source_allele = allele_from_column(source_column, variant_id)

    if source_allele == alt:
        effect_dosage = raw
        action = "retained"
    elif source_allele == ref:
        effect_dosage = 2.0 - raw
        action = "flipped_2_minus_dosage"
    elif source_name == "new_targeted_extraction":
        effect_dosage = raw
        action = "retained_manifest_effect_aligned"
        source_allele = alt
    else:
        raise RuntimeError(
            f"Could not prove allele alignment for {variant_id} from "
            f"column {source_column!r}; REF={ref}, ALT={alt}."
        )

    finite = np.isfinite(effect_dosage)
    if finite.sum() == 0:
        raise RuntimeError(f"No finite dosages for {variant_id}.")
    observed_min = float(np.nanmin(effect_dosage))
    observed_max = float(np.nanmax(effect_dosage))
    if observed_min < -1e-6 or observed_max > 2.000001:
        raise RuntimeError(
            f"Dosage outside [0,2] for {variant_id}: "
            f"{observed_min}, {observed_max}"
        )

    aligned_column = f"effect_dosage_{lead_index + 1:02d}"
    aligned[aligned_column] = effect_dosage
    variant_to_aligned_column[variant_id] = aligned_column

    gwas_a1freq = (
        float(row["A1FREQ"])
        if "A1FREQ" in row.index and pd.notna(row["A1FREQ"])
        else np.nan
    )
    dosage_af = float(np.nanmean(effect_dosage) / 2.0)
    af_difference = (
        dosage_af - gwas_a1freq
        if np.isfinite(gwas_a1freq)
        else np.nan
    )
    if np.isfinite(af_difference) and abs(af_difference) > 0.01:
        raise RuntimeError(
            f"Effect-allele frequency validation failed for {variant_id}: "
            f"dosage AF={dosage_af:.6f}, GWAS A1FREQ={gwas_a1freq:.6f}"
        )

    qc_rows.append({
        "physical_group": row.get("physical_group", f"G{lead_index + 1:02d}"),
        "variant_id": variant_id,
        "ref": ref,
        "effect_allele": alt,
        "source": source_name,
        "source_column": str(source_column),
        "source_counted_allele": source_allele,
        "alignment_action": action,
        "n_total": int(len(effect_dosage)),
        "n_nonmissing": int(finite.sum()),
        "n_missing": int((~finite).sum()),
        "missing_rate": float((~finite).mean()),
        "dosage_mean": float(np.nanmean(effect_dosage)),
        "dosage_sd": float(np.nanstd(effect_dosage, ddof=1)),
        "dosage_min": observed_min,
        "dosage_max": observed_max,
        "dosage_effect_allele_frequency": dosage_af,
        "corrected_gwas_a1freq": gwas_a1freq,
        "frequency_difference": af_difference,
        "aligned_column": aligned_column,
    })

qc = pd.DataFrame(qc_rows)
qc_path = AGG / "corrected_25_lead_dosage_qc.tsv"
qc.to_csv(qc_path, sep="\t", index=False)

dosage_path = PRIVATE / "corrected_25_lead_dosages_effect_allele.tsv.gz"
dosage_export = aligned.rename(columns={base_id_col: "research_id"})
dosage_export.to_csv(
    dosage_path, sep="\t", index=False, compression="gzip"
)

print("All 25 lead dosages aligned and validated.", flush=True)
print(
    qc[[
        "variant_id", "effect_allele", "source",
        "alignment_action", "missing_rate", "frequency_difference"
    ]].to_string(index=False),
    flush=True,
)

# ----------------------------
# Exact OLS inference via FWL
# ----------------------------
trait_specs = [
    ("Systemic Immune-Inflammation Index", "log_inflammatory_index_z"),
    ("HDL cholesterol", "HDL_z"),
    ("LDL cholesterol", "LDL_z"),
    ("Triglycerides", "log_TG_z"),
]
base_covariates = [
    "age", "sex",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PC6", "PC7", "PC8", "PC9", "PC10",
]
full_covariates = base_covariates + ["psoriasis_case"]

required_model_columns = {
    outcome for _, outcome in trait_specs
}.union(full_covariates)
missing_model_columns = required_model_columns.difference(base.columns)
if missing_model_columns:
    raise RuntimeError(
        f"Biomarker base lacks model columns: {sorted(missing_model_columns)}"
    )

def ols_variant(y_series, g_array, covariate_frame, analysis_mask):
    y = pd.to_numeric(y_series, errors="coerce").to_numpy(dtype=float)
    cov = covariate_frame.apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    valid = (
        analysis_mask
        & np.isfinite(y)
        & np.isfinite(g_array)
        & np.isfinite(cov).all(axis=1)
    )
    n = int(valid.sum())
    x_cov = np.column_stack([
        np.ones(n, dtype=float),
        cov[valid, :],
    ])
    yv = y[valid]
    gv = g_array[valid]
    cov_rank = int(np.linalg.matrix_rank(x_cov))
    if n <= cov_rank + 2:
        raise RuntimeError("Insufficient complete observations for OLS.")
    coef_y, _, _, _ = np.linalg.lstsq(x_cov, yv, rcond=None)
    coef_g, _, _, _ = np.linalg.lstsq(x_cov, gv, rcond=None)
    y_resid = yv - x_cov @ coef_y
    g_resid = gv - x_cov @ coef_g
    denominator = float(g_resid @ g_resid)
    if denominator <= 0:
        raise RuntimeError("Residualized genotype variance is zero.")
    beta = float((g_resid @ y_resid) / denominator)
    residual = y_resid - beta * g_resid
    df_resid = int(n - cov_rank - 1)
    sigma2 = float((residual @ residual) / df_resid)
    se = float(math.sqrt(sigma2 / denominator))
    t_stat = float(beta / se)
    p_value = float(2.0 * stats.t.sf(abs(t_stat), df_resid))
    critical = float(stats.t.ppf(0.975, df_resid))
    return {
        "n": n,
        "beta": beta,
        "se": se,
        "ci_lower": beta - critical * se,
        "ci_upper": beta + critical * se,
        "t_stat": t_stat,
        "df_resid": df_resid,
        "p_value": p_value,
    }

def add_multiplicity(results):
    output = results.copy()
    p = output["p_value"].to_numpy(dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m, dtype=float)
    q[order] = np.clip(adjusted_sorted, 0.0, 1.0)
    output["fdr_q_value"] = q
    output["bonferroni_adjusted_p"] = np.minimum(p * m, 1.0)
    output["bonferroni_threshold"] = 0.05 / m
    output["bonferroni_significant"] = p < (0.05 / m)
    output["fdr_significant"] = q < 0.05
    output["significance_tier"] = np.where(
        output["bonferroni_significant"],
        "Bonferroni significant",
        np.where(
            output["fdr_significant"],
            "FDR significant only",
            "Not significant after correction",
        ),
    )
    return output

def run_analysis(analysis_name, sample_mask, covariate_columns):
    records = []
    covariate_frame = base[covariate_columns]
    for trait_label, outcome_column in trait_specs:
        print(
            f"Running {analysis_name}: {trait_label} (25 variants)",
            flush=True,
        )
        for _, lead in groups.iterrows():
            variant_id = str(lead["ID"])
            aligned_column = variant_to_aligned_column[variant_id]
            model = ols_variant(
                base[outcome_column],
                aligned[aligned_column].to_numpy(dtype=float),
                covariate_frame,
                sample_mask,
            )
            records.append({
                "analysis": analysis_name,
                "trait": trait_label,
                "outcome_column": outcome_column,
                "physical_group": lead.get("physical_group"),
                "chromosome": lead.get("CHROM"),
                "position": lead.get("POS", lead.get("GENPOS")),
                "variant_id": variant_id,
                "reference_allele": str(lead["REF"]),
                "effect_allele": str(lead["ALT"]),
                "corrected_gwas_beta": lead.get("BETA"),
                "corrected_gwas_or": lead.get("OR"),
                "corrected_gwas_p": lead.get("P"),
                "covariates": ",".join(covariate_columns),
                **model,
            })
    return add_multiplicity(pd.DataFrame(records))

all_samples = np.ones(len(base), dtype=bool)
case_status = pd.to_numeric(
    base["psoriasis_case"], errors="coerce"
).to_numpy(dtype=float)
cases_only = case_status == 1

full_results = run_analysis(
    "full_cohort_adjusted_for_psoriasis_status",
    all_samples,
    full_covariates,
)
cases_results = run_analysis(
    "psoriasis_cases_only",
    cases_only,
    base_covariates,
)

full_path = AGG / "corrected_25_lead_biomarker_results_full_cohort.tsv"
cases_path = AGG / "corrected_25_lead_biomarker_results_cases_only.tsv"
full_fdr_path = (
    AGG / "corrected_25_lead_biomarker_fdr_significant_full_cohort.tsv"
)
cases_fdr_path = (
    AGG / "corrected_25_lead_biomarker_fdr_significant_cases_only.tsv"
)

full_results.to_csv(full_path, sep="\t", index=False)
cases_results.to_csv(cases_path, sep="\t", index=False)
full_results.loc[full_results["fdr_significant"]].to_csv(
    full_fdr_path, sep="\t", index=False
)
cases_results.loc[cases_results["fdr_significant"]].to_csv(
    cases_fdr_path, sep="\t", index=False
)

summary = {
    "analysis_date": "2026-07-31",
    "corrected_lead_variants": 25,
    "traits": 4,
    "tests_per_analysis": 100,
    "bonferroni_threshold": 0.0005,
    "full_cohort_n_total": int(len(base)),
    "psoriasis_cases_n_total": int(cases_only.sum()),
    "full_cohort_fdr_significant": int(full_results["fdr_significant"].sum()),
    "full_cohort_bonferroni_significant": int(
        full_results["bonferroni_significant"].sum()
    ),
    "cases_only_fdr_significant": int(cases_results["fdr_significant"].sum()),
    "cases_only_bonferroni_significant": int(
        cases_results["bonferroni_significant"].sum()
    ),
    "full_cohort_covariates": full_covariates,
    "cases_only_covariates": base_covariates,
    "outcomes": {label: column for label, column in trait_specs},
    "transformation_note": (
        "HDL and LDL were standardized. Triglycerides and the Systemic "
        "Immune-Inflammation Index were natural-log transformed before "
        "standardization, matching the frozen prior biomarker pipeline."
    ),
    "dosage_alignment": (
        "All 25 dosages count the corrected primary-GWAS effect allele "
        "(REGENIE ALLELE1/ALT)."
    ),
    "private_participant_level_dosage_file": str(dosage_path),
}
summary_json_path = (
    AGG / "corrected_25_lead_biomarker_analysis_summary.json"
)
summary_tsv_path = (
    AGG / "corrected_25_lead_biomarker_analysis_summary.tsv"
)
summary_json_path.write_text(json.dumps(summary, indent=2) + "\n")
pd.DataFrame([{
    key: json.dumps(value) if isinstance(value, (list, dict)) else value
    for key, value in summary.items()
}]).to_csv(summary_tsv_path, sep="\t", index=False)

print("\n" + "=" * 72, flush=True)
print("CORRECTED 25-LEAD BIOMARKER ANALYSIS COMPLETE", flush=True)
print("=" * 72, flush=True)
print(json.dumps(summary, indent=2), flush=True)
print("\nFull-cohort corrected-significant associations:", flush=True)
display_columns = [
    "trait", "physical_group", "variant_id", "effect_allele",
    "beta", "ci_lower", "ci_upper", "p_value",
    "fdr_q_value", "bonferroni_adjusted_p", "significance_tier",
]
print(
    full_results.loc[
        full_results["fdr_significant"], display_columns
    ].sort_values(
        ["trait", "bonferroni_significant", "p_value"],
        ascending=[True, False, True],
    ).to_string(index=False),
    flush=True,
)
print("\nCases-only corrected-significant associations:", flush=True)
print(
    cases_results.loc[
        cases_results["fdr_significant"], display_columns
    ].sort_values(
        ["trait", "bonferroni_significant", "p_value"],
        ascending=[True, False, True],
    ).to_string(index=False),
    flush=True,
)
print("\nAggregate outputs:", flush=True)
for output_path in [
    qc_path, full_path, cases_path, full_fdr_path, cases_fdr_path,
    summary_json_path, summary_tsv_path,
]:
    print(output_path, flush=True)
print("\nDo not download the participant-level dosage file:", flush=True)
print(dosage_path, flush=True)
