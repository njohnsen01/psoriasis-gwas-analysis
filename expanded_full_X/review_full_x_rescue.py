"""Protected final recovery review; only disclosure-limited aggregate files are staged.

Run after RESCUE_AND_LEADS_COMPUTED. This does not mark scientific review complete
or download anything. Original scan, checkpoints and unsuccessful fits are untouched.
"""
from pathlib import Path
import datetime
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

R = Path('/home/jupyter/x_followup_20261013/continuation_20261013')
O = R / 'full_x_rescue'
V = R / 'full_x_rescue_review'
FAMILY = 334926


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t')


def holm_full_family(p, family):
    """Unestimated tests stay NaN; family size includes all planned hypotheses."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    assert ok.sum() <= family and np.all((p[ok] >= 0) & (p[ok] <= 1))
    out = np.full(p.shape, np.nan)
    ix = np.flatnonzero(ok)
    order = np.argsort(p[ix], kind='stable')
    sorted_ix = ix[order]
    out[sorted_ix] = np.minimum(1., np.maximum.accumulate(
        p[sorted_ix] * (family - np.arange(len(sorted_ix)))))
    # Independent library comparison. Padding is solely multiplicity bookkeeping;
    # it never enters an output or converts an unestimated test into a null result.
    padded = np.r_[p[ok], np.ones(family - ok.sum())]
    reference = multipletests(padded, method='holm')[1][:ok.sum()]
    assert np.allclose(out[ok], reference, rtol=1e-12, atol=1e-15)
    return out


def verify_stats(frame, statistic, probability):
    assert np.isfinite(frame[statistic]).all() and (frame[statistic] >= 0).all()
    expected = chi2.sf(frame[statistic], 3)
    assert np.isfinite(frame[probability]).all()
    assert np.allclose(expected, frame[probability], rtol=1e-10, atol=1e-15)
    assert ((frame[probability] >= 0) & (frame[probability] <= 1)).all()


def verify_solver(frame, max_iterations):
    assert np.isfinite(frame[['max_score', 'max_decrement', 'max_condition']]).all().all()
    assert (frame.max_score <= 1e-6).all()
    assert (frame.max_decrement <= 1e-10).all()
    assert ((frame.max_condition > 0) & (frame.max_condition < 1e10)).all()
    assert ((frame.iterations >= 0) & (frame.iterations < max_iterations)).all()


def summary(frame, stat, p):
    return {
        'estimated': len(frame),
        'minimum_p': float(frame[p].min()) if len(frame) else None,
        'genome_wide_p_below_5e8': int((frame[p] < 5e-8).sum()),
        'full_family_Holm_below_005': int((frame.p_holm_full_family < .05).sum()),
        'full_family_Bonferroni_below_005': int((frame.p_bonferroni_full_family < .05).sum()),
        'threshold_counts': {str(t): int((frame[p] < t).sum())
                             for t in [.05, .01, .001, .0001, .00001, 5e-8]},
    }


def main():
    required = [
        O / 'RESCUE_AND_LEADS_COMPUTED.json',
        O / 'RECOVERY_COMPUTED.json',
        O / 'FIRTH_RECOVERY_COMPUTED.json',
        O / 'FINAL_LEAD_CHECKS_COMPUTED.json',
        O / 'BENCHMARK_VALIDATED.json',
        O / 'FIRTH_BENCHMARK_VALIDATED.json',
        O / 'R_ML_BENCHMARK_VALIDATED.json',
        R / 'full_x_estimability_audit/DIAGNOSTICS_COMPLETE.json',
        R / 'full_x_ploidy_audit/GLOBAL_PLOIDY_AUDIT.json',
    ]
    assert all(p.is_file() for p in required), 'Computation/validation not complete'
    receipts = {p.name: json.loads(p.read_text()) for p in required}
    for name in ['BENCHMARK_VALIDATED.json', 'FIRTH_BENCHMARK_VALIDATED.json',
                 'R_ML_BENCHMARK_VALIDATED.json', 'GLOBAL_PLOIDY_AUDIT.json']:
        assert receipts[name]['status'] == 'passed', name
    lead_receipt = receipts['FINAL_LEAD_CHECKS_COMPUTED.json']
    assert lead_receipt['status'] == 'lead_checks_passed_pending_final_review'
    assert lead_receipt['max_ML_stat_difference'] < 1e-4
    assert lead_receipt['max_Firth_stat_difference'] < 1e-4

    base = read(R / 'full_x/full_X_joint_Wald.tsv.gz')
    audit = read(R / 'full_x_estimability_audit/diagnostics.private.tsv')
    recovery = read(O / 'recovery.private.tsv')
    all_recovery = read(O / 'all_rescue_outcomes.private.tsv')
    firth = read(O / 'firth_recovery.private.tsv')
    for frame in [base, audit, recovery, all_recovery, firth]:
        assert frame.ID.is_unique
    assert len(base) == len(audit) == FAMILY and set(base.ID) == set(audit.ID)
    assert audit.full_genotype_rank.dtype == bool and audit.empty_outcome_cell.dtype == bool
    assert audit.full_genotype_rank.all(), 'Review rank-deficient tests separately'
    a = base.merge(audit, on='ID', validate='one_to_one')
    assert (a.ERRCODE == a.original_error).all()
    assert ((a.POS > 2781479) & (a.POS < 155701383)).all()
    assert len(all_recovery) == 21730
    assert set(all_recovery.ID) == set(a.loc[a.ERRCODE != '.', 'ID'])
    expected_ml = a.loc[(a.ERRCODE != '.') & ~a.empty_outcome_cell, 'ID']
    assert set(recovery.ID) == set(expected_ml)
    assert len(recovery) == receipts['RECOVERY_COMPUTED.json']['attempted']
    good_rec = recovery[recovery.status == 'estimated'].copy()
    assert set(good_rec.ID) == set(all_recovery.loc[all_recovery.status == 'estimated', 'ID'])
    verify_stats(good_rec, 'WALD_CHISQ_3DF', 'P_CHISQ_3DF')
    verify_solver(good_rec, 80)

    original = a[(a.ERRCODE == '.') & ~a.empty_outcome_cell].copy()
    verify_stats(original, 'WALD_CHISQ_3DF', 'P_CHISQ_3DF')
    original = original[['ID', 'WALD_CHISQ_3DF', 'P_CHISQ_3DF']]
    original['method'] = 'original_PLINK_joint_Wald'
    recovered = good_rec[['ID', 'WALD_CHISQ_3DF', 'P_CHISQ_3DF']].copy()
    recovered['method'] = 'equivalent_stable_ML_joint_Wald'
    wald = pd.concat([original, recovered], ignore_index=True)
    assert wald.ID.is_unique
    expected_firth = set(a.loc[a.full_genotype_rank, 'ID']) - set(wald.ID)
    assert set(firth.ID) == expected_firth, 'Incomplete or overlapping sensitivity scope'
    assert len(firth) == receipts['FIRTH_RECOVERY_COMPUTED.json']['attempted']
    good_firth = firth[firth.status == 'estimated'].copy()
    verify_stats(good_firth, 'FIRTH_CHISQ_3DF', 'P_FIRTH_3DF')
    verify_solver(good_firth, 500)
    assert len(wald) + len(good_firth) + int((firth.status != 'estimated').sum()) == FAMILY

    # Derive call counts and allele frequency only within Verily. Do not export
    # cells, carrier counts, frequencies or exact per-variant sample sizes.
    female = np.column_stack([a[f'cell_0_0_{k}'] + a[f'cell_0_1_{k}'] for k in range(3)])
    male = np.column_stack([a[f'cell_1_0_{k}'] + a[f'cell_1_1_{k}'] for k in range(3)])
    assert (male[:, 1] == 0).all()
    n = female.sum(axis=1) + male.sum(axis=1)
    assert ((n > 0) & (n <= 357530)).all()
    known_n = pd.to_numeric(a.OBS_CT, errors='coerce')
    assert np.array_equal(known_n[known_n.notna()], n[known_n.notna()])
    an = 2 * female.sum(axis=1) + male.sum(axis=1)
    ac = female[:, 1] + 2 * female[:, 2] + male[:, 2]
    maf = np.minimum(ac / an, 1 - ac / an)
    alt_carriers = female[:, 1] + female[:, 2] + male[:, 2]
    ref_carriers = female[:, 0] + female[:, 1] + male[:, 0]
    participant_cells = np.column_stack([alt_carriers, ref_carriers,
                                        n - alt_carriers, n - ref_carriers])
    no_small_derived_carrier_cell = ~((participant_cells >= 1) & (participant_cells <= 20)).any(axis=1)
    reportable = (maf >= .01) & no_small_derived_carrier_cell
    meta = a[['ID', 'POS']].copy()
    meta['n_complete_private'] = n
    meta['per_variant_export_allowed'] = reportable
    meta['N_rounded_to_100'] = (np.round(n / 100) * 100).astype(int)
    for frame in [good_rec, good_firth]:
        counts = frame[['ID', 'n_complete']].merge(meta, on='ID', validate='one_to_one')
        assert (counts.n_complete == counts.n_complete_private).all()

    # Separate estimators and keep the full planned family in both corrections.
    for frame, p in [(wald, 'P_CHISQ_3DF'), (good_firth, 'P_FIRTH_3DF')]:
        frame['p_holm_full_family'] = holm_full_family(frame[p], FAMILY)
        frame['p_bonferroni_full_family'] = np.minimum(1., frame[p] * FAMILY)
    wald = wald.merge(meta, on='ID', validate='one_to_one')
    good_firth = good_firth.merge(meta, on='ID', validate='one_to_one')

    # Independently recompute lead selection instead of trusting the saved marker.
    def spaced_ids(frame, p):
        positions, ids = [], []
        for row in frame.sort_values(p).itertuples():
            if all(abs(row.POS - pos) > 1000000 for pos in positions):
                positions.append(row.POS)
                ids.append(row.ID)
            if len(ids) == 10:
                break
        return set(ids)

    prior = {'chrX:5995921:G:A', 'chrX:5996494:GCT:G', 'chrX:5998724:A:G'}
    wald_selection = spaced_ids(wald, 'P_CHISQ_3DF') | set(wald.loc[wald.P_CHISQ_3DF < 5e-8, 'ID']) | prior
    firth_selection = spaced_ids(good_firth, 'P_FIRTH_3DF') | set(good_firth.loc[good_firth.P_FIRTH_3DF < 5e-8, 'ID'])
    leads = read(O / 'final_lead_checks/verified_leads.private.tsv')
    assert leads.ID.is_unique and set(leads.ID) == wald_selection | firth_selection
    assert (leads.status == 'estimated').all()
    assert leads.warnings.fillna('').eq('').all()
    assert set(leads.loc[leads.ml_eligible, 'ID']) == set(leads.ID) & set(wald.ID)
    assert len(leads) == lead_receipt['selected']
    leads = leads.merge(meta, on='ID', validate='one_to_one')
    leads['independent_R_Wald_P'] = chi2.sf(leads.wald, 3)
    verify_stats(leads, 'firth_stat', 'firth_p')
    for frame, p in [(wald, 'P_CHISQ_3DF'), (good_firth, 'P_FIRTH_3DF')]:
        joined = frame.merge(leads[['ID', 'wald', 'firth_stat', 'firth_p']], on='ID')
        ref_p = chi2.sf(joined.wald, 3) if p == 'P_CHISQ_3DF' else joined.firth_p
        # Every selected result must have the same genome-wide and full-family
        # Bonferroni decision under the independent reference implementation.
        for threshold in [5e-8, .05 / FAMILY]:
            assert np.array_equal(joined[p] < threshold, ref_p < threshold)
    compare = wald.merge(leads[['ID', 'wald']], on='ID')
    compare['stat_difference'] = np.abs(compare.WALD_CHISQ_3DF - compare.wald)
    assert (compare.loc[compare.method == 'equivalent_stable_ML_joint_Wald', 'stat_difference'] < 1e-4).all()
    # Original native approximations are reported, never silently overwritten.
    legacy_max_difference = float(compare.loc[compare.method == 'original_PLINK_joint_Wald', 'stat_difference'].max())

    V.mkdir(exist_ok=True)
    wfields = ['ID', 'POS', 'method', 'WALD_CHISQ_3DF', 'P_CHISQ_3DF',
               'p_holm_full_family', 'p_bonferroni_full_family', 'N_rounded_to_100']
    ffields = ['ID', 'POS', 'FIRTH_CHISQ_3DF', 'P_FIRTH_3DF',
               'p_holm_full_family', 'p_bonferroni_full_family', 'N_rounded_to_100']
    lfields = ['ID', 'POS', 'ml_eligible', 'selected_Wald', 'selected_Firth',
               'wald', 'independent_R_Wald_P', 'firth_stat', 'firth_p', 'N_rounded_to_100']
    for frame, fields, name in [(wald, wfields, 'full_X_eligible_Wald'),
                                (good_firth, ffields, 'full_X_Firth_sensitivity'),
                                (leads, lfields, 'full_X_independent_lead_checks')]:
        safe = frame.loc[frame.per_variant_export_allowed, fields].copy()
        safe.to_csv(V / (name + '.aggregate.tsv.gz'), sep='\t', index=False, compression='gzip')
    wald.loc[wald.per_variant_export_allowed, wfields].sort_values('P_CHISQ_3DF').head(100).to_csv(
        V / 'full_X_top100_Wald.aggregate.tsv', sep='\t', index=False)
    good_firth.loc[good_firth.per_variant_export_allowed, ffields].sort_values('P_FIRTH_3DF').head(100).to_csv(
        V / 'full_X_top100_Firth_sensitivity.aggregate.tsv', sep='\t', index=False)
    leads.loc[leads.ID.isin(prior) & leads.per_variant_export_allowed, lfields].to_csv(
        V / 'NLGN4X_full_X.aggregate.tsv', sep='\t', index=False)

    result = {
        'status': 'numerical_and_disclosure_checks_passed_pending_model_scientific_review',
        'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'planned_variant_family': FAMILY,
        'original_software_returned_results': int((a.ERRCODE == '.').sum()),
        'original_failed_results': len(all_recovery),
        'native_returned_but_separated': int(((a.ERRCODE == '.') & a.empty_outcome_cell).sum()),
        'original_nonseparated_results_retained': len(original),
        'ordinary_ML_recovered': len(good_rec),
        'ordinary_ML_recovery_status_counts': recovery.status.value_counts().to_dict(),
        'ordinary_ML_unestimated_after_review': FAMILY - len(wald),
        'Firth_sensitivity_attempted': len(firth),
        'Firth_sensitivity_status_counts': firth.status.value_counts().to_dict(),
        'unestimated_in_both_methods': int((firth.status != 'estimated').sum()),
        'eligible_Wald': summary(wald, 'WALD_CHISQ_3DF', 'P_CHISQ_3DF'),
        'Firth_sensitivity': summary(good_firth, 'FIRTH_CHISQ_3DF', 'P_FIRTH_3DF'),
        'eligible_Wald_median_chisq3_ratio_descriptive_only': float(wald.WALD_CHISQ_3DF.median() / chi2.ppf(.5, 3)),
        'independent_lead_variants_checked': len(leads),
        'independent_lead_R_vs_Python_max_ML_difference': lead_receipt['max_ML_stat_difference'],
        'independent_lead_R_vs_Python_max_Firth_difference': lead_receipt['max_Firth_stat_difference'],
        'selected_legacy_PLINK_vs_exact_R_max_Wald_difference': legacy_max_difference,
        'independent_selected_leads_Firth_below_5e8': int((leads.firth_p < 5e-8).sum()),
        'disclosure': {
            'individual_records_and_genotype_outcome_cells': 'remain protected; not exported',
            'per_variant_N': 'rounded to nearest 100; exact counts omitted',
            'record_eligibility': 'sex-aware pooled MAF >= 1%; no direct or derived carrier cells 1-20',
            'not_eligible_for_variant_record_export': int((~reportable).sum()),
            'Wald_records_staged': int(wald.per_variant_export_allowed.sum()),
            'Firth_records_staged': int(good_firth.per_variant_export_allowed.sum()),
            'lead_records_staged': int(leads.per_variant_export_allowed.sum()),
            'counts_in_this_receipt': 'numbers of tests or variants, not small participant cells',
        },
        'limitations': [
            'The 3-df ordinary logistic joint Wald scan has no LOCO/mixed-model adjustment; the separately completed rare-gene tests do use LOCO.',
            'The ordinary-ML retained subset and Firth sensitivity subset are selected using outcome support. Their distributions are not unbiased whole-scan null calibration evidence.',
            'The median statistic ratio and numerical implementation validation do not establish extreme-tail calibration.',
            'All 334926 planned hypotheses remain in multiplicity corrections; failed tests are inconclusive, not null results.',
            'Firth penalized likelihood-ratio P values are a separately labeled exploratory sensitivity and are not substituted for ordinary Wald P values.',
            'Independent lead computation uses the same participants and is verification, not independent replication.',
            'The rescue addenda were frozen before recovered association results but after original results; this is not original prospective prespecification.',
        ],
        'input_sha256': {str(p.relative_to(R)): sha(p) for p in required},
        'review_code_sha256': sha(Path(__file__)),
    }
    (V / 'NUMERICAL_DISCLOSURE_REVIEW.json').write_text(json.dumps(result, indent=2) + '\n')
    manifest = {p.name: sha(p) for p in sorted(V.glob('*.aggregate.*'))}
    (V / 'AGGREGATE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    # Only the aggregate receipt is printed. Native files, raw log text, rare
    # variant IDs and participant-level data are never sent to the terminal UI.
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
