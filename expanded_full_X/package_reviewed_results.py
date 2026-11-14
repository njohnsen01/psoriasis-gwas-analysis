"""Package an explicit allowlist of reviewed aggregate outputs only."""
from pathlib import Path
import datetime, hashlib, json, shutil, zipfile
import numpy as np
import pandas as pd

R = Path('/home/jupyter/x_followup_20261013/continuation_20261013')
V = R / 'full_x_rescue_review'
PACKAGE = R / 'reviewed_results_package_20261014'
DEST = Path('/home/jupyter/PSORIASIS_X_REVIEWED_RESULTS_20261014.zip')


def main():
    review = json.loads((V / 'SCIENTIFIC_REVIEW.json').read_text())
    assert review['status'] == 'full_X_results_scientifically_reviewed'
    rare = R / 'rare_gene_tests'
    assert json.loads((rare / 'SCIENTIFIC_REVIEW.json').read_text())['status'] == 'rare_gene_results_scientifically_reviewed'
    assert not DEST.exists(), 'Existing package must be preserved and reviewed before replacement'
    PACKAGE.mkdir(exist_ok=True)
    assert not list(PACKAGE.iterdir()), 'Package directory is not empty'
    names = [
        'full_X_eligible_Wald.aggregate.tsv.gz',
        'full_X_Firth_sensitivity.aggregate.tsv.gz',
        'full_X_independent_lead_checks.aggregate.tsv.gz',
        'full_X_top100_Wald.aggregate.tsv',
        'full_X_top100_Firth_sensitivity.aggregate.tsv',
        'NLGN4X_full_X.aggregate.tsv',
        'full_X_recovery_diagnostics.png',
        'full_X_recovery_diagnostics.pdf',
        'full_X_final_scientific_review_20261014.md',
        'NUMERICAL_DISCLOSURE_REVIEW.json',
        'SCIENTIFIC_REVIEW.json',
    ]
    for name in names:
        assert (V / name).is_file()
        shutil.copy2(V / name, PACKAGE / name)
    rare_names = ['primary_gene_results.aggregate.tsv', 'secondary_gene_results.aggregate.tsv',
                  'top20_primary.aggregate.tsv', 'NLGN4X.aggregate.tsv']
    for name in rare_names:
        df = pd.read_csv(rare / name, sep='\t')
        prohibited = {'FID', 'IID', 'person_id', 'sample_id', 'participant_id', 'variant_id',
                      'pooled_carriers_private', 'carriers', 'carrier_count', 'cell_count'}
        assert not prohibited.intersection(df.columns)
        assert not any(c.startswith('cell_') or c.endswith('_private') for c in df.columns)
        for col in ['combined_p', 'combined_minus_log10_p', 'p_holm', 'p_bonferroni',
                    'p_holm_secondary', 'p_bonferroni_secondary']:
            if col in df:
                assert df.loc[~df.eligible_for_external_report, col].isna().all()
                assert df.loc[~df.both_sex_tests_estimated, col].isna().all()
        for sex in ['female', 'male']:
            if sex + '_result_reportable' in df:
                assert df.loc[~df[sex + '_result_reportable'], sex + '_minus_log10_p'].isna().all()
        shutil.copy2(rare / name, PACKAGE / ('rare_' + name))
    primary = pd.read_csv(rare / 'primary_gene_results.aggregate.tsv', sep='\t')
    secondary = pd.read_csv(rare / 'secondary_gene_results.aggregate.tsv', sep='\t')
    assert len(primary) == 1294 and primary.both_sex_tests_estimated.sum() == 974
    assert len(secondary) == 3882 and secondary.both_sex_tests_estimated.sum() == 2766
    assert not (primary.p_holm < .05).any() and not (secondary.p_holm_secondary < .05).any()
    nr = primary[(primary.gene_symbol == 'NLGN4X') & (primary['mask'] == 'PTV_DAMAGING')]
    assert len(nr) == 1 and nr.eligible_for_external_report.all()
    assert np.isclose(float(nr.combined_p.iloc[0]), .80883731484, rtol=1e-10)
    rr = {'status': 'reviewed_rare_gene_aggregate_export',
          'genes': 647, 'qualifying_variants': 10931,
          'planned_primary_tests': 1294, 'estimated_both_sexes_primary': 974,
          'unestimated_primary': 320, 'externally_reportable_estimated_primary': 577,
          'planned_secondary_tests': 3882, 'estimated_both_sexes_secondary': 2766,
          'primary_Holm_significant': 0, 'secondary_Holm_significant': 0,
          'NLGN4X_PTV_DAMAGING_p': float(nr.combined_p.iloc[0]),
          'NLGN4X_PTV_DAMAGING_Holm': float(nr.p_holm.iloc[0]),
          'interpretation': 'No corrected association detected; unestimated or suppressed tests are inconclusive. Independent sex-specific LOCO-adjusted tests were combined only when both probabilities were estimated.',
          'disclosure': 'No rare-variant records, carrier counts, participant identifiers, or raw logs are included. Gene-mask probabilities are suppressed unless eligible; sex-specific probabilities require their own carrier checks.',
          'rare_scientific_review_sha256': hashlib.sha256((rare / 'SCIENTIFIC_REVIEW.json').read_bytes()).hexdigest()}
    (PACKAGE / 'rare_review.aggregate.json').write_text(json.dumps(rr, indent=2) + '\n')
    (PACKAGE / 'README.txt').write_text(
        'Reviewed psoriasis chromosome-X follow-up results, 14 October 2026.\n\n'
        'Start with full_X_final_scientific_review_20261014.md.\n'
        'The ordinary Wald and Firth sensitivity tables are separate analyses.\n'
        'Full-X table P values use chi-square with 3 degrees of freedom.\n'
        'Holm/Bonferroni denominators retain all 334926 planned variants.\n'
        'N_rounded_to_100 is intentionally coarsened and is not an exact denominator.\n'
        'Exported variant rows are a disclosure-qualified subset of those analyzed.\n'
        'Rare tables contain gene-mask aggregates only. Blank probabilities are not null results.\n'
        'Rare primary correction covers 1294 gene-mask tests; secondary covers 3882.\n'
        'Lead checks use the same cohort and do not constitute replication.\n'
        'TSV files are tab-delimited; .gz files are gzip-compressed TSV.\n'
        'The numerical receipt precedes the separate final scientific-review receipt.\n'
        'Original inputs, private data, raw logs and unsuccessful validation attempts remain in Verily.\n'
        'No manuscript submission occurred.\n')
    manifest = {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                for p in sorted(PACKAGE.iterdir()) if p.is_file()}
    (PACKAGE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with zipfile.ZipFile(DEST, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(PACKAGE.iterdir()):
            assert p.is_file() and '.private.' not in p.name and not p.name.endswith('.log')
            z.write(p, p.name)
    with zipfile.ZipFile(DEST) as z:
        assert z.testzip() is None
    receipt = {'status': 'reviewed_aggregate_package_created',
               'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'file': DEST.name, 'files': len(manifest) + 1, 'bytes': DEST.stat().st_size,
               'sha256': hashlib.sha256(DEST.read_bytes()).hexdigest()}
    (V / 'PACKAGE_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
