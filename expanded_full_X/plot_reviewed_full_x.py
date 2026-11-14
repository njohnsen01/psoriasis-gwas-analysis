"""Plot only the disclosure-qualified aggregate recovery outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

V = Path('/home/jupyter/x_followup_20261013/continuation_20261013/full_x_rescue_review')


def main():
    review = json.loads((V / 'NUMERICAL_DISCLOSURE_REVIEW.json').read_text())
    assert review['status'] == 'numerical_and_disclosure_checks_passed_pending_model_scientific_review'
    w = pd.read_csv(V / 'full_X_eligible_Wald.aggregate.tsv.gz', sep='\t')
    f = pd.read_csv(V / 'full_X_Firth_sensitivity.aggregate.tsv.gz', sep='\t')
    assert len(w) == review['disclosure']['Wald_records_staged']
    assert len(f) == review['disclosure']['Firth_records_staged']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={'width_ratios': [1, 1.2, 1.2]})
    floor = np.nextafter(0., 1.)
    p = np.sort(w.P_CHISQ_3DF.to_numpy())
    observed = -np.log10(np.maximum(p, floor))
    expected = -np.log10((np.arange(len(p)) + .5) / len(p))
    ix = np.unique(np.r_[np.arange(min(len(p), 1500)), np.linspace(0, len(p)-1, min(len(p), 5000)).astype(int)])
    axes[0].scatter(expected[ix], observed[ix], s=6, alpha=.65, color='#265b82', rasterized=True)
    limit = max(float(expected.max()), float(observed.max()))
    axes[0].plot([0, limit], [0, limit], color='#888888', lw=1)
    axes[0].set(xlabel='Expected −log₁₀(P)', ylabel='Observed −log₁₀(P)', title='Eligible ordinary Wald tests\nDescriptive Q–Q plot')
    for ax, data, col, title, color in [
        (axes[1], w, 'P_CHISQ_3DF', 'Ordinary logistic joint Wald', '#265b82'),
        (axes[2], f, 'P_FIRTH_3DF', 'Separate Firth sensitivity', '#a86627'),
    ]:
        y = -np.log10(np.maximum(data[col].to_numpy(), floor))
        ax.scatter(data.POS / 1e6, y, s=3, alpha=.5, color=color, rasterized=True)
        ax.axhline(-np.log10(5e-8), ls='--', lw=1, color='#af3c3c')
        ax.set(xlabel='Chromosome X position (Mb; GRCh38)', ylabel='−log₁₀(P)',
               title=f'{title}\n{len(data):,} disclosure-qualified variants',
               xlim=(2.78, 155.71), ylim=(0, max(8., float(y.max()) + .4)))
    fig.text(.025, .055, 'Dashed line: P = 5 × 10⁻⁸. The estimable and sensitivity subsets depend on outcome support; these plots do not establish null calibration.', fontsize=8.5)
    fig.text(.025, .022, 'Only reviewed common-variant records are plotted. Independent lead verification uses the same cohort and is not replication.', fontsize=8.5)
    fig.tight_layout(rect=[0, .1, 1, 1])
    fig.savefig(V / 'full_X_recovery_diagnostics.png', dpi=180)
    fig.savefig(V / 'full_X_recovery_diagnostics.pdf')
    plt.close(fig)
    print('AGGREGATE_DIAGNOSTIC_PLOTS_CREATED')


if __name__ == '__main__':
    main()
