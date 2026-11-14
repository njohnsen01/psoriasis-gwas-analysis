# Refit the genome-wide significant joint-scan variants under two nuisance
# specifications, holding participants and genotypes fixed.
#   Model A: shared ancestry (PC) slopes across sexes   <- the specification that produced the scan
#   Model C: Model A + sex x PC1..PC10
# Test in both: the 2 df joint null gf = gm = 0.
# Ordinary maximum likelihood is used for BOTH models so that the ONLY thing
# that differs between them is the nuisance specification.
import os, glob, json
import numpy as np, pandas as pd
WS = chr(92) + 's+'   # whitespace regex, built without a literal backslash
from scipy.stats import chi2

N = os.environ['N']
OUT = os.path.expanduser('~/joint_refit_20261014')

# ---- genotypes -------------------------------------------------------------
raws = sorted(glob.glob(os.path.join(OUT, 'geno_*.raw')))
assert raws, 'no .raw files'
geno = None
for r in raws:
    d = pd.read_csv(r, sep=WS)
    keep = [c for c in d.columns if c not in
            ('FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE')]
    d = d[keep]
    d['IID'] = d['IID'].astype(str)
    geno = d if geno is None else geno.merge(d, on='IID', how='outer', validate='one_to_one')
vcols = [c for c in geno.columns if c != 'IID']
print('raw files', len(raws), 'variant columns', len(vcols), 'rows', len(geno))

# ---- covariates + phenotype ------------------------------------------------
cov = pd.read_csv(os.path.join(N, 'inputs/covariates.tsv'), sep=WS)
cov.columns = [c.lstrip('#') for c in cov.columns]
phe = pd.read_csv(os.path.join(N, 'inputs/phenotype.tsv'), sep=WS)
phe.columns = [c.lstrip('#') for c in phe.columns]
cov['IID'] = cov['IID'].astype(str); phe['IID'] = phe['IID'].astype(str)
d = cov.merge(phe[['IID', 'psoriasis_case']], on='IID', validate='one_to_one')
d = d.merge(geno, on='IID', how='inner', validate='one_to_one')
print('merged rows', len(d))

PCS = [f'PC{i}' for i in range(1, 11)]
base = ['sex', 'age', 'age2', 'age_sex'] + PCS
for k in range(1, 11):
    d[f'sex_PC{k}'] = d['sex'] * d[f'PC{k}']
sexpc = [f'sex_PC{k}' for k in range(1, 11)]

def irls(X, y, maxit=100, tol=1e-10):
    n, p = X.shape
    b = np.zeros(p)
    ll_old = -np.inf
    for it in range(maxit):
        eta = np.clip(X @ b, -35, 35)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(mu * (1 - mu), 1e-10)
        z = eta + (y - mu) / w
        XTW = X.T * w
        try:
            b_new = np.linalg.solve(XTW @ X, XTW @ z)
        except np.linalg.LinAlgError:
            b_new = np.linalg.lstsq(XTW @ X, XTW @ z, rcond=None)[0]
        eta = np.clip(X @ b_new, -35, 35)
        ll = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
        if not np.isfinite(ll):
            return b, ll_old, False
        if abs(ll - ll_old) < tol * (abs(ll) + 1):
            return b_new, ll, True
        b, ll_old = b_new, ll
    return b, ll_old, False

rows = []
for v in vcols:
    g = d[v].to_numpy(float)
    ok = np.isfinite(g)
    z = d.loc[ok].copy()
    gg = g[ok]
    z['gf'] = gg * (1 - z['sex'].to_numpy(float))
    z['gm'] = gg * z['sex'].to_numpy(float) / 2.0
    y = z['psoriasis_case'].to_numpy(float)
    rec = {'variant': v, 'n': int(ok.sum()),
           'n_cases': int(y.sum()), 'maf': float(np.nanmean(gg) / 2.0)}
    for label, extra in (('A_shared', []), ('C_sexPC', sexpc)):
        cols = base + extra
        Xn = np.column_stack([np.ones(len(z))] + [z[c].to_numpy(float) for c in cols])
        Xf = np.column_stack([Xn, z['gf'].to_numpy(float), z['gm'].to_numpy(float)])
        if np.linalg.matrix_rank(Xf) < Xf.shape[1]:
            rec[f'{label}_chisq'] = np.nan; continue
        bn, lln, cn = irls(Xn, y)
        bf, llf, cf = irls(Xf, y)
        stat = max(0.0, 2.0 * (llf - lln))
        rec[f'{label}_chisq'] = stat
        rec[f'{label}_p'] = float(chi2.sf(stat, 2))
        rec[f'{label}_conv'] = bool(cn and cf)
        rec[f'{label}_beta_gf'] = float(bf[-2]); rec[f'{label}_beta_gm'] = float(bf[-1])
    rows.append(rec)
    print('done', v, flush=True)

res = pd.DataFrame(rows)
res['chisq_drop'] = res['A_shared_chisq'] - res['C_sexPC_chisq']
res['pct_drop'] = 100 * res['chisq_drop'] / res['A_shared_chisq']
res.to_csv(os.path.join(OUT, 'REFIT_RESULTS.csv'), index=False)
GW = chi2.isf(5e-8, 2)
print()
print('genome-wide 2df threshold chisq =', round(GW, 4))
print(res[['variant', 'maf', 'n', 'A_shared_chisq', 'C_sexPC_chisq', 'pct_drop']].to_string(index=False))
print()
print('ge 5e-8 under A_shared :', int((res.A_shared_chisq >= GW).sum()), 'of', len(res))
print('ge 5e-8 under C_sexPC  :', int((res.C_sexPC_chisq >= GW).sum()), 'of', len(res))
print('median chisq A_shared:', round(float(res.A_shared_chisq.median()), 4))
print('median chisq C_sexPC :', round(float(res.C_sexPC_chisq.median()), 4))
print('median pct drop      :', round(float(res.pct_drop.median()), 2))
