from pathlib import Path
import json, hashlib, datetime, itertools
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912')
O=I/'sex_heterogeneity_current'; O.mkdir(exist_ok=True)
source=I/'followup_model_audit/with_QKI_source_QC/X_reconciled.private.tsv'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
f=pd.read_csv(source,sep=chr(9))
f=f[(f['condition']=='none') & f['stratum'].isin(['female','male'])].copy()
keys=['subset','male_scale','variant','stratum']
assert not f.duplicated(keys).any()
assert set(f['subset'])=={'all','pruned'} and set(f['male_scale'])=={1,2} and f['variant'].nunique()==4
assert len(f)==32
protocol=dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), source=str(source), source_sha256=sha(source), family='All four historical targets x two male dosage codings x two participant sets; 16 attempted exploratory contrasts, including failed-QC targets as diagnostics', family_size=16, null='Male and female log odds ratios per coded target-ALT dosage unit are equal under the specified male dosage coding', statistic='(beta_male-beta_female)^2/(SE_male^2+SE_female^2), asymptotic chi-square with one degree of freedom', correction='Bonferroni and Benjamini-Hochberg across all 16 attempted contrasts; nonestimable contrasts use P=1 for family correction and remain unavailable in output', timing='Post-hoc descriptive analysis using previously selected historical targets and already available corrected model estimates. Not the original three-test family and not a prespecified replication test.', covariance='Assumes zero covariance of sex-specific estimates. Cross-sex familial relationships can violate this in the all-participant set; provider-pruned sensitivity does not exclude all distant relatedness.', units='Female dosage 0/1/2; male dosage 0/2 or 0/1 as labeled. Contrast is per coded unit. Male carrier OR/CI additionally transformed from source OR/CI using the male scale.', limitations=['Overlapping discovery-selected follow-up is exploratory', 'Historical QC-failed targets remain diagnostic', 'Estimability is not power', 'Wald heterogeneity contrast is separate from source model association P', 'No conditional independence or X-inactivation mechanism inferred'], release='HELD')
(O/'PROTOCOL.private.json').write_text(json.dumps(protocol,indent=2))
rows=[]
for (subset,scale,variant),g in f.groupby(['subset','male_scale','variant'],sort=True):
    assert set(g['stratum'])=={'female','male'} and len(g)==2
    d={r['stratum']:r for _,r in g.iterrows()}; a=d['female']; b=d['male']
    assert a['EA']==b['EA'] and a['primary_QC_status']==b['primary_QC_status']
    rec=dict(subset=subset,male_scale=int(scale),variant=variant,lead=a['lead'],EA=a['EA'],primary_QC_status=a['primary_QC_status'],family_size=16,analysis='post-hoc exploratory Wald sex heterogeneity',relatedness=a['relatedness'],interpretation='diagnostic only; failed primary QC' if a['primary_QC_status']!='passed' else 'exploratory; not replication')
    for sex,r in d.items():
        for col in ['status','converged','n','cases','controls','beta','se','OR','CI_low','CI_high','p','P_method','CI_method','source_output','source_sha256']:
            rec[sex+'_'+col]=r[col]
    ok=all(r['status']=='estimated' and r['converged']==True and np.isfinite(r['beta']) and np.isfinite(r['se']) and r['se']>0 for r in [a,b])
    rec['estimable']=bool(ok)
    if ok:
        delta=float(b['beta']-a['beta']); variance=float(a['se']**2+b['se']**2); z=delta/np.sqrt(variance)
        rec.update(beta_male_minus_female=delta,SE_difference=np.sqrt(variance),Q=z*z,df=1,P_heterogeneity=float(chi2.sf(z*z,1)),male_carrier_OR=float(b['OR']**scale),male_carrier_CI_low=float(b['CI_low']**scale),male_carrier_CI_high=float(b['CI_high']**scale))
        assert np.isclose(rec['P_heterogeneity'],2*norm.sf(abs(z)),rtol=1e-11,atol=1e-15)
    else:
        rec['unavailable_reason']='At least one source sex model is nonestimable or unconverged; source status retained.'
    rows.append(rec)
out=pd.DataFrame(rows); assert len(out)==16
p=out['P_heterogeneity'].fillna(1).to_numpy(); order=np.argsort(p,kind='stable'); ranked=p[order]*len(p)/np.arange(1,len(p)+1); adj=np.minimum.accumulate(ranked[::-1])[::-1]; q=np.empty(len(p)); q[order]=np.minimum(adj,1)
out['P_heterogeneity_Bonferroni_16']=np.minimum(p*16,1); out['P_heterogeneity_BH_16']=q
out.loc[~out['estimable'],['P_heterogeneity_Bonferroni_16','P_heterogeneity_BH_16']]=np.nan
target=O/'sex_heterogeneity.private.tsv'; out.to_csv(target,sep=chr(9),index=False)
check=pd.read_csv(target,sep=chr(9)); assert len(check)==16 and set(check['lead'])==set(f['lead'])
for col in ['beta_male_minus_female','SE_difference','Q','P_heterogeneity','P_heterogeneity_Bonferroni_16','P_heterogeneity_BH_16']:
    assert np.allclose(check[col],out[col],rtol=1e-12,atol=1e-15,equal_nan=True)
receipt=dict(status='computed_and_verified',source_sha256=sha(source),table_sha256=sha(target),script_sha256=sha(__file__),attempted_contrasts=len(out),estimable_contrasts=int(out['estimable'].sum()),family_size=16,source_models_unchanged=True,original_three_test_family_not_redefined=True,release='HELD')
(O/'COMPLETE.private.json').write_text(json.dumps(receipt,indent=2)); print(json.dumps(receipt))