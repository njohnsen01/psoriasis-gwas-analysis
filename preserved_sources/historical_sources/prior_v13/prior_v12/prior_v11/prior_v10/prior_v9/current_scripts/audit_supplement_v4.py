from pathlib import Path
import json, hashlib, math, re, itertools
import numpy as np
import pandas as pd
from scipy.stats import chi2
from openpyxl import load_workbook
from datetime import datetime, timezone
I=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912'); N=I.parent; C=I/'center_refit_v2'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(a,b):
    np.testing.assert_allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),rtol=2e-12,atol=1e-300,equal_nan=True)
def truth(x):
    return str(x).lower() in ('true','1','1.0')
d=pd.read_csv(C/'estimates.private.tsv',sep='	'); pairs=pd.read_csv(C/'paired_comparison.private.tsv',sep='	')
assert len(d)==16 and not d.duplicated(['adjustment','model','ID']).any()
q=pd.read_csv(I/'estimate_audit/target_QC.private.tsv.gz',sep='	')
q=q[(q.group=='pooled') & (q.domain=='X')].set_index('ID')
assert q.index.is_unique
checked_sources={}
for (adj,model),v in d.groupby(['adjustment','model']):
    prefix=C/(adj+'_'+model)
    cmd=json.loads(Path(str(prefix)+'.command.json').read_text())
    assert ('--skip-dosage-comp' in cmd)==(model=='male01')
    assert cmd[cmd.index('--par-region')+1]=='hg38'
    assert cmd[cmd.index('--step')+1]=='2'
    expected_cov=C/'center_covariates.private.tsv' if adj=='center_refit' else N/'inputs/covariates.tsv'
    expected_pred=C/'step1_pred.list' if adj=='center_refit' else N/'step1/model_pred.list'
    assert Path(cmd[cmd.index('--covarFile')+1])==expected_cov
    assert Path(cmd[cmd.index('--pred')+1])==expected_pred
    log=Path(str(prefix)+'.log').read_text()
    assert re.search(r'Number of failed tests\s*:\s*\(0/4\)',log)
    for _,r in v.iterrows():
        source=Path(r.source_file)
        assert sha(source)==r.source_sha256
        checked_sources[str(source)]=sha(source)
        raw=pd.read_csv(source,sep=r'\s+')
        a=raw[(raw.ID==r.ID)&(raw.TEST==r.TEST)];assert len(a)==1
        a=a.iloc[0]
        for col in ('ALLELE0','ALLELE1','TEST'): assert str(a[col])==str(r[col])
        close(a[['BETA','SE','LOG10P','N','N_CASES','N_CONTROLS']],r[['BETA','SE','LOG10P','N','N_CASES','N_CONTROLS']])
        close([r.OR,r.CI95_low,r.CI95_high,r.P_from_REGENIE_LOG10P],[math.exp(r.BETA),math.exp(r.BETA-1.959963984540054*r.SE),math.exp(r.BETA+1.959963984540054*r.SE),10**(-r.LOG10P)])
        assert truth(r.primary_QC_PASS)==truth(q.loc[r.ID,'PASS'])
        for col in ('PASS','pooled_PASS'):
            if pd.notna(q.loc[r.ID,col]): assert truth(r.primary_QC_PASS)==truth(q.loc[r.ID,col])
for _,r in pairs.iterrows():
    a=d[(d.ID==r.ID)&(d.model==r.model)&(d.adjustment=='matched')].iloc[0]
    b=d[(d.ID==r.ID)&(d.model==r.model)&(d.adjustment=='center_refit')].iloc[0]
    assert a.N==b.N and a.ALLELE1==b.ALLELE1 and a.ALLELE0==b.ALLELE0
    for col in ('BETA','SE','LOG10P','N','OR','CI95_low','CI95_high','P_from_REGENIE_LOG10P'):
        close(r[col+'_matched'],a[col]);close(r[col+'_center'],b[col])
    np.testing.assert_allclose(r.relative_OR_change,b.OR/a.OR-1,rtol=2e-12,atol=1e-15)
s=pd.read_csv(I/'sex_heterogeneity_current/sex_heterogeneity.private.tsv',sep='	')
x=pd.read_csv(I/'followup_model_audit/with_QKI_source_QC/X_reconciled.private.tsv',sep='	')
assert len(s)==16 and not s.duplicated(['subset','male_scale','variant']).any()
assert set(map(tuple,s[['subset','male_scale','variant']].values))==set(itertools.product(['all','pruned'],[1,2],d.ID.unique()))
assert (s.family_size==16).all()
for _,r in s.iterrows():
    for sex in ('female','male'):
        source=Path(r[sex+'_source_output']);assert sha(source)==r[sex+'_source_sha256'];checked_sources[str(source)]=sha(source)
        a=x[(x.subset==r.subset)&(x.variant==r.variant)&(x.stratum==sex)&(x.source_output==str(source))&(x.male_scale==r.male_scale)]
        assert len(a)==1,(sex,r.variant,len(a))
        a=a.iloc[0]
        for col in ('beta','se','OR','CI_low','CI_high','p','n','cases','controls'):close(r[sex+'_'+col],a[col])
        assert r[sex+'_P_method']==a.P_method and r[sex+'_CI_method']==a.CI_method
        assert truth(r[sex+'_converged'])==truth(a.converged)
    if truth(r.estimable):
        delta=r.male_beta-r.female_beta; se=np.hypot(r.male_se,r.female_se); Q=(delta/se)**2
        close([r.beta_male_minus_female,r.SE_difference,r.Q,r.P_heterogeneity],[delta,se,Q,chi2.sf(Q,1)])
        close([r.male_carrier_OR,r.male_carrier_CI_low,r.male_carrier_CI_high],np.array([r.male_OR,r.male_CI_low,r.male_CI_high])**r.male_scale)
p=np.where(s.estimable,s.P_heterogeneity,1);order=np.argsort(p);bh=np.empty(16);bh[order]=np.minimum(1,np.minimum.accumulate((p[order]*16/np.arange(1,17))[::-1])[::-1])
close(s.P_heterogeneity_Bonferroni_16,np.minimum(1,p*16));close(s.P_heterogeneity_BH_16,bh)
m=pd.read_csv(I/'supplement_coverage/original57_to_current_sources.private.tsv',sep='	',keep_default_na=False)
assert len(m)==57 and m.original_sheet.is_unique
original=Path('/home/jupyter/S1_current57_preserved_20260910.xlsm'); oh=sha(original)
assert set(m.original_sha256)=={oh};assert set(m.original_workbook)=={str(original)}
w=load_workbook(original,read_only=True,keep_vba=True);assert list(m.sort_values('original_order').original_sheet)==w.sheetnames;w.close()
checked_sources.update({str(p):sha(p) for p in [C/'estimates.private.tsv',C/'paired_comparison.private.tsv',I/'sex_heterogeneity_current/sex_heterogeneity.private.tsv',I/'sex_heterogeneity_current/PROTOCOL.private.json',I/'supplement_coverage/original57_to_current_sources.private.tsv',I/'estimate_audit/target_QC.private.tsv.gz']})
summary=[]
for _,r in pairs.iterrows():
    summary.append({'variant':r.ID,'model':r.model,'primary_QC_PASS':truth(r.primary_QC_PASS_center),'OR_matched':r.OR_matched,'CI_matched':[r.CI95_low_matched,r.CI95_high_matched],'P_matched':r.P_from_REGENIE_LOG10P_matched,'OR_center':r.OR_center,'CI_center':[r.CI95_low_center,r.CI95_high_center],'P_center':r.P_from_REGENIE_LOG10P_center,'relative_OR_change':r.relative_OR_change})
out=I/'supplement_review_v4';out.mkdir(exist_ok=True)
receipt={'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),'checks':['raw REGENIE rows, alleles and N','same-output OR and Wald CI; REGENIE P retained','male scaling and covariate/prediction commands','all four Step2 logs have zero failed tests','paired identical N and effect allele','primary QC flags independently reconciled','sex source models, fixed family, statistics, corrections and carrier scaling','original 57-sheet identity and preserved workbook hash'],'source_hashes':checked_sources,'center_comparison':summary,'release':'HELD; secure review only; not approval to export','limits':['targeted center sensitivity is not genome-wide recalibration','full-cohort finer batch and read mapping validation incomplete','sex tests exploratory with relatedness assumptions; estimability does not establish power','coverage map does not establish final supplement completeness']}
(out/'INDEPENDENT_AUDIT.private.json').write_text(json.dumps(receipt,indent=2)+chr(10))
print('INDEPENDENT_AUDIT_PASS; receipt saved in secure review directory; no export performed')
