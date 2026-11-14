from pathlib import Path
import datetime,hashlib,json,re
import numpy as np,pandas as pd,openpyxl
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';O=I/'conditional_current_v1'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
c=json.loads((O/'COMPLETE.private.json').read_text());assert c['status']=='complete_secure_conditional_sensitivity';assert not c['fit_failures'];assert c['conditional_output_rows']==22
for name,digest in c['outputs'].items():assert sha(O/name)==digest
a=pd.read_csv(O/'conditional_estimates.private.tsv',sep='\t');base=pd.read_csv(O/'frozen_primary_targets.private.tsv',sep='\t').set_index('ID');b=pd.read_csv(O/'unconditioned_concordance.private.tsv',sep='\t');assert b.passed.all() and b.absolute_difference.eq(0).all();assert a.estimable.all()
checks=[]
for label,g in a.groupby('analysis'):
 d=O/'fits'/label;cmd=json.loads((d/'console.command.json').read_text());log=(d/'assoc.log').read_text();conds=(d/'condition.txt').read_text().splitlines();targets=(d/'targets.txt').read_text().splitlines();assert set(g.ID)==set(targets)-set(conds)
 assert re.search(r'n_cov\s*=\s*14\b',log);used=re.search(r'\+conditioning on variants.*?n_used\s*=\s*(\d+)',log,re.S);assert used and int(used.group(1))==len(conds)
 fail=re.search(r'Number of failed tests : \((\d+)/(\d+)\)',log);assert fail and int(fail.group(1))==0;ignored=re.search(r'Number of ignored tests due to low MAC : (\d+)',log);assert ignored and int(ignored.group(1))==0
 warnings=[x for x in log.splitlines() if re.search(r'warning|collinear|singular|not converg|failed to|dropping.*covariate|removed.*covariate',x,re.I)];assert not warnings,(label,warnings)
 raw=pd.read_csv(d/'assoc_psoriasis_case.regenie',sep=r'\s+').set_index('ID')
 for _,r in g.iterrows():
  v=r.ID
  for k in ['BETA','SE','LOG10P','N','ALLELE0','ALLELE1','CHROM','GENPOS','TEST']:assert r[k]==raw.loc[v,k],(label,v,k)
  assert r.ALLELE0==base.loc[v,'ALLELE0'] and r.ALLELE1==base.loc[v,'ALLELE1'];assert r.N==base.loc[v,'N'];assert np.isclose(r.OR,np.exp(r.BETA));assert np.isclose(r.CI95_Wald_low,np.exp(r.BETA-1.959963984540054*r.SE));assert np.isclose(r.CI95_Wald_high,np.exp(r.BETA+1.959963984540054*r.SE));assert np.isclose(r.P_REGENIE,10.**(-r.LOG10P),rtol=1e-12,atol=0)
 checks.append(dict(analysis=label,conditional_variants_requested=len(conds),conditional_variants_used=int(used.group(1)),original_numeric_covariates=14,failed_tests=int(fail.group(1)),low_MAC_ignored=int(ignored.group(1)),warnings='none',output_rows=len(g),log_sha256=sha(d/'assoc.log'),command_sha256=sha(d/'console.command.json')))
pd.DataFrame(checks).to_csv(O/'model_diagnostics.private.tsv',sep='\t',index=False)
ld=pd.read_csv(O/'ancestry_control_LD.private.tsv',sep='\t');assert len(ld)==60;assert ld[ld.available].r2.between(0,1).all();assert ld.loc[~ld.available,'r2'].isna().all()
D=pd.read_pickle('/home/jupyter/psoriasis_reproducibility_20260910/reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');D.index=D.index.astype(str);assert D.unrelated.dropna().isin([True,False]).all()
ph=pd.read_csv(N/'inputs/phenotype.tsv',sep='\t',dtype={'IID':str}).set_index('IID');keep=pd.read_csv(N/'inputs/keep.tsv',sep='\t',header=None,names=['FID','IID'],dtype=str);D=D.loc[D.index.intersection(keep.IID)];assert D.psoriasis_case.eq(ph.loc[D.index,'psoriasis_case']).all();anc=D.ancestry_pred.astype(str).str.upper()
for group in ['EUR','AFR','AMR']:
 members=D.index[D.unrelated.astype(bool)&D.psoriasis_case.eq(0)&anc.eq(group)];reference=pd.read_csv(N/'inputs'/(group+'.controls.keep.tsv'),sep='\t',header=None,names=['FID','IID'],dtype=str);assert set(members).issubset(set(reference.IID));assert ld.loc[ld.ancestry.eq(group),'reference_N'].eq(len(members)).all()
summary=a.groupby('family').agg(attempted=('ID','size'),estimable=('estimable','sum'),family_significant=('family_significant','sum'),genome_wide_significant=('genome_wide_significant','sum')).reset_index();summary.to_csv(O/'conditional_summary.private.tsv',sep='\t',index=False)
ldsummary=ld.groupby(['region','ancestry']).agg(pair_count=('r2','size'),available=('available','sum'),r2_min=('r2','min'),r2_max=('r2','max')).reset_index();ldsummary.to_csv(O/'LD_summary.private.tsv',sep='\t',index=False)
receipt=dict(status='source_model_and_summary_verified',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),baseline_printed_estimates_exact=True,conditional_models=len(checks),conditional_estimates=22,all_conditioning_variants_used=True,all_models_zero_failed_tests=True,LD_rows=60,LD_reference_verified_against_current_primary=True,script_sha256=sha(__file__),source_complete_sha256=sha(O/'COMPLETE.private.json'),outputs={p.name:sha(p) for p in O.glob('*.private.tsv')})
(O/'AUDIT.private.json').write_text(json.dumps(receipt,indent=2));print('CONDITIONAL_AUDIT_VERIFIED');print(summary.to_string(index=False));print(ldsummary.to_string(index=False))
