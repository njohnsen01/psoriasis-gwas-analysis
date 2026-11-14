from pathlib import Path
import json,hashlib,datetime,sys,re
import numpy as np,pandas as pd
from scipy import special,stats
P=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
O=P/'reconciliation_20260912/female_X_source_audit_20260918';O.mkdir(exist_ok=True)
assert not (O/'VERIFIED.json').exists(),'Prior audit exists; inspect instead of overwrite'
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as handle:
  for block in iter(lambda:handle.read(8388608),b''): h.update(block)
 return h.hexdigest()
read=lambda f:json.loads(Path(f).read_text())
checks={}; r=read(P/'results/X_stratified/female/COMPLETE.json')
meta=read(P/'stratified_step1/female/VALIDATED.json')
inp=P/'stratified_step1/female'
checks['source_script_matches_reviewed']=sha(P/'scripts/run_stratified_X.py')=='79ada6faf3e26f1b0cc1d3c2c6436ab7c97686f3ef01633d824d888f880e5c68'
checks['step1_script_matches_reviewed']=sha(P/'scripts/refit_stratified_step1.py')=='bdbe50f7042f9ffab992d28ea94b8eba8c159c93b9a10dc3790ab98dfac16220'
checks['female_covariates']=r['covariates']==['age','age2']+['PC'+str(i) for i in range(1,11)]
for name,key in [('phenotype.tsv','phenotype_sha256'),('covariates.tsv','covariate_sha256'),('model_1.loco','loco_sha256')]: checks['retained_'+name]=sha(inp/name)==meta[key]
ph=pd.read_csv(inp/'phenotype.tsv',sep='\t',dtype={'FID':str,'IID':str}).set_index('IID')
main=pd.read_csv(P/'inputs/phenotype.tsv',sep='\t',dtype={'FID':str,'IID':str}).set_index('IID')
ps=pd.read_csv(P/'inputs/provider_X_aligned.psam',sep='\t',dtype=str).set_index('IID')
checks['unique_female_ids']=ph.index.is_unique
checks['subset_final_primary']=ph.index.isin(main.index).all()
checks['phenotype_matches_primary']=ph.psoriasis_case.equals(main.loc[ph.index,'psoriasis_case'])
checks['all_provider_XX']=ps.loc[ph.index,'SEX'].astype(int).eq(2).all()
checks['female_cohort_equals_all_primary_XX']=set(ph.index)==set(main.index.intersection(ps.index[ps.SEX.astype(int).eq(2)]))
checks['cohort_receipt_counts_match']=len(ph)==meta['n'] and int(ph.psoriasis_case.sum())==meta['cases']
with (inp/'model_1.loco').open() as h:
 header=h.readline().split();ids=set(ph.FID+'_'+ph.index); rows=[];finite=True
 checks['loco_ids_match_female']=set(header[1:])==ids
 for line in h:
  fields=line.split();rows.append(int(fields[0]));a=np.asarray(fields[1:],float);finite=finite and len(a)==len(ph) and np.isfinite(a).all()
checks['loco_23_rows_finite']=rows==list(range(1,24)) and finite and r['LOCO_row']==23
paths=sorted((P/'results/X_stratified/female').glob('chunk_*/COMPLETE.json'))
checks['all_expected_chunks']=len(paths)==r['chunks']==134
wanted=['chrX:5995921:G:A','chrX:5996494:GCT:G','chrX:5998724:A:G']
parts=[];targets=[];source_receipts=[];failed=0; log_failed=0;target_chunks=[]
for f in paths:
 m=read(f);out=f.parent;chunk=out.name
 assert m['status']=='validated'
 if not m['QC_passed']: continue
 s=m['models']['male02'];rp=out/'male02_psoriasis_case.regenie'
 assert sha(rp)==s['results_sha256'],'Association checksum mismatch'
 d=pd.read_csv(rp,sep=r'\s+');assert d.ID.is_unique and len(d)==m['QC_passed'] and d.TEST.eq('ADD').all()
 assert set(d.ID)==set((out/'passed_variants.txt').read_text().split())
 vals=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');ok=np.isfinite(vals).all(axis=1)&vals.SE.gt(0)&vals.LOG10P.ge(0)
 assert int(ok.sum())==s['estimated'] and int((~ok).sum())==s['failed_or_unavailable']
 failed+=int((~ok).sum())
 log=(out/'male02.log').read_text();assert 'LOCO predictions' in log and 'End time' in log and str(inp/'model_pred.list') in log
 ff=re.findall(r'Number of failed tests\s*:\s*\(?(\d+)',log);assert ff,'Missing failed-test diagnostic';log_failed+=sum(map(int,ff))
 assert '--skip-dosage-comp' not in log
 parts.append(d.loc[ok,['ID','GENPOS','LOG10P']])
 sub=d[d.ID.isin(wanted)].copy()
 if len(sub):
  targets.append(sub[['ID','GENPOS','ALLELE0','ALLELE1','BETA','SE','LOG10P','TEST','EXTRA']]);target_chunks.append(chunk)
  src=P/'X_projected_source'/chunk
  for name,key in [('input.pgen','PGEN_sha256'),('input.pvar','PVAR_sha256'),('input.psam','PSAM_sha256')]: assert sha(src/name)==m[key],'Target input checksum mismatch'
  q=pd.read_csv(out/'final_variant_QC.private.tsv.gz',sep='\t').set_index('ID');assert q.loc[sub.ID,'PASS'].all()
  pv=pd.read_csv(src/'input.pvar',sep='\t').set_index('ID')
  assert (pv.loc[sub.ID,'ALT'].to_numpy()==sub.ALLELE1.to_numpy()).all()
 source_receipts.append({'file':str(rp),'sha256':s['results_sha256'],'valid_variants':int(ok.sum()),'failed_variants':int((~ok).sum())})
full=pd.concat(parts,ignore_index=True);assert full.ID.is_unique
v=full.LOG10P.to_numpy(float)
lam=lambda a:float(np.median(special.ndtri_exp(-np.asarray(a)*np.log(10)-np.log(2))**2)/stats.chi2.ppf(.5,1))
cal=read(P/'diagnostics/X_stratified_complete/female/male02_calibration.json')
checks['full_variant_count_matches']=len(full)==r['QC_passed']==cal['N_valid_tests']
checks['lambda_matches_saved']=abs(lam(v)-cal['lambda_GC_from_full_output_P'])<1e-12
checks['no_missing_or_failed_tests']=failed==0 and log_failed==0 and cal['N_failed_or_unavailable']==0
qq=pd.read_csv(P/'diagnostics/X_stratified_complete/female/male02_QQ_full_ranks.tsv.gz',sep='\t')
ranks=qq.full_output_rank.to_numpy(int);sv=np.sort(v)[::-1]
checks['QQ_full_ranks_unit_weight']=np.allclose(qq.observed_log10P,sv[ranks-1],atol=1e-12) and np.allclose(qq.expected_log10P,-np.log10((ranks-.5)/len(v)),atol=1e-12)
t=pd.concat(targets,ignore_index=True).set_index('ID').loc[wanted].reset_index();t['OR']=np.exp(t.BETA);t['CI_low_Wald']=np.exp(t.BETA-1.959963984540054*t.SE);t['CI_high_Wald']=np.exp(t.BETA+1.959963984540054*t.SE);t['P_REGENIE']=10.**(-t.LOG10P)
checks['three_targets_available']=len(t)==3 and t.ID.is_unique
checks={k:bool(v) for k,v in checks.items()};assert all(checks.values()),[k for k,v in checks.items() if not v]
report={'verified_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'valid_variant_tests':len(full),'failed_or_unavailable_variant_tests':failed,'log_failed_variant_tests':log_failed,'lambda_GC':lam(v),'lambda_excluding_5_5_to_6_5_Mb':lam(full.loc[~full.GENPOS.between(5500000,6500000),'LOG10P']),'significant_at_5e_8':int((v>-np.log10(5e-8)).sum()),'minimum_P':float(10.**(-v.max())),'covariates':r['covariates'],'target_chunks':target_chunks,'targets':t.to_dict(orient='records'),'scope':'Post-hoc female-only full-X additive REGENIE with sex-specific refitted LOCO; within-cohort analysis, not external replication; female diploid 0/1/2; CI Wald, P original approximate-Firth REGENIE.'}
t.to_csv(O/'female_scan_target_estimates.tsv',sep='\t',index=False)
(O/'source_receipts.json').write_text(json.dumps(source_receipts,indent=2))
(O/'VERIFIED.json').write_text(json.dumps(report,indent=2))
print('FEMALE_X_SOURCE_AUDIT_VERIFIED')
print(json.dumps({k:report[k] for k in ['valid_variant_tests','failed_or_unavailable_variant_tests','lambda_GC','lambda_excluding_5_5_to_6_5_Mb','significant_at_5e_8','minimum_P']},indent=2))
print(t[['ID','OR','CI_low_Wald','CI_high_Wald','P_REGENIE','EXTRA']].to_string(index=False))

