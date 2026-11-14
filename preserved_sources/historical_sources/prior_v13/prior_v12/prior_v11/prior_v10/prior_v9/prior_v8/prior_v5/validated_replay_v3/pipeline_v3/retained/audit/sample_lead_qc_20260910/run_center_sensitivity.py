from pathlib import Path
import os,json,subprocess,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;X=R.parent/'X_source_reanalysis_20260910'
d=pd.read_csv(R/'sample_audit.private.tsv',sep='	');c=pd.read_csv(X/'private_inputs/pooled.covar.tsv',sep='	');d=d[d.strict_eligible]
c=c.merge(d[['IID','site_id']],on='IID',validate='one_to_one')
for i in range(1,11):c['sex_PC'+str(i)]=c.sex*c['PC'+str(i)]
for name in ['uw','bi']:c['center_'+name]=(c.site_id==name).astype(int)
assert c.site_id.isin(['bcm','uw','bi']).all()
c=c.drop(columns='site_id');c.to_csv(R/'center_covariates.private.tsv',sep='	',index=False)
c.drop(columns=['center_uw','center_bi']).to_csv(R/'matched_covariates.private.tsv',sep='	',index=False)
for fn in ['center_covariates.private.tsv','matched_covariates.private.tsv']:
 v=pd.read_csv(R/fn,sep='	').drop(columns=['FID','IID']).to_numpy(float);assert np.isfinite(v).all();assert np.linalg.matrix_rank(np.column_stack([np.ones(len(v)),v]))==v.shape[1]+1
pred=json.loads((X/'step1/pooled.validated.json').read_text())['prediction_list'];results=[]
for model in ['male02','male01']:
 for adjustment,cf in [('matched','matched_covariates.private.tsv'),('center','center_covariates.private.tsv')]:
  name=model+'_'+adjustment;out=R/name
  cmd=['regenie','--step','2','--bt','--pgen',str(X/'private_provider_leads'),'--keep',str(R/'strict.keep.private.tsv'),'--phenoFile',str(X/'private_inputs/pooled.pheno.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(R/cf),'--bsize','100','--firth','--approx','--pThresh','0.01','--pred',pred,'--af-cc','--par-region','hg38','--minMAC','5','--threads','4','--out',str(out)]
  if model=='male01':cmd+=['--skip-dosage-comp']
  (R/(name+'.command.json')).write_text(json.dumps(cmd,indent=2))
  with open(R/(name+'.console.log'),'w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
  assert p.returncode==0,name
  t=pd.read_csv(str(out)+'_psoriasis_case.regenie',sep=r'\s+');assert len(t)==4 and np.isfinite(t[['BETA','SE','LOG10P']]).all().all()
  for _,z in t.iterrows():results.append({'model':model,'adjustment':adjustment,'ID':z.ID,'N':int(z.N),'OR':float(np.exp(z.BETA)),'CI_low':float(np.exp(z.BETA-1.96*z.SE)),'CI_high':float(np.exp(z.BETA+1.96*z.SE)),'P':float(10**(-z.LOG10P))})
a=pd.DataFrame(results);m=a[a.adjustment=='matched'].merge(a[a.adjustment=='center'],on=['model','ID'],suffixes=('_base','_center'));assert (m.N_base==m.N_center).all()
(R/'center_sensitivity.private.json').write_text(json.dumps(results,indent=2))
public_results=[dict(z,N='S' if 0<len(c)-z['N']<20 else z['N']) for z in results]
(R/'center_sensitivity.public.json').write_text(json.dumps(public_results,indent=2));(R/'center_summary.json').write_text(json.dumps({'status':'complete','n_models':len(results),'min_P':float(a.P.min()),'max_relative_OR_change':float(abs(m.OR_center/m.OR_base-1).max()),'matched_N':True,'cohort_n':len(c),'LOCO':'coverage cohort retained; no Step 1 refit for stricter exclusions or added center covariates'},indent=2))
