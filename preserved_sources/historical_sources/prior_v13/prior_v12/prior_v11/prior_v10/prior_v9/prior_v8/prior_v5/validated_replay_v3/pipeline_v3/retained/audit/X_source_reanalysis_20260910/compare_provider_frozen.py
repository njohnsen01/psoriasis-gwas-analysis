from pathlib import Path
import json,subprocess,os,pandas as pd,numpy as np,hashlib
w=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/X_source_reanalysis_20260910')
old=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}')
p=json.loads((old/'status/step2_plan.json').read_text())
results=[]
for model in ['male02','male01']:
 out=w/('frozen_provider_'+model)
 cmd=['regenie','--step','2','--bt','--pgen',str(w/'private_provider_leads'),'--phenoFile',p['phenotype_file'],'--phenoCol','psoriasis_case','--covarFile',p['covariate_file'],'--bsize','400','--firth','--approx','--pThresh','0.01','--pred',p['prediction_list'],'--af-cc','--par-region','hg38','--threads','2','--out',str(out)]
 if model=='male01':cmd+=['--skip-dosage-comp']
 (w/('frozen_provider_'+model+'_command.json')).write_text(json.dumps(cmd,indent=2))
 with open(str(out)+'.console.log','w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
 assert r.returncode==0
 d=pd.read_csv(str(out)+'_psoriasis_case.regenie',sep=r'\s+');assert len(d)==4
 d['P']=10.**(-d.LOG10P);d['OR']=np.exp(d.BETA);d['CI_low']=np.exp(d.BETA-1.96*d.SE);d['CI_high']=np.exp(d.BETA+1.96*d.SE);d['model']=model
 results.append(d)
pd.concat(results).to_csv(w/'frozen_provider_comparison.tsv',sep='	',index=False)
print('FROZEN_PROVIDER_COMPARISON_COMPLETE')
