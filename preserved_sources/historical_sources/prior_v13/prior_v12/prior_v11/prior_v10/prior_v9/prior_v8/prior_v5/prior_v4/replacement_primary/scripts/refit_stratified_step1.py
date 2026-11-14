from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,os,shutil,subprocess,traceback
import numpy as np,pandas as pd
from reconcile_X_site_filters import sha,save
N=Path(__file__).resolve().parents[1];O=N/'stratified_step1';O.mkdir(exist_ok=True);ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
def prepare():
 ph=pd.read_csv(N/'inputs/phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str}).set_index('IID')
 cv=pd.read_csv(N/'inputs/covariates.tsv',sep=chr(9),dtype={'FID':str,'IID':str}).set_index('IID')
 ca=pd.read_pickle(N/'cohort_audit/cohort_clinical.private.pkl').reindex(ph.index);assert ca.strict.all()
 assert len(ph)==386069 and ph.psoriasis_case.sum()==5705 and ph.index.is_unique and ph.index.equals(cv.index)
 groups={};pcs=['PC'+str(i) for i in range(1,11)]
 for label in ['AMR','AFR','female','male','EUR']:
  mask=ca.complete_provider_ancestry.eq(label) if label in ['EUR','AFR','AMR'] else ca.SEX.eq(2 if label=='female' else 1)
  p=ph.loc[mask];c=cv.loc[mask,['FID']+(['age','age2']+pcs if label in ['female','male'] else ['sex','age','age2','age_sex']+pcs)]
  values=c.drop(columns='FID').to_numpy(float);assert np.isfinite(values).all();sd=values.std(axis=0);assert (sd>0).all()
  design=np.column_stack([np.ones(len(c)),(values-values.mean(axis=0))/sd]);assert np.linalg.matrix_rank(design)==design.shape[1]
  folder=O/label;folder.mkdir(exist_ok=True);p.reset_index()[['FID','IID','psoriasis_case']].to_csv(folder/'phenotype.tsv',sep=chr(9),index=False)
  c.reset_index()[['FID','IID']+list(c.columns[1:])].to_csv(folder/'covariates.tsv',sep=chr(9),index=False)
  p.reset_index()[['FID','IID']].to_csv(folder/'keep.tsv',sep=chr(9),index=False,header=False)
  cases=int(p.psoriasis_case.sum());assert cases>=200 and len(p)-cases>=500
  groups[label]={'n':len(p),'cases':cases,'controls':len(p)-cases,'covariates':list(c.columns[1:]),'rank_with_intercept':int(design.shape[1]),'phenotype_sha256':sha(folder/'phenotype.tsv'),'covariate_sha256':sha(folder/'covariates.tsv'),'group_assignment':'SEX from frozen sex covariate concordant with provider XX/XY' if label in ['female','male'] else 'complete echo_v4_r2 ancestry_pred, exact label, no probability threshold'}
 assert groups['female']['n']+groups['male']['n']==386069
 save(O/'strata.json',groups);return groups

def fit(label,meta):
 folder=O/label;record=folder/'VALIDATED.json'
 if record.exists():
  r=json.loads(record.read_text());assert r['loco_sha256']==sha(folder/'model_1.loco');return r
 assert shutil.disk_usage(N).free>128*1024**3
 cmd=['regenie','--step','1','--bt','--bed',str(N.parent/'stratified_calibration/private_panel/calibration_pruned'),'--keep',str(folder/'keep.tsv'),'--phenoFile',str(folder/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(folder/'covariates.tsv'),'--bsize','2000','--lowmem','--lowmem-prefix',str(folder/'lowmem'),'--threads','4','--out',str(folder/'model')]
 save(folder/'command.json',cmd)
 with (folder/'console.log').open('w') as handle:
  proc=subprocess.Popen(cmd,stdout=handle,stderr=subprocess.STDOUT,env=ENV);save(folder/'pid.json',{'pid':proc.pid,'status':'running'});rc=proc.wait()
 assert rc==0,'Stratum '+label+' Step 1 failed'
 log=(folder/'model.log').read_text()
 assert 'number of individuals used in analysis = '+str(meta['n']) in log
 assert str(meta['cases'])+' cases and '+str(meta['controls'])+' controls' in log
 keep=pd.read_csv(folder/'keep.tsv',sep=chr(9),header=None,dtype=str);ids=set(keep[0]+'_'+keep[1]);loco=folder/'model_1.loco'
 with loco.open() as handle:
  header=handle.readline().split();assert header[0]=='FID_IID' and len(header)==len(ids)+1 and set(header[1:])==ids
  rows=[]
  for line in handle:
   fields=line.split();rows.append(int(fields[0]));v=np.array(fields[1:],float);assert len(v)==len(ids) and np.isfinite(v).all()
 assert rows==list(range(1,24));assert (folder/'model_pred.list').read_text().split()==['psoriasis_case',str(loco)]
 r={'status':'validated','stratum':label,**meta,'loco_rows':rows,'loco_sha256':sha(loco),'log_sha256':sha(folder/'model.log'),'regenie_binary_sha256':sha(Path(shutil.which('regenie'))),'scope':'new strict-cohort prediction fit; association tests and stratum-specific variant QC still required'}
 save(record,r);return r

def main():
 groups=prepare();save(O/'analysis_addendum.json',{'status':'dated post-hoc strict-cohort calibration refits','created_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'purpose':'sex and EUR/AFR/AMR diagnostics for corrected primary analysis','predictor_panel':'same fixed LD-pruned autosomal calibration_pruned BED as strict pooled Step 1','single_sex_covariates':'age, age squared, PC1-PC10; omit sex and age-by-sex','ancestry_covariates':'sex, age, age squared, age-by-sex, PC1-PC10','independent_replication':False,'model_selection':'retain prespecified dosage comparisons independent of significance','concurrency':'at most two Step 1 fits, four threads each','code_sha256':sha(Path(__file__))})
 done={};errors={}
 with ThreadPoolExecutor(max_workers=2) as pool:
  work={pool.submit(fit,label,m):label for label,m in groups.items()}
  for future in as_completed(work):
   label=work[future]
   try:done[label]=future.result();print('STRICT_STRATUM_STEP1_VALIDATED',label,flush=True)
   except Exception as e:errors[label]={'error':str(e),'traceback':traceback.format_exc()}
   save(O/'status.json',{'status':'running','validated':list(done),'errors':errors,'expected':list(groups)})
 assert not errors,errors
 save(O/'COMPLETE.json',{'status':'all_stratum_step1_validated','strata':list(done),'association_and_variant_QC_pending':True})
 save(O/'status.json',{'status':'complete','validated':list(done),'errors':{},'expected':list(groups)})
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
