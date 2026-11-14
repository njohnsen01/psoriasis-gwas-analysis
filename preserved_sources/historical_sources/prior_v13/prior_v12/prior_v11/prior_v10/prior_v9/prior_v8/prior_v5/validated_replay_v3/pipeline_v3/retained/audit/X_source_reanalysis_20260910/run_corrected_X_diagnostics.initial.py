from pathlib import Path
import os,sys,json,time,subprocess,traceback,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np,pandas as pd
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;S=json.loads((R/'strata.json').read_text());D=R/'corrected_cohort_diagnostics';D.mkdir(exist_ok=True)
def wait_model(n):
 while not (R/'step1'/(n+'.validated.json')).exists():
  st=json.loads((R/'pipeline_status.json').read_text())
  if st.get('status')=='failed':raise RuntimeError(st.get('error'))
  time.sleep(10)
 return json.loads((R/'step1'/(n+'.validated.json')).read_text())
def job(n,model,scope='leads',condition=None):
 rec=wait_model(n);label=n+'_'+model+'_'+scope+('_condition_'+condition.split(':')[1] if condition else '')
 out=D/label;receipt=D/(label+'.validated.json')
 if receipt.exists():return json.loads(receipt.read_text())
 pref=R/('private_provider_leads' if scope=='leads' else 'private_corrected_X_original_panel')
 cmd=['regenie','--step','2','--bt','--pgen',str(pref),'--keep',str(R/'private_inputs'/(n+'.keep.tsv')),'--phenoFile',str(R/'private_inputs'/(n+'.pheno.tsv')),'--phenoCol','psoriasis_case','--covarFile',str(R/'private_inputs'/(n+'.covar.tsv')),'--bsize','400','--firth','--approx','--pThresh','0.01','--pred',rec['prediction_list'],'--af-cc','--par-region','hg38','--minMAC','5','--threads','4','--out',str(out)]
 if model=='male01':cmd+=['--skip-dosage-comp']
 if condition:
  cfile=D/(label+'.condition.txt');cfile.write_text(condition+'\n');cmd+=['--condition-list',str(cfile)]
 (D/(label+'.command.json')).write_text(json.dumps(cmd,indent=2))
 with open(str(out)+'.console.log','w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1'))
 if r.returncode:raise RuntimeError(label+': '+Path(str(out)+'.console.log').read_text()[-3000:])
 log=Path(str(out)+'.log').read_text()
 assert 'number of individuals used in analysis = '+str(S[n]['n']) in log
 assert str(S[n]['cases'])+' cases and '+str(S[n]['controls'])+' controls' in log
 f=Path(str(out)+'_psoriasis_case.regenie');d=pd.read_csv(f,sep=r'\s+');assert d.ID.is_unique
 valid=np.isfinite(d.LOG10P)&(d.LOG10P>=0);v=d[valid].copy()
 result={'stratum':n,'model':model,'scope':scope,'condition':condition,'output_rows':len(d),'valid_rows':int(valid.sum()),'failed_or_nonfinite_rows':int((~valid).sum()),'output_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'status':'validated_numerical_output','site_panel':'original tested X sites projected from unaffected provider PGEN; newly eligible sites not assessed','post_hoc':True}
 if scope=='panel':
  pieces={'X_all':v,'X_nonPAR':v[(v.GENPOS>2781479)&(v.GENPOS<155701383)],'X_PAR':v[(v.GENPOS<=2781479)|(v.GENPOS>=155701383)]}
  for region,t in pieces.items():
   lp=t.LOG10P.to_numpy();p=np.maximum(np.power(10.,-lp),np.finfo(float).tiny);lam=float(np.median(chi2.isf(p,1))/chi2.ppf(.5,1)) if len(p) else None;result[region]={'n':len(p),'lambda_gc':lam}
   if len(lp):
    observed=np.sort(lp);expected=-np.log10((np.arange(len(lp),0,-1)-.5)/len(lp));fig,ax=plt.subplots(figsize=(5,5));ax.scatter(expected,observed,s=1,rasterized=True);ax.plot([0,max(expected)],[0,max(expected)],color='gray',lw=.8);ax.set(xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title=n+' '+model+' '+region+'\n'+'n='+str(len(lp))+', lambda='+format(lam,'.4f'));fig.tight_layout();fig.savefig(D/(label+'_'+region+'_QQ.pdf'));plt.close(fig)
  result['qq_method']='All finite valid output P values; expected order statistics (rank-0.5)/N; no sampling weights, no enriched subsample. Lambda from all finite valid tests.'
 receipt.write_text(json.dumps(result,indent=2));print('VALIDATED',label,flush=True);return result
def leads(n):
 for m in (['male02'] if n=='female' else ['male01'] if n=='male' else ['male02','male01']):
  job(n,m)
  if n in ['pooled','EUR','AFR','AMR']:
   for c in ['chrX:67736589:T:C','chrX:67842269:G:T']:job(n,m,condition=c)
try:
 with ThreadPoolExecutor(max_workers=6) as ex:
  for f in as_completed([ex.submit(leads,n) for n in S]):f.result()
 (D/'leads_complete.json').write_text(json.dumps({'status':'complete'}))
 while True:
  if (R/'projection_status.json').exists():
   st=json.loads((R/'projection_status.json').read_text())
   if st['status']=='complete':break
  time.sleep(10)
 tasks=[(n,m) for n in S for m in (['male02'] if n=='female' else ['male01'] if n=='male' else ['male02','male01'])]
 with ThreadPoolExecutor(max_workers=6) as ex:
  for f in as_completed([ex.submit(job,n,m,'panel') for n,m in tasks]):f.result()
 (D/'pipeline_complete.json').write_text(json.dumps({'status':'complete'}))
except Exception as e:
 (D/'pipeline_error.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc()},indent=2));raise
