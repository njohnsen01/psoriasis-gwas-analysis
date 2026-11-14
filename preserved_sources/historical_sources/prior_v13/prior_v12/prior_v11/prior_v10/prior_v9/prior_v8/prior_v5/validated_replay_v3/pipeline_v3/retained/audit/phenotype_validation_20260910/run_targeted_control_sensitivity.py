from pathlib import Path
import json, hashlib, warnings
from datetime import datetime, timezone
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.special import expit
from scipy.stats import norm
from threadpoolctl import threadpool_limits
W=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/phenotype_validation_20260910')
df=pd.read_pickle(W/'private_phenotype_analysis.pkl')
G=np.load(W/'private_all_sample_lead_genotypes.npy');H=np.load(W/'private_all_sample_lead_heterozygotes.npy')
leads=[('rs183018930',5998724,'A','G'),('rs5970952',22917619,'C','T'),('rs5964607',67736589,'T','C'),('rs5919457',67842269,'G','T')]
rows=[]
status={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),'model':'binomial logistic regression with age, age squared, sex, age-by-sex, PC1-PC10 and sex-by-PC1-PC10; no LOCO or kinship correction','scope':'24 post-hoc targeted tests; no genome-wide calibration or independent replication'}
(W/'association_status.json').write_text(json.dumps(status,indent=2))
with threadpool_limits(limits=4):
 for definition in ['A_original','B_observed','C_no_drug']:
  d=df.loc[df[definition]].copy();idx=d.source_index.to_numpy(int);male=d.sex.to_numpy(float)==1;y0=d.y.to_numpy(float)
  cov=d[['age','age2','sex','age_sex']+['PC'+str(i) for i in range(1,11)]].to_numpy(float)
  cov=np.c_[cov,d[['PC'+str(i) for i in range(1,11)]].to_numpy(float)*male[:,None]]
  cov=(cov-cov.mean(0))/cov.std(0);assert np.isfinite(cov).all()
  for j,(rs,pos,ref,ea) in enumerate(leads):
   raw=G[j,idx].astype(float);valid=(raw>=0)&~(male&H[j,idx]);y=y0[valid];sx=male[valid];gt=raw[valid];c=cov[valid]
   biological=gt/np.where(sx,2.,1.);ploidy=np.where(sx,1.,2.)
   acase=biological[y==1].sum();an_case=ploidy[y==1].sum();actrl=biological[y==0].sum();an_ctrl=ploidy[y==0].sum()
   mac_case=min(acase,an_case-acase);mac_control=min(actrl,an_ctrl-actrl);mac_total=min(acase+actrl,an_case+an_ctrl-acase-actrl)
   cr_case=float(valid[y0==1].mean());cr_control=float(valid[y0==0].mean())
   assert mac_case>=20 and mac_control>=20 and mac_total>=100 and min(cr_case,cr_control)>=.98
   for coding in ['male01','male02']:
    xg=gt/np.where(sx,2.,1.) if coding=='male01' else gt.copy()
    X=np.c_[np.ones(len(y)),c,xg];assert np.linalg.matrix_rank(X.T@X)==X.shape[1]
    with warnings.catch_warnings(record=True) as caught:
     warnings.simplefilter('always');fit=sm.GLM(y,X,family=sm.families.Binomial()).fit(maxiter=100,tol=1e-10)
    beta=float(fit.params[-1]);se=float(fit.bse[-1]);mu=expit(X@fit.params);score=float(np.max(np.abs(X.T@(y-mu))))
    assert fit.converged and np.isfinite([beta,se]).all() and se>0 and score<1e-5
    independent_se=float(np.sqrt(np.linalg.inv(X.T@((mu*(1-mu))[:,None]*X))[-1,-1]))
    assert abs(independent_se-se)<1e-7
    out={'definition':definition,'coding':coding,'rsid':rs,'chr':'X','pos':pos,'reference':ref,'effect_allele':ea,'n':len(y),'cases':int(y.sum()),'controls':int((y==0).sum()),'beta':beta,'se':se,'OR':float(np.exp(beta)),'CI_low':float(np.exp(beta-1.96*se)),'CI_high':float(np.exp(beta+1.96*se)),'P':float(2*norm.sf(abs(beta/se))),'case_MAC':float(mac_case),'control_MAC':float(mac_control),'total_MAC':float(mac_total),'case_call_rate':cr_case,'control_call_rate':cr_control,'converged':bool(fit.converged),'score_max_abs':score,'SE_independently_checked':True,'warnings':[str(w.message) for w in caught]}
    rows.append(out);pd.DataFrame(rows).to_csv(W/'targeted_control_sensitivity.private.tsv',sep='\t',index=False)
    print('Completed',definition,coding,rs,flush=True)
status.update(status='complete',completed_utc=datetime.now(timezone.utc).isoformat(),models=len(rows),all_converged=all(x['converged'] for x in rows),all_SE_independently_checked=True,python_package_versions={'numpy':np.__version__,'pandas':pd.__version__,'statsmodels':sm.__version__},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(W/'association_status.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status),flush=True)
