from pathlib import Path
import os,json,sys,hashlib,datetime,subprocess,traceback
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
R=Path('/home/jupyter/psoriasis_reproducibility_20260910')
O=B/'reconciliation_20260912/SII_final_cohort_20260914';O.mkdir(exist_ok=True)
def save(n,x): (O/n).write_text(json.dumps(x,indent=2,default=str))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
try:
 L=pd.read_csv(O/'final_variant_list.tsv',sep='\t');base=pd.read_pickle(O/'matched_inputs.private.pkl');base['ancestry_pred']=base.ancestry_pred.str.upper()
 z=json.loads((O/'scaling.private.json').read_text());scaling=z['scaling'];winsor=z['winsor_bounds'];outcomes=['original_rebuilt','strict_units','platelet_count','ANC','ALC'];cols=['age','age2','SEX','age_sex']+['PC'+str(i) for i in range(1,11)]
 from scipy.stats import norm,chi2
 from statsmodels.stats.multitest import multipletests
 rows=[]
 for mode in ['baseline','winsorized','EUR','AFR','AMR']:
  a=base if mode in ['baseline','winsorized'] else base.loc[base.ancestry_pred.eq(mode)]
  for vid in L.ID:
   d=a.dropna(subset=[vid]);n=len(d);mac=min(d[vid].sum(),2*n-d[vid].sum());nuis=cols+['psoriasis_case'];nuis=[c for c in nuis if d[c].nunique()>1]
   X=np.column_stack([np.ones(n),d[nuis].to_numpy(float),d[vid].to_numpy(float)])
   ok=n>=500 and mac>=20 and np.linalg.matrix_rank(X)==X.shape[1]
   if ok:
    inv=np.linalg.inv(X.T@X);hat=np.einsum('ij,jk,ik->i',X,inv,X);assert hat.max()<1
   for y in outcomes:
    row={'mode':mode,'variant':vid,'effect_allele':str(L.set_index('ID').loc[vid,'ALLELE1']),'outcome':y,'N':n,'MAC':mac,'covariates':','.join(nuis),'status':'estimated' if ok else 'not_estimable'}
    if ok:
     raw=np.log(d[y].to_numpy(float));raw=np.clip(raw,*winsor[y]) if mode=='winsorized' else raw
     yy=(raw-scaling[y]['mean'])/scaling[y]['sd'];beta=inv@X.T@yy;res=yy-X@beta;meat=X.T@((res/(1-hat))[:,None]**2*X);vc=inv@meat@inv;se=np.sqrt(vc[-1,-1]);b=beta[-1];pv=2*norm.sf(abs(b/se));row.update(beta=b,se=se,CI_low=b-1.95996398454*se,CI_high=b+1.95996398454*se,p=pv)
    rows.append(row)
  print('MODE_COMPLETE',mode,flush=True)
 res=pd.DataFrame(rows)
 for label,mask,m in [('main',res['mode'].eq('baseline'),5*len(L)),('robustness',res['mode'].ne('baseline'),20*len(L))]:
  pv=res.loc[mask,'p'].fillna(1).to_numpy();assert len(pv)==m;res.loc[mask,'family']=label;res.loc[mask,'p_bonferroni']=np.minimum(1,pv*m);res.loc[mask,'q_BH']=multipletests(pv,method='fdr_bh')[1]
 res.to_csv(O/'SII_component_results.private.tsv',sep='	',index=False)
 het=[]
 for vid in L.ID:
  for y in outcomes:
   z=res.loc[res['mode'].isin(['EUR','AFR','AMR'])&res.variant.eq(vid)&res.outcome.eq(y)&res.status.eq('estimated')];h={'variant':vid,'outcome':y,'groups':','.join(z['mode']),'status':'estimated' if len(z)>=2 else 'not_estimable'}
   if len(z)>=2:
    w=1/z.se**2;mu=(w*z.beta).sum()/w.sum();Q=(w*(z.beta-mu)**2).sum();h.update(Q=Q,df=len(z)-1,p=chi2.sf(Q,len(z)-1))
   het.append(h)
 h=pd.DataFrame(het);h['p']=h.get('p',pd.Series(np.nan,index=h.index));pv=h.p.fillna(1);h['p_bonferroni']=np.minimum(1,pv*len(h));h['q_BH']=multipletests(pv,method='fdr_bh')[1];h.to_csv(O/'ancestry_heterogeneity.private.tsv',sep='	',index=False)
 save('COMPLETE.json',{'status':'complete','variants':len(L),'planned_model_tests':len(res),'estimated_models':int(res.status.eq('estimated').sum()),'main_bonferroni_significant':int((res.family.eq('main')&res.p_bonferroni.lt(.05)).sum()),'robustness_bonferroni_significant':int((res.family.eq('robustness')&res.p_bonferroni.lt(.05)).sum()),'protocol_sha256':sha(O/'PROTOCOL.json'),'outputs_remain_private':True,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 print('SII_ANALYSIS_COMPLETE',flush=True)
except Exception:
 save('FAILED.json',{'traceback':traceback.format_exc()});traceback.print_exc();raise
