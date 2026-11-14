from pathlib import Path
import json, math, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent
D=pd.read_pickle(R/'X_inputs.private.pkl'); B=pd.read_pickle(R/'biomarker_inputs.private.pkl'); M=pd.read_csv(R/'variant_map.tsv',sep=chr(9))
L=[('rs183018930','chrX:5998724:A:G','G'),('rs5970952','chrX:22917619:C:T','T'),('rs5964607','chrX:67736589:T:C','C'),('rs5919457','chrX:67842269:G:T','T')]
pcs=['PC'+str(i) for i in range(1,11)]
for p in pcs:D['sex_'+p]=D.sex*D[p]
pooled=['sex','age','age2','age_sex']+pcs+['sex_'+p for p in pcs]
sexcols=['age','age2']+pcs
rows=[];qc=[]
def fit_logistic(d,rs,condition,cols,scale,meta):
    names=[rs]+([condition] if condition else [])
    z=d[cols+names+['psoriasis_case','SEX']].replace([np.inf,-np.inf],np.nan).dropna().copy()
    y=z.psoriasis_case.to_numpy(float);n=len(z);counts=[int((y==v).sum()) for v in [1,0]]
    why=[];macs=[];rates=[]
    for name in names:
        g=z[name].to_numpy(float);male=z.SEX.eq(1).to_numpy();ac=g/np.where(male,2.,1.);an=np.where(male,1.,2.)
        for v in [1,0]:
            use=y==v;a=ac[use].sum();total=an[use].sum();macs.append(min(a,total-a));rates.append(float(d.loc[d.psoriasis_case.eq(v),name].notna().mean()))
        if min(macs[-2:])<20 or sum(macs[-2:])<100:why.append('MAC')
    if min(rates)<.98:why.append('call_rate')
    if counts[0]<200 or counts[1]<500:why.append('sample_size')
    rec={**meta,'lead':rs,'condition':condition or 'none','male_scale':scale,'n':n,'cases':counts[0],'controls':counts[1]}
    qc.append({**rec,'MACs':macs,'call_rates':rates,'reason':','.join(sorted(set(why)))})
    if why:rows.append({**rec,'status':'not_estimable','reason':','.join(sorted(set(why)))});return
    x=z[cols].astype(float);sd=x.std(ddof=0);assert (sd>0).all();x=(x-x.mean())/sd
    for name in names:x[name]=z[name]/np.where(z.SEX.eq(1)&(scale=='01'),2.,1.)
    x=sm.add_constant(x);rank=int(np.linalg.matrix_rank(x));assert rank==x.shape[1]
    f=sm.GLM(y,x,family=sm.families.Binomial()).fit(maxiter=100,tol=1e-8)
    b=float(f.params[rs]);se=float(f.bse[rs]);p=float(f.pvalues[rs]);assert f.converged and np.isfinite([b,se,p]).all() and se>0
    rows.append({**rec,'status':'estimated','beta':b,'se':se,'OR':float(np.exp(b)),'CI_low':float(np.exp(b-1.95996398454*se)),'CI_high':float(np.exp(b+1.95996398454*se)),'p':p,'rank':rank,'converged':bool(f.converged)})
D=pd.read_pickle(R/'ancestry_autosomal_inputs.private.pkl');D['SEX']=2
A=M.loc[~M.variant_id.str.startswith('chrX:')].copy()
A=pd.concat([A,pd.DataFrame([{'variant_id':'chr6:163803171:CTTATT:C','aligned_column':'QKI','effect_allele':'C'}])],ignore_index=True)
assert len(A)==22
cols=['sex','age','age2','age_sex']+pcs
with threadpool_limits(limits=2):
 for subset in ['all','pruned']:
  d=D if subset=='all' else D.loc[D.unrelated]
  for stratum in ['EUR','AFR','AMR']:
   z=d.loc[d.ancestry_pred.str.upper().eq(stratum)]
   for _,v in A.iterrows():fit_logistic(z,v.aligned_column,None,cols,'02',{'subset':subset,'stratum':stratum,'variant':v.variant_id,'EA':v.effect_allele})
   pd.DataFrame(rows).to_csv(R/'ancestry_autosomal_results.private.tsv',sep=chr(9),index=False)
   print('ANCESTRY_DONE',subset,stratum,len(rows),flush=True)
(R/'ancestry_autosomal_QC.private.json').write_text(json.dumps(qc,indent=2))
print('ANCESTRY_COMPLETE',len(rows),flush=True)
