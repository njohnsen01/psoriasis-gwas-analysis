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
with threadpool_limits(limits=2):
 for subset in ['all','pruned']:
  d=D if subset=='all' else D.loc[D.unrelated]
  for stratum in ['pooled','female','male','EUR','AFR','AMR']:
   z=d if stratum=='pooled' else d.loc[d.SEX.eq(2 if stratum=='female' else 1)] if stratum in ['female','male'] else d.loc[d.ancestry_pred.str.upper().eq(stratum)]
   cols=sexcols if stratum in ['female','male'] else pooled
   for scale in ['02','01']:
    for rs,vid,ea in L:fit_logistic(z,rs,None,cols,scale,{'subset':subset,'stratum':stratum,'variant':vid,'EA':ea})
    if stratum not in ['female','male']:
     for rs,other in [('rs5964607','rs5919457'),('rs5919457','rs5964607')]:fit_logistic(z,rs,other,cols,scale,{'subset':subset,'stratum':stratum,'variant':next(x[1] for x in L if x[0]==rs),'EA':next(x[2] for x in L if x[0]==rs)})
   pd.DataFrame(rows).to_csv(R/'X_results.private.tsv',sep=chr(9),index=False)
   print('X_DONE',subset,stratum,len(rows),flush=True)
(R/'X_QC.private.json').write_text(json.dumps(qc,indent=2))
print('X_COMPLETE',len(rows),flush=True)
