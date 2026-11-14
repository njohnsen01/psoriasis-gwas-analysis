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
traits={'SII':'log_inflammatory_index_z','HDL':'HDL_z','LDL':'LDL_z','TG':'log_TG_z'}
def ols(y,g,c):
    valid=np.isfinite(y)&np.isfinite(g)&np.isfinite(c).all(axis=1);y=y[valid];g=g[valid];c=c[valid]
    c=np.column_stack([np.ones(len(y)),c]);rank=np.linalg.matrix_rank(c);assert rank==c.shape[1]
    a=np.linalg.lstsq(c,np.column_stack([y,g]),rcond=None)[0];r=np.column_stack([y,g])-c@a;yr=r[:,0];gr=r[:,1];den=float(gr@gr);assert den>0
    b=float(gr@yr/den);err=yr-b*gr;df=len(y)-rank-1;se=float(np.sqrt((err@err/df)/den));p=float(2*stats.t.sf(abs(b/se),df));crit=stats.t.ppf(.975,df)
    assert np.isfinite([b,se,p]).all() and se>0
    return {'n':len(y),'beta':b,'se':se,'CI_low':float(b-crit*se),'CI_high':float(b+crit*se),'p':p,'df':int(df),'rank':int(rank+1)}
br=[]
with threadpool_limits(limits=2):
 for subset in ['all','pruned']:
  for status in ['full','cases']:
   d=B.loc[(B.unrelated if subset=='pruned' else np.ones(len(B),dtype=bool))&(B.psoriasis_case.eq(1) if status=='cases' else True)]
   cols=['age','sex']+pcs+(['psoriasis_case'] if status=='full' else [])
   for scale in ['02','01']:
    family=[]
    for trait,ycol in traits.items():
     for _,v in M.iterrows():
      g=d[v.aligned_column].to_numpy(float).copy()
      if str(v.variant_id).startswith('chrX:') and scale=='01':g[d.sex.eq(1).to_numpy()]/=2
      r=ols(d[ycol].to_numpy(float),g,d[cols].to_numpy(float))
      family.append({'subset':subset,'status':status,'male_scale':scale,'trait':trait,'variant':v.variant_id,'EA':v.effect_allele,**r})
    assert len(family)==100
    p=np.array([r['p'] for r in family]);order=np.argsort(p);q=np.empty(100);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*100/np.arange(1,101))[::-1])[::-1])
    for r,qq in zip(family,q):r['q']=float(qq);r['bonferroni']=r['p']<.0005
    br+=family;pd.DataFrame(br).to_csv(R/'biomarker_results.private.tsv',sep=chr(9),index=False)
    print('BIO_DONE',subset,status,scale,'FDR',int((q<.05).sum()),'BONF',int((p<.0005).sum()),flush=True)
print('BIO_COMPLETE',len(br),flush=True)
