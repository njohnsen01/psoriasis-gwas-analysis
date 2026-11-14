import os
from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import chi2
A=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09')
F=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}/final_aggregate')
p=pd.read_csv(F/'psoriasis_corrected_manhattan_plot_data.tsv.gz',sep='	');v=p.LOG10P.to_numpy();g=v>=-np.log10(5e-8);n=13255796
w=np.where(g,1.,(n-g.sum())/(~g).sum());o=np.argsort(-v);sv=v[o];cw=np.cumsum(w[o]);exact=np.load(A/'all_sorted_logp.npy',mmap_mode='r')
maxdiff=0.;sumsq=0.;sumabs=0.;maxrank=0
for start in range(1,n+1,1000000):
 r=np.arange(start,min(start+1000000,n+1));est=sv[np.minimum(np.searchsorted(cw,r,side='left'),len(sv)-1)];dif=est-exact[r-1];j=np.argmax(np.abs(dif))
 if abs(dif[j])>maxdiff:maxdiff=float(abs(dif[j]));maxrank=int(r[j])
 sumsq+=float(np.square(dif).sum());sumabs+=float(np.abs(dif).sum())
# Reproduce the exact ranks used by the original Figure 1b builder.
r=np.unique(np.rint(np.r_[np.arange(1,2501),np.geomspace(2501,n,8500)]).astype(int));r=r[(r>=1)&(r<=n)]
est=sv[np.minimum(np.searchsorted(cw,r),len(sv)-1)];dif=est-exact[r-1]
result={'scope':'all full-output ranks and original displayed ranks','n_full_ranks':n,'n_displayed_ranks':len(r),'full_rank_max_abs_logp_difference':maxdiff,'full_rank_max_difference_rank':maxrank,'full_rank_mean_abs_logp_difference':sumabs/n,'full_rank_rmse':(sumsq/n)**.5,'displayed_max_abs_logp_difference':float(np.abs(dif).max()),'displayed_mean_abs_logp_difference':float(np.abs(dif).mean()),'background_weight':float(w[~g][0]),'significant_weight':1,'all_GWS_order_statistics_exact':bool(np.array_equal(sv[:int(g.sum())],exact[:int(g.sum())]))}
(A/'weighted_plot_full_rank_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
# Audit the complete pooled X output by PAR region and allele frequency.
x=[]
for d in pd.read_csv(F/'psoriasis_corrected_regenie_all_results.tsv.gz',sep='	',chunksize=500000,dtype={'CHROM':str}):
 z=d[d.CHROM.isin(['X','23'])];x.append(z)
x=pd.concat(x,ignore_index=True);assert len(x)==376171
lp=x.LOG10P.to_numpy();pos=x.GENPOS if 'GENPOS' in x else x.POS;par=(pos<=2781479)|(pos>=155701383);af=x.A1FREQ if 'A1FREQ' in x else x.EAF;maf=np.minimum(af,1-af)
rows=[]
for key,mask in [('X_all',np.ones(len(x),dtype=bool)),('X_nonPAR',~par),('X_PAR',par),('X_MAF_lt_0.01',maf<.01),('X_MAF_0.01_to_0.05',(maf>=.01)&(maf<.05)),('X_MAF_ge_0.05',maf>=.05)]:
 q=lp[np.asarray(mask)];rows.append({'group':key,'n_tests':len(q),'lambda_GC':float(np.median(chi2.isf(np.power(10.,-q),1))/chi2.ppf(.5,1)) if len(q) else None})
pd.DataFrame(rows).to_csv(A/'pooled_X_additional_diagnostics.tsv',sep='	',index=False);print(pd.DataFrame(rows).to_string(index=False))
