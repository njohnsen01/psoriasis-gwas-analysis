from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import mannwhitneyu,fisher_exact
O=Path(__file__).resolve().parent;assert (O/'MODEL_INPUT_VERIFIED.json').exists();d=pd.read_csv(O/'model_data.private.csv',low_memory=False);rows=[];private=[]
for t in ['g0','g1','g2']:
 for ancestry in ['AFR','AMR','EUR']:
  for sex in [0,1]:
   z=d[(d.ancestry4==ancestry)&(d.sex==sex)].copy();tag=dict(target=t,ancestry=ancestry,sex='female' if sex==0 else 'male')
   for metric,flag in [('genotype_missing',z[t].isna()),('provider_quality_failure',~z[t+'_quality_pass'].astype(str).str.lower().isin(['true','1']))]:
    zz=z.copy();zz['flag']=flag;cells=pd.crosstab(zz.psoriasis_case,zz.flag).reindex(index=[0,1],columns=[False,True],fill_value=0).to_numpy();private.append(dict(**tag,metric=metric,cells=cells.tolist()));r=dict(**tag,genotype='all',metric=metric)
    if np.all(cells>=21):
     odds,p=fisher_exact(cells);r.update(status='tested',case_vs_control_odds=float(odds),p=float(p),interpretation='Unadjusted within ancestry/sex; genotype composition can affect quality-failure comparison')
    else:r.update(status='support_gate_not_met')
    rows.append(r)
   for g in [0,1,2]:
    if sex==1 and g==1:continue
    zz=z[z[t]==g].copy();zz['AB_deviation']=(zz[t+'_AD_fraction']-g/2).abs()
    for metric,col in [('GQ',t+'_GQ'),('AD_depth',t+'_AD_depth'),('AB_deviation','AB_deviation')]:
     x=zz.loc[zz.psoriasis_case==1,col].dropna().astype(float).to_numpy();y=zz.loc[zz.psoriasis_case==0,col].dropna().astype(float).to_numpy();r=dict(**tag,genotype=str(g),metric=metric)
     if min(len(x),len(y))>=21:
      res=mannwhitneyu(x,y,alternative='two-sided',method='asymptotic');r.update(status='tested',case_median=float(np.median(x)),control_median=float(np.median(y)),rank_biserial=float(2*res.statistic/(len(x)*len(y))-1),p=float(res.pvalue),interpretation='Within ancestry/sex/genotype descriptive comparison')
     else:r.update(status='support_gate_not_met')
     rows.append(r)
r=pd.DataFrame(rows);ix=r.index[r.status=='tested'];p=r.loc[ix,'p'].to_numpy();order=np.argsort(p);adj=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[order]=np.minimum(1,adj);r.loc[ix,'BH_q']=q
r.to_csv(O/'QUALITY_DIAGNOSTICS.csv',index=False);(O/'quality_cells.private.json').write_text(json.dumps(private));assert r.loc[ix,['p','BH_q']].apply(lambda s:s.between(0,1).all()).all()
(O/'QUALITY_DIAGNOSTICS_VERIFIED.json').write_text(json.dumps(dict(status='complete',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),all_released_tests_pass_support=True,multiple_testing='BH across all eligible provider-quality comparisons',limitations='Descriptive within-cohort checks. Quality-failure Fisher comparison is not genotype-adjusted; metric rank tests condition on genotype. No full-cohort read-level MQ exists in this table.'),indent=2));print('QUALITY_DIAGNOSTICS_VERIFIED')
