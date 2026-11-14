from pathlib import Path
import sys
sys.path.insert(0,'/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/X_source_reanalysis_20260910/python_dependencies')
import pandas as pd,numpy as np,pgenlib,json
from scipy.stats import fisher_exact
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');O=B/'reconciliation_20260912/joint_G_sex_20260914';P=B/'X_projected_source/chunk_0009/input'
ps=pd.read_csv(str(P)+'.psam',sep=r'\s+',dtype={'IID':str});pv=pd.read_csv(str(P)+'.pvar',sep=r'\s+');cv=pd.read_csv(B/'inputs/covariates.tsv',sep=r'\s+',dtype={'IID':str});y=pd.read_csv(B/'inputs/phenotype.tsv',sep=r'\s+',dtype={'IID':str});D=pd.read_pickle('/home/jupyter/psoriasis_reproducibility_20260910/reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');D.index=D.index.astype(str)
z=ps[['IID','SEX']].merge(cv[['IID','sex']],on='IID',validate='one_to_one').merge(y[['IID','psoriasis_case']],on='IID',validate='one_to_one');z['ancestry']=D.ancestry_pred.reindex(z.IID).to_numpy();assert z.ancestry.notna().all();idx=pd.Series(np.arange(len(ps)),index=ps.IID).loc[z.IID].to_numpy();rows=[]
r=pgenlib.PgenReader(bytes(str(P)+'.pgen','utf8'))
for ID in (O/'targets.txt').read_text().split():
 g=np.empty(len(ps),dtype=np.int8);r.read(int(pv.index[pv.ID.eq(ID)][0]),g);z['g']=g[idx];z['valid']=z.g.ge(0)&~(z.SEX.eq(1)&z.g.eq(1))
 for (a,s),d in z.groupby(['ancestry','SEX']):
  t=pd.crosstab(d.psoriasis_case,d.valid).reindex(index=[0,1],columns=[False,True],fill_value=0);p=float(fisher_exact(t.to_numpy()).pvalue)
  for status,e in d.groupby('psoriasis_case'):
   good=e[e.valid];ac=float(good.g.sum()/(2 if s==1 else 1));rows.append(dict(variant=ID,ancestry=a,sex=int(s),case_status=int(status),N=len(e),called=len(good),call_rate=len(good)/len(e),effect_allele_count=ac,differential_missingness_P=p))
r.close();pd.DataFrame(rows).to_csv(O/'target_QC_by_sex_ancestry.private.tsv',sep='\t',index=False);(O/'target_QC_COMPLETE.json').write_text(json.dumps({'complete':True,'targets':3,'scope':'call rate, effect allele counts and case-control differential missingness by sex and ancestry; source PGEN; all counts private','mapping_quality':'not established by PGEN','genotype_quality':'not established by PGEN'}));print('TARGET_QC_COMPLETE')
