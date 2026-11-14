from pathlib import Path
import json,hashlib,traceback
import numpy as np,pandas as pd
from scipy import stats
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');R=Path('/home/jupyter/psoriasis_reproducibility_20260910');I=N/'reconciliation_20260912';O=I/'followup_model_audit';O.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def main():
 assert json.loads((I/'followup_input_audit/COMPLETE.private.json').read_text())['status']=='input_equality_and_provider_pruning_verified'
 Q=pd.read_csv(I/'estimate_audit/target_QC.private.tsv.gz',sep=chr(9));D=pd.read_pickle(R/'reconstructed_x_v2/X_inputs.private.pkl');B=pd.read_pickle(R/'reconstructed_biomarkers_v2/biomarker_inputs.private.pkl');G=pd.read_pickle(R/'reconstructed_ancestry_v2/ancestry_autosomal_inputs.private.pkl');M=pd.read_csv(R/'reconstructed_inputs_v2/variant_map.tsv',sep=chr(9));S=pd.read_pickle(R/'reconstructed_sii_v2/SII_person.private.pkl');S.index=S.index.astype(str);sources={};audit=[];outputs={};pcs=['PC'+str(i) for i in range(1,11)]
 # Persist the already inspected biological-count comparison.
 C=pd.read_csv(I/'followup_input_audit/X_biological_counts.private.tsv',sep=chr(9));comparisons=0
 for row in Q.loc[Q.domain.eq('X')].itertuples():
  for case,label in [(1,'case'),(0,'control')]:
   t=C.loc[C.subset.eq('all')&C.group.eq(row.group)&C.variant.eq(row.ID)&C['case'].eq(case)].iloc[0]
   for a,b in [('AC','AC'),('AN','AN'),('called_n','called'),('call_rate','callrate')]:assert np.isclose(t[a],getattr(row,label+'_'+b),rtol=1e-12,atol=1e-12);comparisons+=1
 def status(vid):
  z=Q.loc[Q.group.eq('pooled')&Q.ID.eq(vid)]
  return 'passed' if len(z)==1 and bool(z.iloc[0].PASS) else 'failed' if len(z)==1 else 'not_reconciled_to_primary_panel'
 for family,folder,file in [('X','reconstructed_x_v2','X_results.private.tsv'),('ancestry','reconstructed_ancestry_v2','ancestry_autosomal_results.private.tsv'),('biomarker','reconstructed_biomarkers_v2','biomarker_results.private.tsv'),('SII','reconstructed_sii_v2','matched_results.private.tsv')]:
  p=R/folder/file;sources[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size};t=pd.read_csv(p,sep=chr(9),dtype={'male_scale':str});estimated=t.status.eq('estimated') if family in ['X','ancestry'] else np.ones(len(t),bool);z=t.loc[estimated]
  assert np.isfinite(z[['beta','se','p']]).all().all() and z.se.gt(0).all();assert z.p.between(0,1).all()
  if family in ['X','ancestry']:
   assert z.converged.all();assert np.allclose(z.OR,np.exp(z.beta),rtol=1e-10,atol=1e-12);assert np.allclose(z.CI_low,np.exp(z.beta-stats.norm.ppf(.975)*z.se),rtol=1e-10,atol=1e-12);assert np.allclose(z.CI_high,np.exp(z.beta+stats.norm.ppf(.975)*z.se),rtol=1e-10,atol=1e-12);assert np.allclose(z.p,2*stats.norm.sf(abs(z.beta/z.se)),rtol=1e-10,atol=1e-14)
   base=D if family=='X' else G
   for row in t.itertuples():
    a=base if row.subset=='all' else base.loc[base.unrelated];a=a if row.stratum=='pooled' else a.loc[a.SEX.eq(2 if row.stratum=='female' else 1)] if row.stratum in ['female','male'] else a.loc[a.ancestry_pred.str.upper().eq(row.stratum)]
    cov=['age','age2']+pcs if row.stratum in ['female','male'] else ['sex','age','age2','age_sex']+pcs
    names=[row.lead]+([row.condition] if row.condition!='none' else []);complete=a[cov+names+['psoriasis_case']].replace([np.inf,-np.inf],np.nan).dropna();assert len(complete)==row.n and int(complete.psoriasis_case.sum())==row.cases and len(complete)-row.cases==row.controls
   t['P_method']='two-sided Wald normal';t['CI_method']='exp(beta +/- 1.95996398454 SE)';t['relatedness']='all versus fixed provider-pruned; no retained released pair above kinship 0.1 in pruned'
  else:
   assert (z['n']-z['rank']).eq(z.df).all();crit=stats.t.ppf(.975,z.df);assert np.allclose(z.CI_low,z.beta-crit*z.se,rtol=1e-10,atol=1e-12);assert np.allclose(z.CI_high,z.beta+crit*z.se,rtol=1e-10,atol=1e-12);assert np.allclose(z.p,2*stats.t.sf(abs(z.beta/z.se),z.df),rtol=1e-10,atol=1e-14)
   famcols=['subset','status','male_scale'] if family=='biomarker' else ['unit_rules','subset','status','male_scale','outcome'];k=100 if family=='biomarker' else 25;qcol='q' if family=='biomarker' else 'q25'
   for _,a in t.groupby(famcols):
    assert len(a)==k;pv=a.p.to_numpy();ix=np.argsort(pv);qq=np.empty(k);qq[ix]=np.minimum(1,np.minimum.accumulate((pv[ix]*k/np.arange(1,k+1))[::-1])[::-1]);assert np.allclose(a[qcol],qq,rtol=1e-10,atol=1e-14)
   cols=['age','sex']+pcs;mapcols=M.set_index('variant_id').aligned_column.to_dict();outcomes={'SII':'log_inflammatory_index_z','HDL':'HDL_z','LDL':'LDL_z','TG':'log_TG_z'}
   for row in t.itertuples():
    a=B.loc[(B.unrelated if row.subset=='pruned' else np.ones(len(B),dtype=bool))&(B.psoriasis_case.eq(1) if row.status=='cases' else True)].copy();cov=cols+(['psoriasis_case'] if row.status=='full' else []);gcol=mapcols[row.variant]
    if family=='SII':
     vals=S.reindex(a.person_id.astype(str));ok=vals[row.unit_rules].notna().to_numpy()&vals.original_rebuilt.notna().to_numpy();a=a.loc[ok]
     n=int(np.isfinite(a[cov+[gcol]]).all(axis=1).sum())
    else:n=int(np.isfinite(a[cov+[gcol,outcomes[row.trait]]]).all(axis=1).sum())
    assert n==row.n,(family,row.variant,n,row.n)
   t['P_method']='two-sided Student t with residual df';t['CI_method']='beta +/- t(0.975,df) SE';t['relatedness']='OLS all versus fixed provider-pruned sensitivity'
  t['primary_QC_status']=t.variant.map(status);t['interpretation']=np.where(t.primary_QC_status.eq('passed'),'exploratory historical-target follow-up','diagnostic historical target; not a primary discovery');t['effect_unit']='per coded target-ALT dosage unit; X male 02 is two units per carrier contrast';t['source_output']=str(p);t['source_sha256']=sources[str(p)]['sha256'];t.to_csv(O/(family+'_reconciled.private.tsv'),sep=chr(9),index=False);outputs[family]=t;audit.append({'family':family,'records':len(t),'estimated_records':int(estimated.sum()),'same_output_arithmetic':True,'model_specific_denominators_reproduced':True})
 # Current exploratory ancestry heterogeneity, with the entire attempted target family retained.
 a=pd.concat([outputs['ancestry'],outputs['X'].loc[outputs['X'].stratum.isin(['EUR','AFR','AMR'])&outputs['X'].condition.eq('none')&outputs['X'].male_scale.eq('02')]],ignore_index=True);h=[]
 for subset,z in a.groupby('subset'):
  family_size=z.variant.nunique()
  for vid,v in z.groupby('variant'):
   z1=v.loc[v.status.eq('estimated')];rec={'subset':subset,'variant':vid,'attempted_family_size':int(family_size),'estimable_groups':len(z1),'primary_QC_status':status(vid),'scope':'exploratory; internal subsets; failed-QC targets diagnostic only'}
   if len(z1)<2:rec['status']='not_estimable'
   else:
    w=1/z1.se.to_numpy()**2;b=z1.beta.to_numpy();mean=float((w*b).sum()/w.sum());q=float((w*(b-mean)**2).sum());p=float(stats.chi2.sf(q,len(b)-1));rec.update(status='estimated',Q=q,df=len(b)-1,p=p,bonferroni_p=min(1.,p*family_size),nominal=p<.05,bonferroni=p<.05/family_size)
   h.append(rec)
 pd.DataFrame(h).to_csv(O/'current_ancestry_heterogeneity.private.tsv',sep=chr(9),index=False);save('input_files.private.json',sources);save('COMPLETE.private.json',{'status':'same-output_arithmetic_and_denominators_verified','families':audit,'X_count_comparisons':comparisons,'new_current_ancestry_family':'includes all attempted current targets; does not substitute for historical 24-test family','release':'HELD','script_sha256':sha(Path(__file__))});print('FOLLOWUP_MODEL_AUDIT_COMPLETE',flush=True)
try:main()
except Exception as e:save('FAILED.private.json',{'error':repr(e),'traceback':traceback.format_exc(),'release':'HELD'});raise
