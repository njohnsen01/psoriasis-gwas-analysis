from pathlib import Path
import json,hashlib,datetime,traceback,importlib.metadata
import numpy as np,pandas as pd
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');R=Path('/home/jupyter/psoriasis_reproducibility_20260910');I=N/'reconciliation_20260912';O=I/'followup_input_audit';O.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,v):(O/name).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def main():
 sources={};checks={};counts=[]
 def remember(p):sources[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
 def read(p):
  remember(p);d=pd.read_csv(p,sep=chr(9),dtype={'FID':str,'IID':str}).set_index('IID');assert d.index.is_unique;return d
 ph=read(N/'inputs/phenotype.tsv');cv=read(N/'inputs/covariates.tsv');ds={}
 for name,folder,file in [('X','reconstructed_x_v2','X_inputs.private.pkl'),('biomarker','reconstructed_biomarkers_v2','biomarker_inputs.private.pkl'),('ancestry','reconstructed_ancestry_v2','ancestry_autosomal_inputs.private.pkl'),('SII','reconstructed_sii_v2','biomarker_inputs.private.pkl')]:
  p=R/folder/file;remember(p);d=pd.read_pickle(p);d.index=d.index.astype(str);assert d.index.is_unique and set(d.index)==set(ph.index),(name,'participant mismatch')
  d=d.reindex(ph.index);pd.testing.assert_series_equal(d.psoriasis_case,ph.psoriasis_case,check_names=False,check_dtype=False,check_exact=True);diff={}
  for col in cv.columns.drop('FID'):
   a=d[col].to_numpy(float);b=cv[col].reindex(d.index).to_numpy(float);assert np.isfinite(a).all() and np.isfinite(b).all();delta=float(np.max(np.abs(a-b)));diff[col]=delta
   assert np.allclose(a,b,rtol=1e-12,atol=1e-12),(name,col,delta)
  checks[name]={'same_participant_keys':True,'same_phenotypes':True,'covariate_max_abs_differences':diff,'float_parse_tolerance':'rtol=1e-12,atol=1e-12'};ds[name]=d
  counts.append({'input':name,'n':len(d),'cases':int(d.psoriasis_case.sum()),'controls':int(d.psoriasis_case.eq(0).sum())})
 D=ds['X'];B=ds['biomarker'];G=ds['ancestry'];S=ds['SII'];M=pd.read_csv(R/'reconstructed_inputs_v2/variant_map.tsv',sep=chr(9));remember(R/'reconstructed_inputs_v2/variant_map.tsv')
 assert D.SEX.eq(np.where(D.sex.eq(1),1,2)).all();checks['SEX_matches_primary_sex_covariate']=True
 ap=N.parent/'phenotype_validation_20260910/private_full_ancestry.tsv';remember(ap);anc=pd.read_csv(ap,sep=chr(9),dtype={'research_id':str},usecols=['research_id','ancestry_pred']).set_index('research_id');assert anc.index.is_unique
 for name in ['X','ancestry']:
  a=ds[name].ancestry_pred.astype(str).str.upper();b=anc.ancestry_pred.reindex(a.index).astype(str).str.upper();assert a.eq(b).all(),(name,'ancestry mismatch')
 checks['complete_provider_ancestry_matches']=True
 pp=N.parent/'relatedness_sensitivity_20260910/samples_relatedness_flagged_samples.tsv';remember(pp);flag=pd.read_csv(pp,sep=chr(9),dtype=str).sample_id;expected=~D.index.isin(flag)
 for name,d in ds.items():assert d.unrelated.to_numpy(bool).tolist()==expected.tolist(),(name,'pruning mismatch')
 pairp=N.parent/'relatedness_sensitivity_20260910/samples_relatedness.tsv';remember(pairp);pairs=pd.read_csv(pairp,sep=chr(9),dtype={'i.s':str,'j.s':str});left=set(D.index[expected]);ret=pairs['i.s'].isin(left)&pairs['j.s'].isin(left)&pairs.kin.gt(.1);assert not ret.any()
 checks['pruning']={'matches_fixed_provider_list':True,'no_retained_released_pair_above_0_1':True,'lower_relatedness_not_excluded_by_this_claim':True}
 leads={'rs183018930':'chrX:5998724:A:G','rs5970952':'chrX:22917619:C:T','rs5964607':'chrX:67736589:T:C','rs5919457':'chrX:67842269:G:T'};geno=[]
 for rs,vid in leads.items():
  v=M.loc[M.variant_id.eq(vid)];assert len(v)==1;col=v.iloc[0].aligned_column
  for name in ['biomarker','ancestry','SII']:pd.testing.assert_series_equal(D[rs],ds[name][col],check_names=False,check_dtype=False,check_exact=True)
  assert D[rs].dropna().isin([0,1,2]).all() and D.loc[D.SEX.eq(1),rs].dropna().isin([0,2]).all()
  for subset,z in [('all',D),('provider_pruned',D.loc[D.unrelated])]:
   for group in ['pooled','female','male','EUR','AFR','AMR']:
    a=z if group=='pooled' else z.loc[z.SEX.eq(2 if group=='female' else 1)] if group in ['female','male'] else z.loc[z.ancestry_pred.str.upper().eq(group)]
    for case,t in a.groupby('psoriasis_case'):
     g=t[rs];valid=g.notna();male=t.SEX.eq(1);ac=float((g/np.where(male,2.,1.)).sum());an=int(np.where(male[valid],1,2).sum())
     geno.append({'variant':vid,'rsID':rs,'subset':subset,'group':group,'case':int(case),'eligible_n':len(t),'called_n':int(valid.sum()),'missing_n':int((~valid).sum()),'AC':ac,'AN':an,'MAC':min(ac,an-ac),'call_rate':float(valid.mean()),'dosage':'female0/1/2 male0/2; biological male AC dosage/2'})
 checks['four_X_vectors_identical_across_followups']=True
 cols=[c for c in B.columns if c in S.columns and c not in ['FID']]
 pd.testing.assert_frame_equal(B[cols],S[cols],check_dtype=False,check_exact=True);checks['SII_original_base_equals_biomarker_base']=True
 pd.testing.assert_frame_equal(B[cols],G[cols],check_dtype=False,check_exact=True);checks['ancestry_base_equals_biomarker_base']=True
 pd.DataFrame(geno).to_csv(O/'X_biological_counts.private.tsv',sep=chr(9),index=False)
 pd.DataFrame(counts).to_csv(O/'input_cohorts.private.tsv',sep=chr(9),index=False)
 keep=M[['physical_group','variant_id','ref','effect_allele','source_counted_allele','alignment_action','aligned_column']].copy();keep['selection']='historical pooled lead panel; exploratory retained targets';keep['denominator_source']='model-specific complete cases; legacy map denominator columns are not current';keep.to_csv(O/'historical_target_provenance.private.tsv',sep=chr(9),index=False)
 save('input_files.private.json',sources);save('COMPLETE.private.json',{'status':'input_equality_and_provider_pruning_verified','checks':checks,'release':'HELD','model_estimate_reconciliation':'separate pending step','versions':{x:importlib.metadata.version(x) for x in ['numpy','pandas']},'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});print('FOLLOWUP_INPUT_AUDIT_COMPLETE',flush=True)
try:main()
except Exception as e:
 save('FAILED.private.json',{'error':repr(e),'traceback':traceback.format_exc(),'release':'HELD'});raise
