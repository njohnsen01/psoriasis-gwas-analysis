from pathlib import Path
import json,hashlib,traceback
import pandas as pd,numpy as np
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');Q=N.parent/'sample_lead_qc_20260910';I=N/'reconciliation_20260912';O=I/'read_mapping_QC';R=Path('/home/jupyter/psoriasis_reproducibility_20260910')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
try:
 d=pd.read_csv(Q/'sample_audit.private.tsv',sep=chr(9),dtype={'IID':str});d['gindex']=np.arange(len(d));g=np.load(Q/'lead_G.private.npy');h=np.load(Q/'lead_H.private.npy');D=pd.read_pickle(R/'reconstructed_x_v2/X_inputs.private.pkl');D.index=D.index.astype(str);s=d.loc[d.strict_eligible].copy();assert set(s.IID)==set(D.index);s=s.set_index('IID').loc[D.index];assert np.array_equal(s.psoriasis_case,D.psoriasis_case)
 E=pd.read_csv(Q/'four_sites_entries.private.tsv.bgz',sep=chr(9),compression='gzip',dtype={'s':str,'GT':str});E=E.loc[E.s.isin(s.index)].merge(s.reset_index(names='IID'),left_on='s',right_on='IID',validate='many_to_one')
 leads=[('rs183018930',5998724,'G'),('rs5970952',22917619,'T'),('rs5964607',67736589,'C'),('rs5919457',67842269,'T')];summaries=[];validations=[]
 for j,(rs,pos,ea) in enumerate(leads):
  expected=g[j,s.gindex.to_numpy()].astype(float);invalid=(expected<0)|((s.SEX.to_numpy()==1)&h[j,s.gindex.to_numpy()]);expected[invalid]=np.nan;assert np.allclose(expected,D[rs].to_numpy(),rtol=0,atol=0,equal_nan=True),rs+' current vector'
  e=E.loc[E.locus.eq('chrX:'+str(pos))].copy();alleles=json.loads(e.alleles.iloc[0]);idx=alleles.index(ea);e['dosage']=g[j,e.gindex.to_numpy()];valid=(e.dosage>=0)&~(e.SEX.eq(1)&h[j,e.gindex.to_numpy()]);e=e.loc[valid].copy()
  calls=e.GT.str.replace('|','/',regex=False).str.split('/').map(lambda x:list(map(int,x)));gt=calls.map(lambda x:sum(a==idx for a in x)*(2 if len(x)==1 else 1));assert np.array_equal(gt,e.dosage)
  def depths(value):
   if not isinstance(value,str):return (np.nan,np.nan,np.nan)
   a=json.loads(value)
   if a is None or len(a)!=len(alleles) or any(v is None or v<0 for v in a):return (np.nan,np.nan,np.nan)
   total=sum(a);return (total,a[idx],a[idx]/total if total>0 else np.nan)
  ad=e.AD.map(depths);e[['ADsum','AD_target','target_balance']]=pd.DataFrame(ad.tolist(),index=e.index);e['genotype_category']=np.where(e.SEX.eq(1),np.where(e.dosage.eq(2),'male_target','male_other'),np.where(e.dosage.eq(0),'female_no_target',np.where(e.dosage.eq(1),'female_target_heterozygote','female_target_homozygote')))
  for fields,family in [(['ancestry_pred','SEX','psoriasis_case','genotype_category'],'ancestry_sex_status_genotype'),(['site_id','SEX','genotype_category'],'center_sex_genotype')]:
   for key,z in e.groupby(fields,dropna=False):
    row=dict(lead=rs,family=family,**dict(zip(fields,key)),called_n=len(z))
    for metric in ['GQ','RGQ','ADsum','AD_target','target_balance']:
     v=pd.to_numeric(z[metric],errors='coerce');v=v.loc[np.isfinite(v)];row[metric+'_available_n']=len(v)
     for name,value in zip(['q25','median','q75'],v.quantile([.25,.5,.75]).tolist()):row[metric+'_'+name]=value
    row['GQ_below20_n']=int(pd.to_numeric(z.GQ,errors='coerce').lt(20).sum());row['ADsum_below10_n']=int(z.ADsum.lt(10).sum());row['AD_zero_n']=int(z.ADsum.eq(0).sum());row['FT_counts']=json.dumps(z.FT.fillna('missing').value_counts().to_dict());summaries.append(row)
  validations.append(dict(lead=rs,primary_dosage_match=True,source_GT_match=True,alleles=alleles,effect_allele=ea,AD_is_complete_DP=False))
 pd.DataFrame(summaries).to_csv(O/'genotype_quality_by_genotype.private.tsv',sep=chr(9),index=False)
 cm=pd.read_csv(O/'cram_manifest.private.csv',dtype={'person_id':str});eligible=s.reset_index(names='IID').merge(cm,left_on='IID',right_on='person_id',validate='one_to_one');eligible['NLGN4X_dosage']=D.rs183018930.reindex(eligible.IID).to_numpy();eligible=eligible.loc[eligible.NLGN4X_dosage.notna()];eligible['carrier']=eligible.NLGN4X_dosage.gt(0);eligible['selection_hash']=eligible.IID.map(lambda x:hashlib.sha256(('20260912:'+x).encode()).hexdigest());selected=[];counts=[]
 for anc in ['EUR','AFR','AMR']:
  for sex in [1,2]:
   for y in [0,1]:
    for carrier in [False,True]:
     a=eligible.loc[eligible.ancestry_pred.str.upper().eq(anc)&eligible.SEX.eq(sex)&eligible.psoriasis_case.eq(y)&eligible.carrier.eq(carrier)].sort_values('selection_hash');selected.append(a.head(2));counts.append(dict(ancestry=anc,sex=sex,case_status=y,carrier=carrier,eligible=len(a),selected=min(2,len(a))))
 selected=pd.concat(selected,ignore_index=True);selected['pilot_id']=['p'+str(i+1).zfill(3) for i in range(len(selected))];selected.to_csv(O/'CRAM_pilot_selection.private.tsv',sep=chr(9),index=False);pd.DataFrame(counts).to_csv(O/'CRAM_selection_strata.private.tsv',sep=chr(9),index=False)
 (O/'genotype_quality_COMPLETE.private.json').write_text(json.dumps(dict(status='verified',lead_validations=validations,rows=len(summaries),source_sha256=sha(Q/'four_sites_entries.private.tsv.bgz'),selection_n=len(selected),selection='Two available samples per ancestry/sex/case-control/carrier stratum, deterministic SHA256 ordering with seed 20260912. Pilot only; not cohort-wide read QC.',all_new_results_held=True),indent=2));print('GENOTYPE_QUALITY_AND_PILOT_SELECTION_COMPLETE',flush=True)
except Exception as e:(O/'genotype_quality_FAILED.private.json').write_text(json.dumps(dict(error=repr(e),traceback=traceback.format_exc()),indent=2));raise
