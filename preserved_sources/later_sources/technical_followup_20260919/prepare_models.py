from pathlib import Path
import json,hashlib,collections,re,shutil
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;B=O.parent;N=B.parent;D=B/'deeper_analysis_20260919';S=N/'stratified_step1/female';V=B/'joint_G_sex_20260914/independent_validation_20260916'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
 x=pd.read_csv(p,sep=chr(92)+'s+',dtype=str);x.columns=x.columns.str.lstrip('#');return x
assert (O/'PROTOCOL.json').exists() and not (O/'MODEL_INPUT_VERIFIED.json').exists()
d=pd.read_csv(D/'annotated_model_data.private.csv',dtype={'IID':str});original=d.copy();assert d.IID.is_unique
meta=json.loads((S/'VALIDATED.json').read_text());assert sha(S/'model_1.loco')==meta['loco_sha256'];assert sha(S/'phenotype.tsv')==meta['phenotype_sha256'];assert sha(S/'covariates.tsv')==meta['covariate_sha256']
ph=read(S/'phenotype.tsv');assert ph.IID.is_unique;ids=(ph.FID+'_'+ph.IID).tolist();loco={}
with (S/'model_1.loco').open() as f:
 head=f.readline().split()[1:];assert len(head)==len(set(head)) and set(head)==set(ids)
 for line in f:
  z=line.split();values=np.array(z[1:],float);assert len(values)==len(head) and np.isfinite(values).all();loco[int(z[0])]=values
assert set(loco)==set(range(1,24));mapper=dict(zip(head,loco[23]));ph['loco']=[mapper[x] for x in ids]
female=d[d.sex==0];assert set(female.IID)<=set(ph.IID);pp=ph.set_index('IID');assert np.array_equal(female.psoriasis_case.to_numpy(),pd.to_numeric(pp.loc[female.IID,'psoriasis_case']).to_numpy())
d=d.merge(ph[['IID','FID','loco']],on='IID',how='left',validate='one_to_one');assert d.loc[d.sex==0,'loco'].notna().all()
q=read(N.parent/'sample_lead_qc_20260910/sample_audit.private.tsv');assert q.IID.is_unique
qcols=['mean_coverage','genome_coverage','dragen_contamination','verify_bam_id2_contamination'];d=d.merge(q[['IID']+qcols],on='IID',how='left',validate='one_to_one')
for c in qcols:d[c]=pd.to_numeric(d[c],errors='coerce')
d['site_id']=d.site_id.fillna('MISSING').astype(str);f=d[d.sex==0];tab=pd.crosstab(f.site_id,f.psoriasis_case).reindex(columns=[0,1],fill_value=0);eligible=tab.index[(tab>=21).all(axis=1)].tolist();labels={s:'S'+str(i+1).zfill(3) for i,s in enumerate(sorted(eligible))};d['site_cat']=d.site_id.map(labels).fillna('OTHER')
(O/'site_mapping.private.json').write_text(json.dumps(dict(mapping=labels,rule='Original site case/control counts each>=21 in pruned women; else OTHER. LOSO retained cohort and omitted site must each meet case/carrier >=21 per cell for reported target.'),indent=2))
valid=set(d.IID);records=[]
targets=(B/'joint_G_sex_20260914/targets.txt').read_text().splitlines()
for e in pd.read_csv(B/'joint_G_sex_20260914/provider_quality/three_sites_entries.private.tsv.bgz',sep=chr(9),compression='gzip',dtype={'s':str,'GT':str},chunksize=100000):
 e=e[e.s.isin(valid)]
 for row in e.itertuples():
  matches=[i for i,t in enumerate(targets) if row.locus==':'.join(t.split(':')[:2])];assert len(matches)==1;j=matches[0];alleles=json.loads(row.alleles);alt=targets[j].split(':')[-1];aidx=alleles.index(alt);gt=[int(x) for x in re.split(r'[/|]',row.GT)] if isinstance(row.GT,str) and re.fullmatch(chr(92)+'d+([/|]'+chr(92)+'d+)?',row.GT) else None
  try:ad=json.loads(row.AD) if isinstance(row.AD,str) else None
  except Exception:ad=None
  adok=isinstance(ad,list) and len(ad)==len(alleles) and all(x is not None and x>=0 for x in ad);depth=sum(ad) if adok else np.nan;frac=ad[aidx]/depth if adok and depth>0 else np.nan
  records.append(dict(IID=row.s,target='g'+str(j),GQ=pd.to_numeric(row.GQ,errors='coerce'),AD_depth=depth,AD_fraction=frac,provider_called=gt is not None,provider_fraction=sum(x==aidx for x in gt)/len(gt) if gt else np.nan,provider_biallelic=bool(gt) and all(x in [0,aidx] for x in gt)))
e=pd.DataFrame(records);assert not e.duplicated(['IID','target']).any();e.to_csv(O/'provider_quality.private.csv',index=False)
for j in range(3):
 t='g'+str(j);z=e[e.target==t].drop(columns='target').set_index('IID').add_prefix(t+'_');d=d.join(z,on='IID',validate='one_to_one')
 expected=d[t]/2;got=d[t+'_provider_fraction'];called=d[t].notna();assert got[called].notna().all() and np.allclose(got[called],expected[called],rtol=0,atol=1e-12)
 gq=d[t+'_GQ'];dep=d[t+'_AD_depth'];af=d[t+'_AD_fraction'];balance=np.where(expected.eq(0),af.le(.1),np.where(expected.eq(1),af.ge(.9),af.between(.25,.75)))
 d[t+'_quality_pass']=called & gq.ge(20) & dep.ge(10) & d[t+'_provider_biallelic'].eq(True) & balance
for c in original.columns:
 if c!='site_id':assert original[c].fillna('NA').equals(d[c].fillna('NA'))
d.to_csv(O/'model_data.private.csv',index=False)
d.loc[d.sex==0,['FID','IID']].to_csv(O/'pruned_female.keep.private.tsv',sep=chr(9),index=False,header=False)
receipt=dict(status='verified',prior_input_unchanged=sha(D/'annotated_model_data.private.csv')==json.loads((O/'PROTOCOL.json').read_text())['prior_input_sha256'],all_female_ids_and_phenotype_verified=True,loco_hash_verified=True,loco_chromosome=23,loco_scale='REGENIE binary-trait log-odds linear predictor, fixed offset',provider_target_fractions_match_accepted_calls=True,rare_site_grouping_frozen_before_fits=True,sequencing_batch_field_available=False,available_technical_covariates=qcols,source_sha256=sha(__file__),output_sha256=sha(O/'model_data.private.csv'))
assert receipt['prior_input_unchanged'];(O/'MODEL_INPUT_VERIFIED.json').write_text(json.dumps(receipt,indent=2));print('MODEL_INPUT_VERIFIED',flush=True)
