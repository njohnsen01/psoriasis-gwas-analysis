from pathlib import Path
import json,hashlib,numpy as np,pandas as pd
p=Path(__file__).resolve().parent.parent;o=p/'deeper_analysis_20260919';v=p/'joint_G_sex_20260914/independent_validation_20260916';m=json.loads((v/'INPUT_MANIFEST.private.json').read_text())
def read(f):
 d=pd.read_csv(f,sep=r'\s+',dtype=str);d.columns=d.columns.str.lstrip('#');return d
raw=read(v/'targets_private.raw');cv=read(m['--covarFile']['path']);ph=read(m['--phenoFile']['path']);names=[x for x in cv if x not in ['FID','IID']]
d=raw.merge(cv,on=['FID','IID'],validate='one_to_one').merge(ph,on=['FID','IID'],validate='one_to_one')
for x in names+['SEX','psoriasis_case']:d[x]=pd.to_numeric(d[x],errors='coerce')
d=d.dropna(subset=names+['SEX','psoriasis_case']).copy();assert np.array_equal(d.SEX.to_numpy(),np.where(d.sex==1,1,2))
targets=Path(m['--extract']['path']).read_text().splitlines()
for i,t in enumerate(targets):d['g'+str(i)]=pd.to_numeric(d[t+'_'+t.split(':')[-1]],errors='coerce')
for x in names:
 if x!='sex':d[x]=(d[x]-d[x].mean())/d[x].std()
cols=['psoriasis_case']+names+['g0','g1','g2'];actual=pd.read_csv(v/'verified_irls_20260918/model_data.private.csv');assert list(actual.columns)==cols and len(actual)==len(d);assert np.allclose(d[cols].to_numpy(float),actual.to_numpy(float),rtol=1e-12,atol=1e-12,equal_nan=True)
a=p.parent.parent/'sample_lead_qc_20260910/sample_audit.private.tsv';meta=pd.read_csv(a,sep='\t',dtype={'IID':str});assert meta.IID.is_unique
z=d[['IID']+cols].merge(meta[['IID','ancestry_pred','site_id']],on='IID',how='left',validate='one_to_one');assert z.ancestry_pred.notna().all();z['ancestry']=z.ancestry_pred.str.upper();z['ancestry4']=z.ancestry.where(z.ancestry.isin(['AFR','AMR','EUR']),'OTHER');z.to_csv(o/'annotated_model_data.private.csv',index=False)
receipt={'status':'exact_row_order_genotypes_phenotype_and_scaled_covariates_reproduced','all_ancestry_labels_available':True,'model_input_sha256':hashlib.sha256((v/'verified_irls_20260918/model_data.private.csv').read_bytes()).hexdigest(),'metadata_sha256':hashlib.sha256(a.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256((o/'annotated_model_data.private.csv').read_bytes()).hexdigest(),'privacy':'Identifiers and participant-level data retained in Workbench only'}
(o/'ANCESTRY_INPUT_VERIFIED.json').write_text(json.dumps(receipt,indent=2));print('ANCESTRY_INPUT_VERIFIED')
