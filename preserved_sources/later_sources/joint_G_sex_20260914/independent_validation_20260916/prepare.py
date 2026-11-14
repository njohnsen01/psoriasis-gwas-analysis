from pathlib import Path
import json,pandas as pd,numpy as np,hashlib
p=Path(__file__).resolve().parent;b=p.parent;c=json.loads((b/'male02_pruned.command.json').read_text());v=lambda k:c[c.index(k)+1]
def read(f):
 d=pd.read_csv(f,sep=chr(92)+'s+',dtype=str);d.columns=d.columns.str.lstrip('#');assert not d.duplicated(['FID','IID']).any();return d
r=read(p/'targets_private.raw');cov=read(v('--covarFile'));ph=read(v('--phenoFile'));d=r.merge(cov,on=['FID','IID'],validate='one_to_one').merge(ph,on=['FID','IID'],validate='one_to_one');names=[x for x in cov.columns if x not in ['FID','IID']]
for x in names+['SEX','psoriasis_case']:d[x]=pd.to_numeric(d[x],errors='coerce')
d=d[d[names+['SEX','psoriasis_case']].notna().all(axis=1)].copy();assert set(d.sex.unique())=={0,1};assert (d.SEX==np.where(d.sex==1,1,2)).all();assert set(d.psoriasis_case.unique())=={0,1}
meta={};targets=Path(v('--extract')).read_text().split()
for i,t in enumerate(targets):
 col=t+'_'+t.rsplit(':',1)[1];assert col in d.columns;g=pd.to_numeric(d[col],errors='coerce');assert g.dropna().between(0,2).all();d['g'+str(i)]=g;meta['g'+str(i)]={'variant':t,'counted_allele':t.rsplit(':',1)[1],'input_scale':'0-2 both sexes; male divided by two for native carrier scale'}
for x in names:
 if x!='sex':
  sd=d[x].std();assert sd>0;d[x]=(d[x]-d[x].mean())/sd
out=d[['psoriasis_case']+names+list(meta)];out.to_csv(p/'model_data.private.csv',index=False)
(p/'MODEL_SPEC.json').write_text(json.dumps({'covariates':names,'targets':meta,'model':'Full nuisance-covariate refit on provider-pruned subset; no LOCO offset in this independent sensitivity, not exact REGENIE approximation replication','missingness':'variant-specific complete cases; same participants for each encoding','pruning_limitation':'Prior audit excludes released kinship pairs >0.1; distant/unreported relationships not ruled out'},indent=2));print('PREPARATION_COMPLETE: allele orientation and sex mapping verified; participant data retained internally',flush=True)
