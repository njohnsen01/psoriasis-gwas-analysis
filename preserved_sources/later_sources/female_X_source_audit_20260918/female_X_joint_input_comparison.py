from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,subprocess
P=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');B=P/'reconciliation_20260912';J=B/'joint_G_sex_20260914';V=J/'independent_validation_20260916';O=B/'female_X_source_audit_20260918'
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as z:
  for block in iter(lambda:z.read(8388608),b''): h.update(block)
 return h.hexdigest()
m=json.loads((V/'INPUT_MANIFEST.private.json').read_text());checks={}
for k,v in m.items(): checks['manifest_'+k]=sha(v['path'])==v['sha256']
def read(f):
 d=pd.read_csv(f,sep=r'\s+',dtype=str);d.columns=d.columns.str.lstrip('#');return d
raw=read(V/'targets_private.raw');cv=read(m['--covarFile']['path']);ph=read(m['--phenoFile']['path']);keep=read(m['--keep']['path']) if False else None
names=[x for x in cv if x not in ['FID','IID']]
d=raw.merge(cv,on=['FID','IID'],validate='one_to_one').merge(ph,on=['FID','IID'],validate='one_to_one')
for x in names+['SEX','psoriasis_case']:d[x]=pd.to_numeric(d[x],errors='coerce')
d=d[d[names+['SEX','psoriasis_case']].notna().all(axis=1)].copy()
checks['sex_labels_consistent']=bool((d.SEX==np.where(d.sex==1,1,2)).all())
source=P/'X_projected_source/chunk_0009';ps=read(source/'input.psam');pv=read(source/'input.pvar')
ix=ps.set_index('IID').index.get_indexer(d.IID);assert (ix>=0).all()
cmd=json.loads((V/'plink.command.json').read_text());cmd[cmd.index('--out')+1]=str(O/'current_targets.private')
assert not (O/'current_targets.private.raw').exists()
with (O/'current_targets.private.console.log').open('w') as h: rc=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT).returncode
assert rc==0,'Private PLINK diagnostic failed'
fresh=read(O/'current_targets.private.raw').set_index('IID').loc[d.IID]
targets=Path(m['--extract']['path']).read_text().split();meta={}
for i,t in enumerate(targets):
 assert int((pv.ID==t).sum())==1
 col=t+'_'+t.split(':')[-1];d['g'+str(i)]=pd.to_numeric(d[col],errors='coerce');a=pd.to_numeric(fresh[col],errors='coerce').to_numpy(float)
 checks['genotype_matches_source_'+str(i)]=bool(np.allclose(a,d['g'+str(i)],equal_nan=True))
 meta['g'+str(i)]={'variant':t,'effect_allele':t.split(':')[-1]}

female_ids=set(read(P/'stratified_step1/female/phenotype.tsv').IID)
checks['joint_females_subset_scan']=set(d.loc[d.sex.eq(0),'IID']).issubset(female_ids)
for x in names:
 if x!='sex': d[x]=(d[x]-d[x].mean())/d[x].std()
expected=d[['psoriasis_case']+names+list(meta)];actual=pd.read_csv(V/'verified_irls_20260918/model_data.private.csv')
checks['reconstructed_joint_model_data_identical']=expected.shape==actual.shape and list(expected)==list(actual) and np.allclose(expected.to_numpy(float),actual.to_numpy(float),rtol=1e-12,atol=1e-12,equal_nan=True)
assert all(checks.values()),'Comparison failed; retain private diagnostics'
report={'checks':checks,'all_passed':True,'covariates_joint':names,'genotype_source':'same corrected provider PGEN chunk_0009 as female-only scan','joint_cohort':'provider-pruned subset of corrected primary; all joint-model women belong to the female-only scan','scaling':'same ALT dosage; canonical female per-allele and male carrier contrasts','differences':['participant subset','female-only refitted LOCO versus no LOCO','female-only nuisance fitting versus pooled nuisance fitting','REGENIE approximate Firth versus full-covariate Firth PLRT'],'conclusion':'Same genotype source and allele orientation verified; model differences are not isolated by this comparison.'}
(O/'JOINT_INPUT_COMPARISON.json').write_text(json.dumps(report,indent=2));print('JOINT_INPUT_COMPARISON_ALL_CHECKS_PASSED');print(json.dumps(report,indent=2))

