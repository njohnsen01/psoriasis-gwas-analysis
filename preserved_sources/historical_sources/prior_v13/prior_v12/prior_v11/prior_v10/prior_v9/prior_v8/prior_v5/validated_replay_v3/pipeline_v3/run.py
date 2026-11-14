#!/usr/bin/env python3
"""Secure reconstruction and replay v3. Replacement genome-wide run is incomplete."""
from pathlib import Path
import argparse,contextlib,hashlib,importlib.metadata as im,json,os,runpy,shutil,subprocess,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('stage',choices=['check','cohort','labs','lab-processing','x','biomarkers','sii','ancestry'])
p.add_argument('--out',type=Path,required=True)
p.add_argument('--audit-root',type=Path,default=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09'))
p.add_argument('--input-root',type=Path,default=Path(f'{os.environ["PSO_WORKSPACE_MOUNT"].rstrip("/")}'))
p.add_argument('--cohort-dir',type=Path)
p.add_argument('--classified-labs',type=Path)
p.add_argument('--ancestry-file',type=Path)
p.add_argument('--model-input-dir',type=Path)
a=p.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=True)
if any(a.out.iterdir()):raise RuntimeError('Output directory must be empty; use a fresh directory.')
versions={'python':sys.version,'packages':{}}
for name in ['numpy','pandas','scipy','statsmodels','pandas-gbq','google-cloud-bigquery','IPython','threadpoolctl']:versions['packages'][name]=im.version(name)
(a.out/'runtime_versions.json').write_text(json.dumps(versions,indent=2))
inputs=[]
manifest=json.loads((ROOT/'manifest.json').read_text())
for rel,expected in manifest['files'].items():
 if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=expected:raise RuntimeError('Source checksum mismatch: '+rel)
def sha(src):
 h=hashlib.sha256()
 with open(src,'rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def link(src,name):
 src=Path(src).resolve();assert src.is_file(),str(src)
 dst=a.out/name;dst.symlink_to(src);inputs.append({'input':str(src),'size':src.stat().st_size,'sha256':sha(src)});return dst
def cell(i,ns):
 source=ROOT/'cells'/f'cell_{i:03}.py';exec(compile(source.read_text(),str(source),'exec'),ns)
ns={'__name__':'__main__','pd':pd,'np':np,'os':os,'display':lambda *args,**kwargs:None}
with open(a.out/'execution.private.log','w') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
 os.chdir(a.out)
 if a.stage=='check':
  for f in (ROOT/'cells').glob('*.py'):compile(f.read_text(),str(f),'exec')
 elif a.stage=='cohort':
  assert os.environ['WORKSPACE_CDR']=='wb-silky-artichoke-2408.C2024Q3R8','CDR mismatch'
  for i in [3,4,6]:cell(i,ns)
  if a.ancestry_file is None:raise ValueError('Provide --ancestry-file with the versioned provider ancestry TSV; see ancestry_source.json for the verified provider URI and generation.')
  link(a.ancestry_file,'ancestry_preds.tsv')
  ns['ancestry_df']=pd.read_csv('ancestry_preds.tsv',sep=chr(9))
  for i in [11,12,13,14,15,16]:cell(i,ns)
  for name,key in [('psoriasis_pheno.tsv','phenotype_out'),('psoriasis_covariates.tsv','covar_out')]:
   assert not ns[key][['FID','IID']].duplicated().any();ns[key].to_csv(name,sep=chr(9),index=False)
 elif a.stage in ['labs','lab-processing']:
  cd=(a.cohort_dir or a.input_root).resolve()
  for name in ['psoriasis_pheno.tsv','psoriasis_covariates.tsv']:link(cd/name,name)
  ph=pd.read_csv('psoriasis_pheno.tsv',sep=chr(9),dtype={'FID':str,'IID':str});cv=pd.read_csv('psoriasis_covariates.tsv',sep=chr(9),dtype={'FID':str,'IID':str})
  ns['analysis_base']=ph.merge(cv,on=['FID','IID'],validate='one_to_one');ns['analysis_base']['person_id']=ns['analysis_base'].IID.astype(str)
  if a.stage=='labs':
   assert os.environ['WORKSPACE_CDR']=='wb-silky-artichoke-2408.C2024Q3R8'
   for i in [95,97]:cell(i,ns)
  else:link(a.classified_labs or a.input_root/'psoriasis_biomarker_labs_classified.tsv','psoriasis_biomarker_labs_classified.tsv')
  for i in [98,99,100]:cell(i,ns)
 elif a.stage in ['x','biomarkers','ancestry']:
  if a.model_input_dir is None:raise ValueError('Provide --model-input-dir from build_followup_inputs.py; no implicit saved model matrix is used.')
  rr=a.model_input_dir.resolve()
  for name in ['X_inputs.private.pkl','biomarker_inputs.private.pkl','variant_map.tsv']:link(rr/name,name)
  if a.stage=='ancestry':link(rr/'ancestry_autosomal_inputs.private.pkl','ancestry_autosomal_inputs.private.pkl')
  source=ROOT/'retained/audit/relatedness_sensitivity_20260910'/({'x':'run_X_relatedness.py','biomarkers':'run_biomarkers.py','ancestry':'run_autosomal_ancestry.py'}[a.stage])
  target=a.out/source.name;shutil.copy2(source,target);runpy.run_path(str(target),run_name='__main__')
 elif a.stage=='sii':
  if a.model_input_dir is None:raise ValueError('Provide --model-input-dir from build_followup_inputs.py; no implicit saved model matrix is used.')
  rr=a.model_input_dir.resolve()
  for name in ['biomarker_inputs.private.pkl','variant_map.tsv']:link(rr/name,name)
  shutil.copy2(ROOT/'retained/audit/relatedness_sensitivity_20260910/run_biomarkers.py',a.out/'run_biomarkers.py')
  labs=link(a.classified_labs or a.input_root/'psoriasis_biomarker_labs_classified.tsv','classified.private.tsv')
  orig=link(a.input_root/'psoriasis_biomarker_person_level.tsv','original_person.private.tsv')
  ns.update({'Path':Path,'json':json,'A':a.out,'R':a.out,'threading':__import__('threading'),'labs':pd.read_csv(labs,sep=chr(9),dtype={'person_id':str},low_memory=False),'orig':pd.read_csv(orig,sep=chr(9),dtype={'person_id':str},low_memory=False),'bio':pd.read_pickle(a.out/'biomarker_inputs.private.pkl')})
  source=ROOT/'retained/audit/SII_timing_sensitivity_20260910/run_sii.py'
  exec(compile(source.read_text(),str(source),'exec'),ns);ns['worker'].join()
  if ns['state']['stage']!='complete':raise RuntimeError(ns['state'])
(a.out/'inputs.private.json').write_text(json.dumps(inputs,indent=2))
(a.out/'COMPLETE.json').write_text(json.dumps({'stage':a.stage,'completed':True,'source_manifest_sha256':sha(ROOT/'manifest.json')},indent=2))
print('Completed',a.stage,'in',a.out)
