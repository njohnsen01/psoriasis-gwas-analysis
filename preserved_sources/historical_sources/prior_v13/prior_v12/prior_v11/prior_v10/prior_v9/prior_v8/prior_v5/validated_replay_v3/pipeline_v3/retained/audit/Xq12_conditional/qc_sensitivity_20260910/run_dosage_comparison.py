import os
from pathlib import Path
from datetime import datetime, timezone
import subprocess, json, hashlib, pandas as pd
base=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/Xq12_conditional')
work=base/'qc_sensitivity_20260910'
private=base/'private_genotypes'
plan=json.loads(Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}/status/step2_plan.json').read_text())
reg='/opt/workbench-tools/binaries/bin/regenie'
assert subprocess.run([reg,'--version'],capture_output=True,text=True).stdout.strip()==plan['regenie_version']
original=[private/'chrX_pso_chunk_0022.pgen',private/'chrX_pso_chunk_0022.pvar.zst',private/'chrX_pso_step2.psam']
assert all(p.exists() for p in original)
status={'started_utc':datetime.now(timezone.utc).isoformat(),'state':'running','primary_inputs_unchanged':True,'arms':{}}
(work/'status.json').write_text(json.dumps(status,indent=2))
extract=['/opt/workbench-tools/binaries/bin/plink2','--pgen',str(original[0]),'--pvar',str(original[1]),'--psam',str(original[2]),'--extract',str(base/'two_ids.txt'),'--make-pgen','--threads','2','--memory','4000','--out',str(private/'retained_same_two')]
r=subprocess.run(extract,capture_output=True,text=True)
(work/'retained_extraction.log').write_text(r.stdout+'\n'+r.stderr)
if r.returncode: raise RuntimeError('Retained two-variant extraction failed; see private log')
ps1=pd.read_csv(private/'retained_same_two.psam',sep=r'\s+',dtype=str)
ps2=pd.read_csv(private/'source_qc'/'native_model.psam',sep=r'\s+',dtype=str)
assert ps1.IID.is_unique and ps2.IID.is_unique and set(ps1.IID)==set(ps2.IID)
sx=ps1[['IID','SEX']].merge(ps2[['IID','SEX']],on='IID',suffixes=('_a','_b'),validate='one_to_one')
assert (sx.SEX_a==sx.SEX_b).all()
status['sample_set_identical']=True; status['sex_assignments_identical']=True
status['input_hashes']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in original+[Path(plan['phenotype_file']),Path(plan['covariate_file']),Path(plan['prediction_list'])]}
commands=[]
for name,prefix in [('retained',private/'retained_same_two'),('native_ploidy',private/'source_qc'/'native_model')]:
 cmd=[reg,'--step','2','--pgen',str(prefix),'--phenoFile',plan['phenotype_file'],'--phenoCol',plan['phenotype'],'--covarFile',plan['covariate_file'],'--covarColList',','.join(plan['covariates']),'--bt','--firth','--approx','--pThresh','0.01','--af-cc','--pred',plan['prediction_list'],'--par-region','hg38','--bsize','400','--threads','2','--out',str(work/name)]
 commands.append(cmd)
 (work/'commands.json').write_text(json.dumps(commands,indent=2))
 with (work/f'{name}.console.log').open('w') as log:
  ret=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 status['arms'][name]={'returncode':ret.returncode,'completed_utc':datetime.now(timezone.utc).isoformat()}
 (work/'status.json').write_text(json.dumps(status,indent=2))
 if ret.returncode: raise RuntimeError('REGENIE diagnostic failed; see private log')
status['state']='completed'
status['finished_utc']=datetime.now(timezone.utc).isoformat()
(work/'status.json').write_text(json.dumps(status,indent=2))
