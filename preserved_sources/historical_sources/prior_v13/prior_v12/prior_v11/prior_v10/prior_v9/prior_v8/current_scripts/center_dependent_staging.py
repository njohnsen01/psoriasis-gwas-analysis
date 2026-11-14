from pathlib import Path
import json,hashlib,datetime,time,subprocess,sys,os,traceback
import pandas as pd,numpy as np,psutil
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';O=I/'center_dependent_staging';O.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,value):
 p=O/name;t=Path(str(p)+'.partial');t.write_text(json.dumps(value,indent=2,default=str));t.replace(p)
def status(stage,**kw):save('STATUS.private.json',dict(stage=stage,updated_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),pid=os.getpid(),release='HELD',**kw))
def run(script,label):
 with (O/(label+'.private.log')).open('w') as f:r=subprocess.run([sys.executable,str(script)],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
 if r.returncode:raise RuntimeError(label+' failed; inspect secure log')
try:
 contract=json.loads((O/'CONTRACT.private.json').read_text());center=I/'center_refit_v2';done=center/'COMPLETE.private.json';deadline=time.monotonic()+24*3600
 status('waiting_for_existing_center_refit',center_pid=contract['center_pid'])
 while not done.exists():
  if (center/'FAILED.private.json').exists():raise RuntimeError('Center refit failed; dependent staging not run')
  try:
   p=psutil.Process(contract['center_pid']);assert abs(p.create_time()-contract['center_create_time'])<.01;assert str(I/'scripts/center_refit_v2.py') in p.cmdline()
   if p.status()==psutil.STATUS_ZOMBIE:raise RuntimeError('Center process exited without validated completion')
  except psutil.NoSuchProcess:
   if done.exists():break
   raise RuntimeError('Center process absent without validated completion')
  if time.monotonic()>deadline:raise TimeoutError('Center-dependent staging deadline reached')
  time.sleep(20)
 z=json.loads(done.read_text());assert z['status']=='verified' and z['exact_primary_cohort'] and z['center_refit_both_stages'] and z['matched_primary_NLGN4X_reproduced'];assert z['loco_rows']==list(range(1,24));assert z['script_sha256']==contract['center_script_sha256']==sha(I/'scripts/center_refit_v2.py')
 d=pd.read_csv(center/'estimates.private.tsv',sep='\t');assert len(d)==16 and not d.duplicated(['adjustment','model','ID']).any();assert np.isfinite(d[['BETA','SE','LOG10P','OR','CI95_low','CI95_high']]).all().all();assert d.SE.gt(0).all()
 for path,expected in contract['builder_hashes'].items():assert sha(Path(path))==expected,'Builder changed after dependency launch'
 status('building_center_integrated_secure_supplements',center_receipt_sha256=sha(done))
 # Build a new candidate; preserve v2 and frozen primary sources.
 src=(I/'scripts/build_S1_secure_staging_v2.py').read_text().replace("I/'supplement_staging_v2'","I/'supplement_staging_v3'")
 src=src.replace('wb=openpyxl.Workbook();',"specs += [('Center_refit',I/'center_refit_v2/estimates.private.tsv','Both-stage center refit and matched baseline; QC-failed leads diagnostic only'),('Center_comparison',I/'center_refit_v2/paired_comparison.private.tsv','Matched participant contrast; not independent replication or genome-wide recalibration')]\nwb=openpyxl.Workbook();")
 src=src.replace('Center refit; final retained-locus','Final retained-locus')
 p=I/'scripts/build_S1_secure_staging_v3.py';p.write_text(src);compile(src,str(p),'exec');run(p,'S1_v3')
 src=(I/'scripts/build_S2_secure_integration_v1.py').read_text().replace("I/'supplement_staging_v2'","I/'supplement_staging_v3'")
 src=src.replace('Its completion and matched baseline concordance must be verified before numerical integration.','Its completion receipt verifies both stages, the same eligible cohort, all chromosome prediction rows and matched primary NLGN4X baseline concordance. Targeted center adjustment does not demonstrate genome-wide recalibration.')
 p=I/'scripts/build_S2_secure_integration_v2.py';p.write_text(src);compile(src,str(p),'exec');run(p,'S2_v2')
 target=I/'supplement_staging_v3/S2_reconciliation_integration.private.md'
 keep=['adjustment','model','ID','ALLELE1','primary_QC_PASS','BETA','SE','OR','CI95_low','CI95_high','LOG10P','P_from_REGENIE_LOG10P','N'];t=d[keep]
 def cell(v):return str(v).replace('|','/').replace(chr(10),' ')
 table='| '+' | '.join(keep)+' |\n| '+' | '.join(['---']*len(keep))+' |\n'+'\n'.join('| '+' | '.join(cell(v) for v in row)+' |' for row in t.itertuples(index=False,name=None))
 target.write_text(target.read_text()+'\n## Verified center-refit estimates\n\nPrivate, source-verified results. Male0/2 and male0/1 denote coded dosage units. Rows with primary_QC_PASS=False are historical diagnostic targets, not restored primary discoveries. Exact variant N remains subject to disclosure review.\n\n'+table+'\n')
 src=(I/'scripts/render_S1_secure_previews.py').read_text().replace('reconciliation_20260912/supplement_staging_v2','reconciliation_20260912/supplement_staging_v3')
 p=I/'scripts/render_S1_secure_previews_v3.py';p.write_text(src);compile(src,str(p),'exec');run(p,'previews_v3')
 save('COMPLETE.private.json',dict(status='center_integrated_private_staging_built_and_numerically_verified',center_receipt_sha256=sha(done),S1=str(I/'supplement_staging_v3/S1_reconciled_tables_HELD.xlsx'),S1_sha256=sha(I/'supplement_staging_v3/S1_reconciled_tables_HELD.xlsx'),S2=str(target),S2_sha256=sha(target),visual_review='new center rows and regenerated previews require human/model review',native_Excel_Word_QA='pending',final_submission_ready=False,release='HELD'))
 status('completed_secure_staging_review_pending')
except Exception as e:
 save('FAILED.private.json',dict(error=repr(e),traceback=traceback.format_exc(),release='HELD'));status('failed');raise
