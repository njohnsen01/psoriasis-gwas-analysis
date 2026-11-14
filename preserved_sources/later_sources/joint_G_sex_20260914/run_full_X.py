from pathlib import Path
import json,subprocess,concurrent.futures,time,hashlib
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911')
O=B/'reconciliation_20260912/joint_G_sex_20260914/full_X';O.mkdir(exist_ok=True)
def run(p):
 d=json.loads(p.read_text())
 if not d.get('QC_passed'):return {'chunk':p.parent.name,'status':'prior_QC_not_passed'}
 out=O/p.parent.name;out.mkdir(exist_ok=True)
 if (out/'COMPLETE.json').exists():return json.loads((out/'COMPLETE.json').read_text())
 c=json.loads((p.parent/'male02.console.command.json').read_text());c[c.index('--out')+1]=str(out/'joint');c[c.index('--threads')+1]='2';c[c.index('--pThresh')+1]='0.999999';c+=['--interaction','sex','--catCovarList','sex']
 (out/'command.json').write_text(json.dumps(c,indent=2));t=time.time()
 with (out/'console.log').open('w') as f:r=subprocess.run(c,stdout=f,stderr=subprocess.STDOUT)
 z={'chunk':p.parent.name,'returncode':r.returncode,'seconds':time.time()-t}
 (out/('COMPLETE.json' if r.returncode==0 else 'FAILED.json')).write_text(json.dumps(z));return z
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
 for z in pool.map(run,sorted((B/'results/X_source').glob('*/COMPLETE.json'))):print(json.dumps(z),flush=True)
(O/'DRIVER_FINISHED.json').write_text(json.dumps({'finished':time.time()}))
