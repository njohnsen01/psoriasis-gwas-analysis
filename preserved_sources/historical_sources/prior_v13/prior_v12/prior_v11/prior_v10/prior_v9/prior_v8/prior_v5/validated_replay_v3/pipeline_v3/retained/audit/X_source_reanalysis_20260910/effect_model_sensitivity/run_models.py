from pathlib import Path
import subprocess,json,traceback,os
from concurrent.futures import ThreadPoolExecutor,as_completed
R=Path(__file__).resolve().parent
jobs=json.loads((R/'all_jobs.json').read_text())
def run(j):
 label=j['label']; cmd=j['cmd'];keep=Path(cmd[cmd.index('--keep')+1]);assert keep.exists() and keep.stat().st_size>0
 with open(R/(label+'.console.log'),'w') as f: r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1'))
 if r.returncode:raise RuntimeError(label+': '+(R/(label+'.console.log')).read_text()[-2500:])
 (R/(label+'.complete.json')).write_text(json.dumps({'status':'complete','command':cmd}))
 print('COMPLETE',label,flush=True)
try:
 with ThreadPoolExecutor(max_workers=4) as ex:
  for f in as_completed([ex.submit(run,j) for j in jobs]):f.result()
 (R/'complete.json').write_text(json.dumps({'status':'complete','jobs':len(jobs)}))
except Exception as e:
 (R/'error.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc()}));raise
