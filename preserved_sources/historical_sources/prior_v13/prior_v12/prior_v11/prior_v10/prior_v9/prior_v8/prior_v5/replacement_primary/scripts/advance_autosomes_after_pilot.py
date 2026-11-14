from pathlib import Path
import json,os,subprocess,time,traceback,psutil
N=Path(__file__).resolve().parents[1]
try:
 pid=json.loads((N/'status/autosomal_pid.json').read_text())['pid']
 while psutil.pid_exists(pid):
  try:
   if psutil.Process(pid).status()==psutil.STATUS_ZOMBIE:break
  except psutil.NoSuchProcess:break
  time.sleep(30)
 status=json.loads((N/'status/autosomal_status.json').read_text())
 assert status['status']=='complete' and status['completed']==status['expected']==2,status
 records=[json.loads(p.read_text()) for p in (N/'results/autosomes').glob('*/COMPLETE.json')]
 assert len(records)==2 and all(r['status']=='validated' and r.get('failed_or_unavailable',0)==0 for r in records)
 (N/'status/autosomal_pilot_validated.json').write_text(json.dumps({'status':'validated','chunks':2,'genotype_counts_independently_checked':True,'sample_counts_alleles_and_result_rows_checked':True,'failed_tests':0},indent=2))
 with (N/'logs/autosomal_full.log').open('w') as log:r=subprocess.run(['python',str(N/'scripts/run_autosomes.py')],stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1'))
 if r.returncode:raise RuntimeError('Full autosomal run failed; inspect its status and log')
 print('AUTOSOMAL_FULL_RUN_COMPLETE',flush=True)
except Exception as e:
 (N/'status/autosomal_advance_failed.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc()}));raise
