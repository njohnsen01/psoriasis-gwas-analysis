from pathlib import Path
import time,subprocess,os
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013')
print('WAITING_FOR_FULL_X_FOLLOWUP',flush=True)
while not (R/'full_x_followup/FOLLOWUP_COMPLETE.json').exists():
 pid=int((R/'chain_full_x_followup.pid').read_text());stat=Path('/proc')/str(pid)/'stat'
 if not stat.exists() or stat.read_text().split()[2]=='Z':raise RuntimeError('Lead follow-up ended without completion marker; inspect logs')
 time.sleep(15)
with (R/'review_full_x.console.log').open('w') as f:p=subprocess.run(['/opt/conda/envs/jupyter/bin/python3',str(R/'review_full_x.py')],stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
print('REVIEW_RETURN_CODE',p.returncode,flush=True);assert p.returncode==0
