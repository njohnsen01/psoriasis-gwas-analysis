from pathlib import Path
import subprocess,time,os,datetime,json,hashlib
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013')
p=R/'FULL_X_FOLLOWUP_PROTOCOL.json'
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selection':'Top ten valid joint-Wald variants with successive leads more than 1 Mb apart, plus three prior NLGN4X targets; fixed before full scan completion','analyses':'Independent full-model maximum-likelihood Wald/LRT and constrained full-penalty Firth LRT for ADD, female heterozygosity, and ADDxsex; same expanded nuisance span','scope':'Follow-up checks of selected scan signals, not independent replication. Report all selected results and failures. Retain genome-wide 5e-8 threshold for candidate discovery; no selected-set multiple-testing shortcut.','limitations':'One-megabase spacing does not establish LD independence; implementation checks do not prove whole-scan tail calibration','sha256':{f:hashlib.sha256((R/f).read_bytes()).hexdigest() for f in ['followup_full_x.py','followup_full_x.R']}}
if not p.exists():p.write_text(json.dumps(protocol,indent=2)+'\n')
else:assert json.loads(p.read_text())['sha256']==protocol['sha256']
print('WAITING_FOR_FULL_X',flush=True)
while not (R/'full_x/SCAN_COMPLETE.json').exists():
 pid=int((R/'run_full_x.production.pid').read_text());stat=Path('/proc')/str(pid)/'stat'
 if not stat.exists() or stat.read_text().split()[2]=='Z':raise RuntimeError('Full scan ended without completion marker')
 time.sleep(15)
assert (R/'FIRTH_LIBRARY_READY').exists(),'Firth dependency restoration failed; inspect install log'
with (R/'followup_full_x.console.log').open('w') as f:p=subprocess.run(['/opt/conda/envs/jupyter/bin/python3',str(R/'followup_full_x.py')],stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
print('FOLLOWUP_RETURN_CODE',p.returncode,flush=True);assert p.returncode==0
