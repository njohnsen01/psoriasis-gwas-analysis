"""Numerical refinement of failed Firth lead checks; scientific criteria unchanged."""
from pathlib import Path
import os,json,subprocess,datetime,hashlib
import pandas as pd
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013')
F=R/'full_x_followup';I=R/'full_x_followup_initial';A=R/'aggregate_review'
def alive(pid):
 p=Path('/proc')/str(pid)/'stat'
 return p.exists() and p.read_text().split()[2]!='Z'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert all(not alive(int((R/(n+'.pid')).read_text())) for n in ['chain_full_x_followup','chain_full_x_review'])
 assert (F/'FOLLOWUP_COMPLETE.json').exists() and (A/'FULL_X_REVIEW.json').exists()
 assert not I.exists() and not (R/'aggregate_review_initial').exists()
 old=pd.read_csv(F/'followup_summary.csv');assert len(old)==13
 failures=old.loc[old.status!='estimated','status'];assert set(failures)<=set(['fs <= 1e-06 is not TRUE','ns <= 1e-06 is not TRUE'])
 assert len(failures)>0
 note={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Eight selected Firth fits failed independently computed full/null modified-score convergence gates','change':'Tighten native IRLS tolerances to lconv1e-12,gconv1e-10,xconv1e-12; maxit500,maxhs100; preserve all prior outputs; reuse five accepted checks','unchanged':'Same selected variants, cohort, nuisance design, likelihood/penalty, tests, thresholds and independent score<=1e-6 acceptance criterion','initial_R_sha256':sha(R/'followup_full_x.R'),'refinement_R_sha256':sha(R/'followup_full_x_refine.R')}
 (R/'FIRTH_NUMERICAL_REFINEMENT.json').write_text(json.dumps(note,indent=2)+'\n')
 F.rename(I);A.rename(R/'aggregate_review_initial');F.mkdir()
 (F/'selected.private.csv').symlink_to(I/'selected.private.csv')
 (F/'selected_variants.json').symlink_to(I/'selected_variants.json')
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
 with (F/'console.log').open('w') as log:
  p=subprocess.run(['/opt/conda/envs/jupyter/bin/Rscript',str(R/'followup_full_x_refine.R'),str(F),str(R.parent/'r_library'),'13',str(I)],stdout=log,stderr=subprocess.STDOUT,env=env)
 assert p.returncode==0,'Refinement R process failed; preserve all outputs and inspect'
 audit=pd.DataFrame(json.loads((I/'selected_variants.json').read_text()))
 result=pd.read_csv(F/'followup_summary.csv').merge(audit,on='key',validate='many_to_one');result.to_csv(F/'followup_named.csv',index=False)
 receipt={'status':'independent_lead_checks_completed_pending_review','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selected':13,'failed':int((result.status!='estimated').sum()),'refinement':'Stricter native convergence controls; independent acceptance criteria unchanged','limitation':'Selected follow-up is not independent replication'}
 (F/'FOLLOWUP_COMPLETE.json').write_text(json.dumps(receipt,indent=2)+'\n')
 with (R/'review_full_x.refined.console.log').open('w') as log:
  p=subprocess.run(['/opt/conda/envs/jupyter/bin/python3',str(R/'review_full_x.py')],stdout=log,stderr=subprocess.STDOUT,env=env)
 assert p.returncode==0,'Refined aggregate review failed'
 (R/'FIRTH_REFINEMENT_COMPLETE.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
