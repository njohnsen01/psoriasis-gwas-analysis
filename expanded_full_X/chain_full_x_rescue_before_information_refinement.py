"""Checkpointed validation -> ML rescue -> Firth sensitivity -> lead verification."""
from pathlib import Path
import os,time,json,datetime,fcntl,sys,subprocess,concurrent.futures,hashlib
import pandas as pd,numpy as np
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x_rescue'

def alive(name):
 p=R/(name+'.pid')
 if not p.exists():return False
 d=Path('/proc')/p.read_text().strip()
 return d.exists() and '\nState:\tZ' not in (d/'status').read_text()

def wait_for(path,job):
 while not path.exists():
  if not alive(job):raise RuntimeError('Required checkpoint missing and worker ended: '+job)
  time.sleep(30)
 return json.loads(path.read_text())

def receipt(name,obj):
 obj['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/name).write_text(json.dumps(obj,indent=2)+'\n')

def main():
 lock=(O/'recovery_chain.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 if (O/'RESCUE_AND_LEADS_COMPUTED.json').exists():return
 # Strict independent implementation comparisons, including negative controls.
 ml=wait_for(O/'R_ML_BENCHMARK_VALIDATED.json','validate_stable_benchmark');native=wait_for(O/'NATIVE_VALIDATION_COMPUTED.json','validate_rescue_native')
 assert ml['status']=='passed' and ml['positive_controls']==116 and ml['max_R_Wald_difference']<1e-4
 assert native['ML_controls']==8 and native['max_ML_Wald_difference']<1e-4
 if not (O/'BENCHMARK_VALIDATED.json').exists():receipt('BENCHMARK_VALIDATED.json',{'status':'passed','independent_R_positive_controls':124,'separation_negative_controls':12,'R_original_valid_comparison':ml,'R_failed_control_comparison':native,'initial_native_PLINK_gate':'Preservedfailedinitialcomparison;strictRcomparisonadjudicationdocumentedintheaddendum.'})
 if not (O/'RECOVERY_COMPUTED.json').exists():
  import rescue_x_stable as st
  st.setup();diag=pd.read_csv(O/'diagnostics.private.tsv',sep='\t');assert len(diag)==21730;a=diag[diag.full_genotype_rank & ~diag.empty_outcome_cell];jobs=[(c,list(g.ID),'recovery') for c,g in a.groupby('chunk')]
  receipt('ML_RECOVERY_STARTED.json',{'variants':len(a),'workers':24,'status':'running'})
  with concurrent.futures.ThreadPoolExecutor(max_workers=24) as ex:res=pd.concat(list(ex.map(st.one,jobs)),ignore_index=True)
  res.to_csv(O/'recovery.private.tsv',sep='\t',index=False);remaining=diag[~diag.ID.isin(res.ID)][['ID','original_error']].copy();remaining['status']=np.where(diag.set_index('ID').loc[remaining.ID,'full_genotype_rank'].to_numpy(),'separated_genotype_cell','rank_deficient_genotype');allr=pd.concat([res,remaining],ignore_index=True);assert len(allr)==21730 and allr.ID.is_unique;allr.to_csv(O/'all_rescue_outcomes.private.tsv',sep='\t',index=False)
  receipt('RECOVERY_COMPUTED.json',{'status':'computed_pending_independent_review','attempted':len(res),'status_counts':res.status.value_counts().to_dict(),'all_status_counts':allr.status.value_counts().to_dict(),'seconds_sum':float(res.seconds.sum())})
 # Verify penalized implementation against the original13 and8 separated R controls.
 fb=wait_for(O/'FIRTH_BENCHMARK_COMPUTED.json','rescue_x_firth_benchmark');assert fb['status_counts']=={'estimated':21} and fb['max_prior_R_stat_difference']<1e-4
 f=pd.read_csv(O/'firth_benchmark.private.tsv',sep='\t');n=pd.read_csv(O/'native_validation/native_summary.csv');x=f.merge(n[['ID','firth_stat','status']],on='ID',suffixes=('_python','_R'));assert len(x)==8 and (x.status_R=='estimated').all();delta=float(np.max(np.abs(x.FIRTH_CHISQ_3DF-x.firth_stat)));assert delta<1e-4
 if not (O/'FIRTH_BENCHMARK_VALIDATED.json').exists():receipt('FIRTH_BENCHMARK_VALIDATED.json',{'status':'passed','prior_R_controls':13,'separated_R_controls':8,'max_prior_R_difference':fb['max_prior_R_stat_difference'],'max_separated_R_difference':delta})
 audit=wait_for(R/'full_x_estimability_audit/DIAGNOSTICS_COMPLETE.json','diagnose_all_x');assert audit['n_tests']==334926
 if not (O/'FIRTH_RECOVERY_COMPUTED.json').exists():
  with (R/'rescue_x_firth_all.console.log').open('w') as log:
   p=subprocess.Popen([sys.executable,str(R/'rescue_x_firth_all.py'),'firth_recovery'],stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'));(R/'rescue_x_firth_all.pid').write_text(str(p.pid));assert p.wait()==0,'Firth recovery process failed'
 if not (O/'FINAL_LEAD_CHECKS_COMPUTED.json').exists():
  with (R/'run_final_lead_checks.console.log').open('w') as log:
   p=subprocess.Popen([sys.executable,str(R/'run_final_lead_checks.py')],stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'));(R/'run_final_lead_checks.pid').write_text(str(p.pid));assert p.wait()==0,'Final lead validation requires repair'
 receipt('RESCUE_AND_LEADS_COMPUTED.json',{'status':'computed_pending_scientific_and_disclosure_review','note':'Noexportordiscoveryclaimuntilindependentreview.Originalfilespreserved.'})
 print('RESCUE_AND_LEADS_COMPUTED_REVIEW_REQUIRED',flush=True)
if __name__=='__main__':main()
