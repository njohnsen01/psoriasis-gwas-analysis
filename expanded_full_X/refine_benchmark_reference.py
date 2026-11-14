from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,datetime
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x_rescue';N=O/'R_ml_benchmark_exact';N.mkdir(exist_ok=True)
for n in ['selected.private.csv','selection.private.csv']:
 p=N/n
 if not p.exists():p.symlink_to(O/'R_ml_benchmark'/n)
with (N/'console.log').open('w') as f:p=subprocess.run(['/opt/conda/envs/jupyter/bin/Rscript',str(R/'validate_stable_benchmark_exact_information.R'),str(N)],stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0;a=pd.read_csv(O/'benchmark.private.tsv',sep='\t');a=a[a.status=='estimated'];x=a.merge(pd.read_csv(N/'R_summary.csv'),on='ID',suffixes=('_python','_R'));assert len(x)==116 and (x.status_R=='estimated').all();delta=float(np.max(np.abs(x.WALD_CHISQ_3DF-x.wald)));assert delta<1e-4,delta
old=pd.read_csv(O/'R_ml_benchmark/R_summary.csv');y=x.merge(old[['ID','wald']],on='ID',suffixes=('_exact','_old'));change=float(np.max(np.abs(y.wald_exact-y.wald_old)))
r={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'passed','positive_controls':116,'separated_negative_controls':12,'max_R_Wald_difference':delta,'max_R_score':float(x.score.max()),'max_covariance_reference_refinement_change':change,'preservation':'InitialPLINKandRQRreferencegatefailuresretained;exactfinalFisherreferenceusedwithoutrelaxinganytolerance.'};(O/'R_ML_BENCHMARK_VALIDATED.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
