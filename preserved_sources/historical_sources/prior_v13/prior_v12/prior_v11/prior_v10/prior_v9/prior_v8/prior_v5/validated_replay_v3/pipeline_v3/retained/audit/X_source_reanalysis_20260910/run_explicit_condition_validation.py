from pathlib import Path
import json,subprocess,os,hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent
jobs=json.loads((R/'private_explicit_condition_inputs/jobs.json').read_text())
D=R/'explicit_condition_validation'
def run(j):
 cmd=j['command']; out=Path(cmd[cmd.index('--out')+1]); (D/(j['label']+'.command.json')).write_text(json.dumps(cmd,indent=2))
 with open(str(out)+'.console.log','w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
 assert r.returncode==0,Path(str(out)+'.console.log').read_text()[-1800:]
 d=pd.read_csv(str(out)+'_psoriasis_case.regenie',sep=r'\s+'); q=d[d.ID==j['target']]; assert len(q)==1 and np.isfinite(q.LOG10P.iloc[0]); z=q.iloc[0]
 result={'analysis':j['label'],'target':j['target'],'BETA':float(z.BETA),'SE':float(z.SE),'P':float(10**(-z.LOG10P)),'scope':'explicit model-matched conditioning dosage; complete Xq12 pair QC subset; corrected coverage-cohort LOCO retained'}
 (D/(j['label']+'.validated.json')).write_text(json.dumps(result,indent=2)); return result
with ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(run,jobs))
(D/'complete.json').write_text(json.dumps(res,indent=2))
print('EXPLICIT_CONDITIONAL_COMPLETE',len(res),flush=True)
