from pathlib import Path
import os,json,re,time,hashlib,subprocess,shutil,traceback,threading,sys
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np,pandas as pd,zstandard as zstd
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parent
OLD=Path(f'/home/jupyter/psoriasis_regenie_loco_step2/attempt_{os.environ["PSO_PRIMARY_RUN_ID"]}')
PLAN=json.loads((OLD/'status/step2_plan.json').read_text())
STRATA=json.loads((ROOT/'strata.json').read_text())
PANEL=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/stratified_calibration/private_panel/calibration_pruned')
LOCK=threading.Lock()
STATUS={'status':'waiting_for_panel','started_utc':datetime.now(timezone.utc).isoformat(),'step1':{},'chunks_completed':0,'expected_chunks':1149,'strata':list(STRATA),'failed_chunks':{}}
THREAD_ENV=dict(os.environ,OMP_NUM_THREADS='8',MKL_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1')
for p in ['step1','step2','scratch','receipts','diagnostics']:(ROOT/p).mkdir(exist_ok=True)
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def atomic(p,d):
 q=p.with_suffix(p.suffix+'.partial');q.write_text(json.dumps(d,indent=2));q.replace(p)
def status(**kw):
 with LOCK:
  STATUS.update(kw);STATUS['updated_utc']=now();atomic(ROOT/'pipeline_status.json',STATUS)
def run(cmd,log,threads=8):
 env=dict(THREAD_ENV,OMP_NUM_THREADS=str(threads),MKL_NUM_THREADS=str(threads))
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
 if r.returncode:raise RuntimeError('Command failed: '+str(cmd[0])+' '+log.read_text()[-5000:])
def validate_samples(name,log):
 d=STRATA[name]
 for text in ['number of individuals used in analysis = '+str(d['n']),str(d['cases'])+' cases and '+str(d['controls'])+' controls']:
  if text not in log:raise RuntimeError(name+' sample-count validation failed: '+text)
def step1(name):
 rec=ROOT/'step1'/(name+'.validated.json')
 if rec.exists():return json.loads(rec.read_text())
 folder=ROOT/'step1'/name;folder.mkdir(exist_ok=True);out=folder/'model'
 inp=ROOT/'private_inputs';ph=inp/(name+'.pheno.tsv');co=inp/(name+'.covar.tsv');keep=inp/(name+'.keep.tsv')
 cov=pd.read_csv(co,sep='\t');assert cov.columns.tolist()==['FID','IID']+STRATA[name]['covariates']
 assert len(cov)==STRATA[name]['n'] and np.isfinite(cov.iloc[:,2:].to_numpy()).all()
 c=cov.iloc[:,2:].to_numpy();c=(c-c.mean(0))/c.std(0);assert np.isfinite(c).all() and np.linalg.matrix_rank(np.c_[np.ones(len(c)),c])==c.shape[1]+1
 cmd=['regenie','--step','1','--bt','--bed',str(PANEL),'--keep',str(keep),'--phenoFile',str(ph),'--phenoCol','psoriasis_case','--covarFile',str(co),'--bsize','2000','--lowmem','--lowmem-prefix',str(folder/'lowmem'),'--threads',str(STRATA[name]['step1_threads']),'--out',str(out)]
 atomic(folder/'command.json',{'command':cmd,'phenotype_sha256':sha(ph),'covariates_sha256':sha(co),'keep_sha256':sha(keep)})
 print('STEP1_START',name,flush=True);run(cmd,folder/'console.log',threads=STRATA[name]['step1_threads']);log=Path(str(out)+'.log').read_text();validate_samples(name,log)
 lf=Path(str(out)+'_1.loco');expected=pd.read_csv(keep,sep='\t',header=None,dtype=str);ids=set(expected[0]+'_'+expected[1])
 with lf.open() as f:
  header=f.readline().split();assert header[0]=='FID_IID' and len(header)-1==len(ids) and set(header[1:])==ids
  rows=[]
  for line in f:
   tokens=line.split();rows.append(int(tokens[0]));v=np.array(tokens[1:],dtype=float);assert len(v)==len(ids) and np.isfinite(v).all()
 assert rows==list(range(1,24))
 d={'status':'validated','completed_utc':now(),'stratum':name,'n':len(ids),'cases':STRATA[name]['cases'],'controls':STRATA[name]['controls'],'loco_rows':rows,'loco_sha256':sha(lf),'prediction_list':str(out)+'_pred.list','log_sha256':sha(Path(str(out)+'.log'))}
 atomic(rec,d)
 with LOCK:STATUS['step1'][name]='validated'
 status();print('STEP1_VALIDATED',name,flush=True);return d

try:
 status(status='refitting_after_provider_coverage_exclusion',expected_chunks=0)
 with ThreadPoolExecutor(max_workers=6) as ex:
  for f in as_completed([ex.submit(step1,n) for n in STRATA]): f.result()
 status(status='step1_complete',completed_utc=now())
except Exception as e:
 status(status='failed',error=str(e),traceback=traceback.format_exc());raise
