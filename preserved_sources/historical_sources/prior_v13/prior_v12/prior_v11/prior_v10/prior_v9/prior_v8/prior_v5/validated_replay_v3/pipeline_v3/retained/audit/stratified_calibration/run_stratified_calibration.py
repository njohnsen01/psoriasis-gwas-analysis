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
PANEL=ROOT/'private_panel/calibration_pruned'
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
def gs_copy(uri,path):
 run(['gsutil','-u',os.environ['GOOGLE_PROJECT'],'cp',uri,str(path)],Path(str(path)+'.copy.log'))
def tasks():
 out=[]
 for cr in PLAN['chunk_records']:
  pv={Path(u).name.removesuffix('.pvar.zst'):u for u in cr['pvar_uris']}
  for u in cr['pgen_uris']:
   stem=Path(u).name.removesuffix('.pgen');out.append({'chromosome':cr['chromosome'],'stem':stem,'pgen_uri':u,'pvar_uri':pv[stem],'psam_uri':cr['psam_uri']})
 assert len(out)==1149
 return sorted(out,key=lambda d:(0 if d['chromosome']=='X' else int(d['chromosome']),d['stem']))
TASKS=tasks()
def chunk(task):
 stem=task['stem'];ch=task['chromosome'];receipt=ROOT/'receipts'/(stem+'.json')
 if receipt.exists():return json.loads(receipt.read_text())
 scratch=ROOT/'scratch'/stem;scratch.mkdir(exist_ok=True);pref=scratch/stem
 try:
  if shutil.disk_usage(ROOT).free<150*1024**3:raise RuntimeError('Protected free disk threshold reached')
  gs_copy(task['pgen_uri'],Path(str(pref)+'.pgen'));gs_copy(task['pvar_uri'],Path(str(pref)+'.pvar.zst'));gs_copy(task['psam_uri'],Path(str(pref)+'.psam'))
  with Path(str(pref)+'.pvar.zst').open('rb') as f,Path(str(pref)+'.pvar').open('wb') as g:zstd.ZstdDecompressor().copy_stream(f,g)
  pv=pd.read_csv(Path(str(pref)+'.pvar'),sep='\t',usecols=['ID']);assert pv.ID.is_unique;source_ids=set(pv.ID)
  ps=pd.read_csv(Path(str(pref)+'.psam'),sep=r'\s+',dtype=str).rename(columns={'#FID':'FID','#IID':'IID'})
  assert len(ps)==392155 and ps.IID.is_unique
  covall=pd.read_csv(PLAN['covariate_file'],sep='\t',usecols=['IID','sex'],dtype={'IID':str})
  chk=ps[['IID','SEX']].merge(covall,on='IID',validate='one_to_one');assert len(chk)==392155 and (chk.SEX.astype(int)==np.where(chk.sex==1,1,2)).all()
  result={'chromosome':ch,'stem':stem,'n_source_variants':len(pv),'source_pgen_uri':task['pgen_uri'],'source_pgen_sha256':sha(Path(str(pref)+'.pgen')),'psam_sex_matches_frozen_covariate':True,'strata':{}}
  for name,d in STRATA.items():
   od=ROOT/'step2'/name/('chr'+ch);od.mkdir(parents=True,exist_ok=True);out=od/stem;inp=ROOT/'private_inputs'
   pred=json.loads((ROOT/'step1'/(name+'.validated.json')).read_text())['prediction_list']
   cmd=['regenie','--step','2','--bt','--pgen',str(pref),'--keep',str(inp/(name+'.keep.tsv')),'--phenoFile',str(inp/(name+'.pheno.tsv')),'--phenoCol','psoriasis_case','--covarFile',str(inp/(name+'.covar.tsv')),'--bsize','400','--firth','--approx','--pThresh','0.01','--minMAC','5','--pred',pred,'--af-cc','--threads','12','--gz','--out',str(out)]
   if ch=='X':cmd+=['--par-region','hg38']
   run(cmd,Path(str(out)+'.console.log'),threads=12);log=Path(str(out)+'.log').read_text();validate_samples(name,log)
   rf=Path(str(out)+'_psoriasis_case.regenie.gz');df=pd.read_csv(rf,sep=r'\s+',dtype={'ID':str})
   assert df.ID.is_unique and set(df.ID)<=source_ids and set(df.TEST)<=set(['ADD'])
   lp=pd.to_numeric(df.LOG10P,errors='coerce');valid=np.isfinite(lp)&(lp>=0)
   fail=re.search(r'Number of failed tests\s*:\s*\((\d+)/(\d+)\)',log)
   ignored=re.search(r'Number of ignored tests due to low MAC\s*:\s*(\d+)',log)
   result['strata'][name]={'n_rows':len(df),'n_valid':int(valid.sum()),'n_invalid':int((~valid).sum()),'n_source_without_row':len(pv)-len(df),'failed_firth_tests':int(fail.group(1)) if fail else 0,'firth_tests':int(fail.group(2)) if fail else 0,'ignored_low_MAC':int(ignored.group(1)) if ignored else 0,'result_sha256':sha(rf),'command':cmd,'result_file':str(rf)}
  result['completed_utc']=now();atomic(receipt,result);return result
 finally:
  shutil.rmtree(scratch,ignore_errors=True)
def summarize(chromosome_x_only=False):
 label='X' if chromosome_x_only else 'genome';diag=ROOT/'diagnostics';rows=[]
 for name in STRATA:
  fs=sorted((ROOT/'step2'/name).glob(('chrX' if chromosome_x_only else 'chr*')+'/*.regenie.gz'))
  expected=sum(1 for t in TASKS if not chromosome_x_only or t['chromosome']=='X')
  assert len(fs)==expected,(name,len(fs),expected)
  pieces=[]
  for f in fs:
   d=pd.read_csv(f,sep=r'\s+',usecols=['CHROM','GENPOS','ID','A1FREQ','N','LOG10P','CHISQ','EXTRA'],dtype={'ID':str,'EXTRA':str})
   pieces.append(d)
  d=pd.concat(pieces,ignore_index=True);assert d.ID.is_unique
  lp=pd.to_numeric(d.LOG10P,errors='coerce');valid=np.isfinite(lp)&(lp>=0);d=d.loc[valid].copy();d.LOG10P=lp[valid].to_numpy()
  ch=d.CHROM.astype(str).replace({'23':'X'});x=ch=='X';pos=d.GENPOS;mhc=(ch=='6')&pos.between(25000000,34000000);par=x&((pos<=2781479)|(pos>=155701383))
  groups={'X':x,'X_nonPAR':x&~par,'X_PAR':par}
  if not chromosome_x_only:groups={'all':np.ones(len(d),dtype=bool),'autosomes':~x,'autosomes_no_MHC':~x&~mhc,**groups}
  qq=[]
  for key,mask in groups.items():
   v=np.sort(d.loc[mask,'LOG10P'].to_numpy())[::-1];n=len(v)
   if not n:
    rows.append({'stratum':name,'group':key,'n_tests':0,'lambda_GC':None});continue
   lam=float(np.median(chi2.isf(np.power(10.,-v),1))/chi2.ppf(.5,1));direct=float(np.median(d.loc[mask,'CHISQ'])/chi2.ppf(.5,1))
   rows.append({'stratum':name,'group':key,'n_tests':n,'lambda_GC':lam,'lambda_from_CHISQ':direct,'n_p_lt_5e8':int((v>=-np.log10(5e-8)).sum()),'n_P_lt_0_05':int((v>-np.log10(.05)).sum())})
   ranks=np.unique(np.r_[np.arange(1,min(n,5000)+1),np.rint(np.geomspace(1,n,10000)).astype(int)])
   qq.append(pd.DataFrame({'stratum':name,'group':key,'rank':ranks,'n_tests':n,'expected':-np.log10((ranks-.5)/n),'observed':v[ranks-1]}))
  pd.concat(qq,ignore_index=True).to_csv(diag/(name+'_'+label+'_exact_QQ.tsv'),sep='\t',index=False)
  candidates=d[d.ID.isin(['chrX:22917619:C:T','chrX:5998724:A:G','chrX:67842269:G:T','chrX:67736589:T:C'])]
  candidates.to_csv(diag/(name+'_'+label+'_X_leads.tsv'),sep='\t',index=False)
  del d,pieces
 pd.DataFrame(rows).to_csv(diag/(label+'_inflation.tsv'),sep='\t',index=False)
 atomic(diag/(label+'_calibration.json'),{'completed_utc':now(),'scope':label,'full_output':True,'sampling_weights':False,'plot_ranks_computed_before_display_thinning':True,'groups':rows})
 print('DIAGNOSTICS_COMPLETE',label,flush=True)
def process_batch(ts,label):
 status(status='running_'+label)
 errors={}
 with ThreadPoolExecutor(max_workers=3) as ex:
  futures={ex.submit(chunk,t):t for t in ts}
  for f in as_completed(futures):
   t=futures[f]
   try:
    r=f.result();STATUS['chunks_completed']+=1;print('CHUNK_COMPLETE',t['stem'],STATUS['chunks_completed'],flush=True)
   except Exception as e:
    errors[t['stem']]=str(e);print('CHUNK_FAILED',t['stem'],str(e),flush=True)
   status(failed_chunks=errors)
 if errors:raise RuntimeError(label+' chunk failures: '+str(errors))
try:
 status()
 while True:
  p=ROOT/'panel_status.json'
  if p.exists():
   s=json.loads(p.read_text())
   if s['status']=='validated':break
   if s['status']=='failed':raise RuntimeError('Panel reconstruction failed: '+s['error'])
  time.sleep(20)
 versions={b:subprocess.run([b,'--version'],capture_output=True,text=True).stdout.strip() for b in ['regenie','plink2']}
 atomic(ROOT/'analysis_protocol.json',{'created_utc':now(),'versions':versions,'strata':STRATA,'source_plan_sha256':sha(OLD/'status/step2_plan.json'),'script_sha256':sha(Path(__file__)),'model':'stratum-specific REGENIE Step 1 and Step 2','step1_marker_set':'frozen globally pruned 247255 autosomal markers','step2_variant_set':'all 13255796 original quality-controlled variants; stratum minMAC 5','sampling_weights':False,'sex_dosage':'PGEN native male haploid genotypes automatically scaled 0/2 by REGENIE; female 0/1/2','PAR':'hg38 default boundaries; report nonPAR and PAR separately','changes_to_primary_analysis':False})
 status(status='fitting_stratum_models')
 with ThreadPoolExecutor(max_workers=5) as ex:
  for f in as_completed([ex.submit(step1,n) for n in STRATA]):f.result()
 xt=[t for t in TASKS if t['chromosome']=='X'];at=[t for t in TASKS if t['chromosome']!='X']
 status(status='validating_pilot_chunks')
 chunk(xt[0]);chunk(at[0]);status(pilot_chunks_validated=True)
 process_batch(xt,'chromosome_X');summarize(True)
 process_batch(at,'autosomes');summarize(False)
 status(status='complete',completed_utc=now())
except Exception as e:
 status(status='failed',error=str(e),traceback=traceback.format_exc());raise
