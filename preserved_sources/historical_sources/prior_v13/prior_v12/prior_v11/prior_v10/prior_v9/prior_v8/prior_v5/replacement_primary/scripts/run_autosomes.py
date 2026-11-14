from pathlib import Path
import argparse,datetime,hashlib,json,os,shutil,subprocess,time,traceback
import numpy as np,pandas as pd,psutil,zstandard
from google.cloud import storage
N=Path(__file__).resolve().parents[1];C=storage.Client();ENV=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='4')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,d):
 q=p.with_suffix(p.suffix+'.partial');q.write_text(json.dumps(d,indent=2));q.replace(p)
def run(cmd,log):
 save(log.with_suffix('.command.json'),cmd)
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=ENV)
 if r.returncode:raise RuntimeError(str(log)+' command failed '+str(r.returncode))
def step1():
 pid=json.loads((N/'status/step1_pid.json').read_text())['pid']
 while psutil.pid_exists(pid):
  try:
   p=psutil.Process(pid)
   if p.status()==psutil.STATUS_ZOMBIE:break
   if 'regenie' not in p.name():break
  except psutil.NoSuchProcess:break
  save(N/'status/autosomal_status.json',{'status':'waiting_for_strict_cohort_step1','pid':os.getpid()});time.sleep(30)
 log=(N/'step1/model.log').read_text()
 for s in ['number of individuals used in analysis = 386069','5705 cases and 380364 controls']:assert s in log,s
 lf=N/'step1/model_1.loco';keep=pd.read_csv(N/'inputs/keep.tsv',sep='\t',header=None,dtype=str);ids=set(keep[0]+'_'+keep[1])
 with lf.open() as f:
  hh=f.readline().split();assert hh[0]=='FID_IID' and len(hh)==len(ids)+1 and set(hh[1:])==ids
  rows=[]
  for line in f:
   v=line.split();rows.append(int(v[0]));v=np.array(v[1:],float);assert len(v)==len(ids) and np.isfinite(v).all()
 assert rows==list(range(1,24))
 pl=N/'step1/model_pred.list';tokens=pl.read_text().split();assert tokens==['psoriasis_case',str(lf)]
 save(N/'step1/VALIDATED.json',{'status':'validated','n':len(ids),'cases':5705,'controls':380364,'loco_rows':rows,'loco_sha256':sha(lf),'log_sha256':sha(N/'step1/model.log'),'prediction_list':str(pl)})
 return pl

def download(uri,path):
 bucket,name=uri[5:].split('/',1);blob=C.bucket(bucket,user_project=os.environ['GOOGLE_PROJECT']).get_blob(name);assert blob is not None
 if not path.exists():
  temp=Path(str(path)+'.partial');blob.download_to_filename(str(temp),if_generation_match=blob.generation,checksum='crc32c',timeout=600);temp.replace(path)
 assert path.stat().st_size==blob.size
 return {'uri':uri,'generation':blob.generation,'crc32c':blob.crc32c,'bytes':blob.size,'sha256':sha(path)}

def qc(folder):
 ca=pd.read_csv(folder/'cases.acount',sep=r'\s+').set_index('ID');co=pd.read_csv(folder/'controls.acount',sep=r'\s+').set_index('ID')
 cm=pd.read_csv(folder/'cases.vmiss',sep=r'\s+').set_index('ID');om=pd.read_csv(folder/'controls.vmiss',sep=r'\s+').set_index('ID')
 pv=pd.read_csv(folder/'input.pvar',sep='\t').set_index('ID');assert pv.index.is_unique and ca.index.equals(pv.index) and co.index.equals(pv.index) and cm.index.equals(pv.index) and om.index.equals(pv.index)
 for t in [ca,co]:assert t.REF.equals(pv.REF) and t.ALT.equals(pv.ALT) and not t.ALT.str.contains(',').any()
 mac=np.minimum(ca.ALT_CTS+co.ALT_CTS,ca.OBS_CT+co.OBS_CT-ca.ALT_CTS-co.ALT_CTS)
 q=pd.DataFrame({'case_MAC':np.minimum(ca.ALT_CTS,ca.OBS_CT-ca.ALT_CTS),'control_MAC':np.minimum(co.ALT_CTS,co.OBS_CT-co.ALT_CTS),'total_MAC':mac,'MAF':mac/(ca.OBS_CT+co.OBS_CT),'case_callrate':1-cm.MISSING_CT/cm.OBS_CT,'control_callrate':1-om.MISSING_CT/om.OBS_CT})
 for an in ['EUR','AFR','AMR']:
  ht=pd.read_csv(folder/(an+'.controls.hardy'),sep=r'\s+').set_index('ID');assert ht.index.equals(q.index);q[an+'_control_HWE_P']=ht.P
 q['count_call_pass']=(q.case_MAC>=20)&(q.control_MAC>=20)&(q.total_MAC>=100)&(q.MAF>=.01)&(q.case_callrate>=.98)&(q.control_callrate>=.98)
 q['HWE_pass']=q[[an+'_control_HWE_P' for an in ['EUR','AFR','AMR']]].ge(1e-15).all(axis=1)
 q['PASS']=q.count_call_pass & q.HWE_pass
 return pv,q

def process(t,pred):
 key='chr'+t['chromosome']+'_'+t['stem'];out=N/'results/autosomes'/key;out.mkdir(parents=True,exist_ok=True)
 done=out/'COMPLETE.json'
 if done.exists():
  rr=json.loads(done.read_text());assert rr['status']=='validated'
  if rr['passed']:
   assert sha(out/'assoc_psoriasis_case.regenie')==rr['results_sha256']
  return rr
 assert shutil.disk_usage(N).free>128*1024**3,'scratch floor reached'
 scratch=N/'scratch'/key;scratch.mkdir(exist_ok=True);sources={}
 for ext in ['pgen','pvar','psam']:
  path=scratch/('input.'+('pvar.zst' if ext=='pvar' else ext));sources[ext]=download(t[ext],path)
  if ext=='pvar':
   with path.open('rb') as f,(scratch/'input.pvar').open('wb') as g:zstandard.ZstdDecompressor().copy_stream(f,g)
 save(out/'source_receipts.json',sources)
 for group in ['cases','controls','EUR.controls','AFR.controls','AMR.controls']:
  cmd=['plink2','--pfile',str(scratch/'input'),'--keep',str(N/'inputs'/(group+'.keep.tsv')),'--threads','4','--memory','4096','--out',str(scratch/group)]
  cmd+=['--freq','counts','--missing','variant-only'] if group in ['cases','controls'] else ['--hardy']
  run(cmd,scratch/(group+'.console.log'))
 pv,q=qc(scratch);q.to_csv(out/'variant_QC.private.tsv.gz',sep='\t');selected=pv.loc[q.PASS];selected.index.to_series().to_csv(out/'passed_variants.txt',index=False,header=False)
 result={'status':'validated','task':key,'source_variants':len(pv),'passed':len(selected),'count_call_pass':int(q.count_call_pass.sum()),'HWE_pass':int(q.HWE_pass.sum()),'n':386069,'cases':5705,'input_kind':'historical fixed-panel autosomal checkpoints re-QC in strict cohort','updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 if len(selected):
  cmd=['regenie','--step','2','--bt','--pgen',str(scratch/'input'),'--keep',str(N/'inputs/keep.tsv'),'--extract',str(out/'passed_variants.txt'),'--phenoFile',str(N/'inputs/phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(N/'inputs/covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(pred),'--af-cc','--threads','4','--out',str(out/'assoc')]
  run(cmd,out/'assoc.console.log')
  res=pd.read_csv(out/'assoc_psoriasis_case.regenie',sep=r'\s+');assert res.ID.is_unique and len(res)==len(selected) and set(res.ID)==set(selected.index) and res.TEST.eq('ADD').all()
  matched=selected.reindex(res.ID)
  assert np.array_equal(res.ALLELE0.to_numpy(str),matched.REF.to_numpy(str)) and np.array_equal(res.ALLELE1.to_numpy(str),matched.ALT.to_numpy(str)) and np.array_equal(res.GENPOS.to_numpy(int),matched.POS.to_numpy(int))
  for term in ['number of individuals used in analysis = 386069','5705 cases and 380364 controls','LOCO predictions']:assert term in (out/'assoc.log').read_text(),term
  estimates=res[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');valid=np.isfinite(estimates).all(axis=1);assert estimates.loc[valid,'SE'].gt(0).all() and estimates.loc[valid,'LOG10P'].ge(0).all();assert res.N.le(386069).all()
  result.update(results_sha256=sha(out/'assoc_psoriasis_case.regenie'),estimated=int(valid.sum()),failed_or_unavailable=int((~valid).sum()),effect_allele='ALLELE1 = PVAR ALT; per coded dosage unit')
 for p in scratch.iterdir():
  if p.suffix in ['.log','.json']:shutil.copy2(p,out/p.name)
 save(done,result)
 shutil.rmtree(scratch)
 return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int);a=parser.parse_args();done=[]
 try:
  pred=step1();tasks=json.loads((N/'autosomal_tasks.json').read_text());assert len(tasks)==1104
  if a.limit:tasks=tasks[:a.limit]
  for t in tasks:
   save(N/'status/autosomal_status.json',{'status':'running','completed':len(done),'expected':len(tasks),'current':'chr'+t['chromosome']+'_'+t['stem'],'pid':os.getpid()})
   rr=process(t,pred);done.append(rr);print('PSORIASIS_AUTOSOME_CHUNK_COMPLETE',rr['task'],rr['passed'],flush=True)
  save(N/'status/autosomal_status.json',{'status':'complete','completed':len(done),'expected':len(tasks),'source_variants':sum(r['source_variants'] for r in done),'passed':sum(r['passed'] for r in done),'scope':'autosomal fixed panel only; X and final diagnostics separate'})
 except Exception as e:
  save(N/'status/autosomal_status.json',{'status':'failed','completed':len(done),'error':str(e),'traceback':traceback.format_exc()});raise
