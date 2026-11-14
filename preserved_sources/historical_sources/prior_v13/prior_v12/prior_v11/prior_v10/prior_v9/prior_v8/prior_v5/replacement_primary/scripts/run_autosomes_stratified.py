from pathlib import Path
import argparse,concurrent.futures,datetime,json,multiprocessing,os,shutil,time,traceback
import numpy as np,pandas as pd,psutil,zstandard
from run_autosomes import download,run,sha,save
from run_stratified_X import validate_prediction
N=Path(__file__).resolve().parents[1]
LABELS=['AMR','AFR','female','male','EUR']
ROOT=N/'results/autosomes_stratified'
SCRATCH=N/'scratch_autosomes_stratified'
CFG=json.loads((N/'stratified_step1/strata.json').read_text())
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def preflight():
 out={}
 for label in LABELS:
  inp=N/'stratified_step1'/label;pred,meta=validate_prediction(label)
  ph=pd.read_csv(inp/'phenotype.tsv',sep='\t',dtype={'FID':str,'IID':str})
  cv=pd.read_csv(inp/'covariates.tsv',sep='\t',dtype={'FID':str,'IID':str})
  kp=pd.read_csv(inp/'keep.tsv',sep='\t',header=None,dtype=str)
  assert len(ph)==meta['n']==CFG[label]['n'] and not ph[['FID','IID']].duplicated().any()
  assert list(cv.columns[2:])==CFG[label]['covariates']
  assert ph[['FID','IID']].equals(cv[['FID','IID']])
  assert set(map(tuple,kp.values))==set(map(tuple,ph[['FID','IID']].values))
  assert set(ph.psoriasis_case)=={0,1} and int(ph.psoriasis_case.sum())==meta['cases']
  if label in ['female','male']:assert not {'sex','age_sex'}.intersection(cv.columns)
  for case,name in [(1,'cases'),(0,'controls')]:
   keep=inp/(name+'.autosomal.keep.tsv');data=ph.loc[ph.psoriasis_case.eq(case),['FID','IID']]
   if keep.exists():assert pd.read_csv(keep,sep='\t',header=None,dtype=str).values.tolist()==data.values.tolist()
   else:data.to_csv(keep,sep='\t',header=False,index=False)
  out[label]={'n':meta['n'],'cases':meta['cases'],'controls':CFG[label]['controls'],'covariates':CFG[label]['covariates'],'prediction_sha256':meta['loco_sha256'],'keep_sha256':sha(inp/'keep.tsv'),'case_keep_sha256':sha(inp/'cases.autosomal.keep.tsv'),'control_keep_sha256':sha(inp/'controls.autosomal.keep.tsv')}
 save(N/'status/autosomal_stratified_preflight.json',{'status':'validated','strata':out,'script_sha256':sha(Path(__file__)),'utc':stamp()})
 return out

def parent_ready(key):
 p=N/'results/autosomes'/key/'COMPLETE.json'
 while not p.exists():
  s=json.loads((N/'status/autosomal_status.json').read_text())
  if s['status']=='failed':raise RuntimeError('Pooled autosomal parent failed')
  if s['status']=='complete':raise RuntimeError('Pooled parent marked complete but chunk missing '+key)
  pid=s.get('pid')
  if pid and (not psutil.pid_exists(pid) or psutil.Process(pid).status()==psutil.STATUS_ZOMBIE):raise RuntimeError('Pooled parent process unavailable before chunk '+key)
  time.sleep(30)
 meta=json.loads(p.read_text());assert meta['status']=='validated';return p.parent,meta

def process(t):
 key='chr'+t['chromosome']+'_'+t['stem'];parent,pm=parent_ready(key)
 chunk=ROOT/key;chunk.mkdir(parents=True,exist_ok=True);complete=chunk/'COMPLETE.json'
 if complete.exists():
  rr=json.loads(complete.read_text());assert rr['status']=='validated' and rr['pooled_complete_sha256']==sha(parent/'COMPLETE.json')
  for label,m in rr['strata'].items():
   if m['passed']:assert sha(chunk/label/'assoc_psoriasis_case.regenie')==m['results_sha256']
  return rr
 assert shutil.disk_usage(N).free>128*1024**3,'Scratch free-space floor reached'
 folder=SCRATCH/key;folder.mkdir(parents=True,exist_ok=True)
 sources={};frozen=json.loads((parent/'source_receipts.json').read_text())
 for ext in ['pgen','pvar','psam']:
  path=folder/('input.'+('pvar.zst' if ext=='pvar' else ext));sources[ext]=download(t[ext],path)
  for k in ['generation','bytes','sha256']:assert sources[ext][k]==frozen[ext][k],(key,ext,k,'source changed')
  if ext=='pvar':
   with path.open('rb') as f,(folder/'input.pvar').open('wb') as g:zstandard.ZstdDecompressor().copy_stream(f,g)
 save(chunk/'source_receipts.json',sources)
 pv=pd.read_csv(folder/'input.pvar',sep='\t').set_index('ID')
 pq=pd.read_csv(parent/'variant_QC.private.tsv.gz',sep='\t').set_index('ID')
 assert pv.index.is_unique and pv.index.equals(pq.index) and pq.PASS.dtype.kind=='b'
 assert int(pq.PASS.sum())==pm['passed']
 rr={'status':'validated','task':key,'source_variants':len(pv),'pooled_passed':pm['passed'],'pooled_complete_sha256':sha(parent/'COMPLETE.json'),'strata':{},'source':'historical fixed autosomal panel, re-QC in strict cohort, then subgroup count/call/MAF gates','utc':stamp()}
 for label in LABELS:
  inp=N/'stratified_step1'/label;out=chunk/label;out.mkdir(exist_ok=True)
  proof=out/'COMPLETE.json'
  if proof.exists():
   old=json.loads(proof.read_text());assert old['status']=='validated'
   if old['passed']:assert sha(out/'assoc_psoriasis_case.regenie')==old['results_sha256']
   rr['strata'][label]=old;continue
  for group in ['cases','controls']:
   prefix=out/group
   run(['plink2','--pfile',str(folder/'input'),'--keep',str(inp/(group+'.autosomal.keep.tsv')),'--freq','counts','--missing','variant-only','--threads','4','--memory','4096','--out',str(prefix)],out/(group+'.console.log'))
  ca=pd.read_csv(out/'cases.acount',sep=r'\s+').set_index('ID');co=pd.read_csv(out/'controls.acount',sep=r'\s+').set_index('ID')
  cm=pd.read_csv(out/'cases.vmiss',sep=r'\s+').set_index('ID');om=pd.read_csv(out/'controls.vmiss',sep=r'\s+').set_index('ID')
  for d in [ca,co,cm,om]:assert d.index.equals(pv.index)
  for d in [ca,co]:assert d.REF.equals(pv.REF) and d.ALT.equals(pv.ALT) and not d.ALT.str.contains(',').any()
  for ac,mis,n in [(ca,cm,CFG[label]['cases']),(co,om,CFG[label]['controls'])]:
   assert mis.OBS_CT.eq(n).all() and ((mis.OBS_CT-mis.MISSING_CT)*2).equals(ac.OBS_CT)
   assert ac.ALT_CTS.between(0,ac.OBS_CT).all()
  q=pd.DataFrame(index=pv.index)
  q['case_AC']=ca.ALT_CTS;q['case_AN']=ca.OBS_CT;q['control_AC']=co.ALT_CTS;q['control_AN']=co.OBS_CT
  q['case_MAC']=np.minimum(ca.ALT_CTS,ca.OBS_CT-ca.ALT_CTS);q['control_MAC']=np.minimum(co.ALT_CTS,co.OBS_CT-co.ALT_CTS)
  q['total_MAC']=np.minimum(ca.ALT_CTS+co.ALT_CTS,ca.OBS_CT+co.OBS_CT-ca.ALT_CTS-co.ALT_CTS)
  q['MAF']=q.total_MAC/(ca.OBS_CT+co.OBS_CT)
  q['case_called']=cm.OBS_CT-cm.MISSING_CT;q['control_called']=om.OBS_CT-om.MISSING_CT
  q['case_callrate']=q.case_called/CFG[label]['cases'];q['control_callrate']=q.control_called/CFG[label]['controls']
  q['pooled_PASS']=pq.PASS
  q['subgroup_PASS']=(q.case_MAC>=20)&(q.control_MAC>=20)&(q.total_MAC>=100)&(q.MAF>=.01)&(q.case_callrate>=.98)&(q.control_callrate>=.98)
  q['PASS']=q.pooled_PASS&q.subgroup_PASS
  q.to_csv(out/'variant_QC.private.tsv.gz',sep='\t')
  selected=pv.loc[q.PASS];selected.index.to_series().to_csv(out/'passed_variants.txt',header=False,index=False)
  m={'status':'validated','stratum':label,'passed':len(selected),'estimated':0,'failed_or_unavailable':0,'n_total_stratum':CFG[label]['n'],'n_total_primary_cohort':386069,'cases':CFG[label]['cases'],'controls':CFG[label]['controls'],'QC_sha256':sha(out/'variant_QC.private.tsv.gz'),'effect_allele':'ALLELE1 = PVAR ALT; OR per one diploid coded dosage unit','CI':'same-model BETA +/- 1.959963984540054 SE on log scale','P':'REGENIE LOG10P; Firth likelihood-ratio where invoked, not Wald reconstructed P'}
  if len(selected):
   run(['regenie','--step','2','--bt','--pgen',str(folder/'input'),'--keep',str(inp/'keep.tsv'),'--extract',str(out/'passed_variants.txt'),'--phenoFile',str(inp/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(inp/'covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(inp/'model_pred.list'),'--af-cc','--threads','4','--out',str(out/'assoc')],out/'assoc.console.log')
   rp=out/'assoc_psoriasis_case.regenie';d=pd.read_csv(rp,sep=r'\s+');assert d.ID.is_unique and len(d)==len(selected) and set(d.ID)==set(selected.index) and d.TEST.eq('ADD').all()
   match=selected.reindex(d.ID)
   assert np.array_equal(d.ALLELE0.to_numpy(str),match.REF.to_numpy(str)) and np.array_equal(d.ALLELE1.to_numpy(str),match.ALT.to_numpy(str)) and np.array_equal(d.GENPOS.to_numpy(int),match.POS.to_numpy(int))
   assert d.CHROM.astype(str).eq(str(t['chromosome'])).all() and d.N.between(1,CFG[label]['n']).all()
   log=(out/'assoc.log').read_text()
   for term in ['number of individuals used in analysis = '+str(CFG[label]['n']),str(CFG[label]['cases'])+' cases and '+str(CFG[label]['controls'])+' controls','LOCO predictions']:assert term in log,term
   values=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');valid=np.isfinite(values).all(axis=1)&values.SE.gt(0)&values.LOG10P.ge(0)
   m.update(estimated=int(valid.sum()),failed_or_unavailable=int((~valid).sum()),results_sha256=sha(rp),threshold_counts={str(t):int((valid & values.LOG10P.gt(-np.log10(t))).sum()) for t in [5e-8,3e-8,5e-9,1e-9]})
  save(proof,m);rr['strata'][label]=m
 for group in ['cases','controls']:
  total=None
  for label in ['female','male']:
   d=pd.read_csv(chunk/label/(group+'.acount'),sep=r'\s+').set_index('ID')[['ALT_CTS','OBS_CT']]
   total=d if total is None else total+d
  assert np.array_equal(np.minimum(total.ALT_CTS,total.OBS_CT-total.ALT_CTS).to_numpy(),pq[group[:-1]+'_MAC'].to_numpy()),'Sex subgroup count partition disagrees with pooled QC'
 rr['sex_partition_count_validation']='all candidate variants matched pooled MAC'
 save(complete,rr);shutil.rmtree(folder);return rr

def calculate_all(tasks):
 from finalize_pooled_calibration import quantiles
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 diag=N/'diagnostics/autosomes_stratified_complete';diag.mkdir(parents=True,exist_ok=True)
 all_reports={}
 for label in LABELS:
  parts=[];nonmhc=[];seen=set();bad=0;total=0;source=[]
  for t in tasks:
   key='chr'+t['chromosome']+'_'+t['stem'];out=ROOT/key/label;m=json.loads((out/'COMPLETE.json').read_text());assert m['status']=='validated'
   total+=m['passed'];bad+=m['failed_or_unavailable']
   if not m['passed']:continue
   rp=out/'assoc_psoriasis_case.regenie';assert sha(rp)==m['results_sha256']
   d=pd.read_csv(rp,sep=r'\s+');assert not seen.intersection(d.ID);seen.update(d.ID)
   v=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');ok=np.isfinite(v).all(axis=1)&v.SE.gt(0)&v.LOG10P.ge(0)
   assert int(ok.sum())==m['estimated'];parts.append(v.loc[ok,'LOG10P'].to_numpy(float))
   mhc=d.CHROM.eq(6)&d.GENPOS.between(25000000,34000000)
   nonmhc.append(v.loc[ok&~mhc,'LOG10P'].to_numpy(float));source.append({'path':str(rp),'sha256':m['results_sha256'],'tested':m['passed'],'estimated':m['estimated'],'unavailable':m['failed_or_unavailable']})
  del seen
  for suffix,arrays in [('all',parts),('MHC_excluded',nonmhc)]:
   vals=np.concatenate(arrays);lam,qq=quantiles(vals);out=diag/label;out.mkdir(exist_ok=True)
   qq.to_csv(out/(suffix+'_QQ_full_ranks.tsv.gz'),sep='\t',index=False)
   report={'stratum':label,'subset':suffix,'lambda_GC_from_full_output_P':lam,'valid_tests':len(vals),'total_tested_before_MHC_exclusion':total,'failed_or_unavailable_before_MHC_exclusion':bad,'weights':'unit weights for all full-output tests, for both lambda and QQ ranks','QQ':'full ranks; deterministic display thinning never reranks; pointwise 95% beta order-statistic band','MHC_interval_GRCh38':'chr6:25000000-34000000 inclusive','scope':'fixed autosomal panel with pooled QC and subgroup count/call/MAF gates; within-cohort calibration, not replication','threshold_counts':{str(t):int((vals>-np.log10(t)).sum()) for t in [5e-8,3e-8,5e-9,1e-9]}}
   save(out/(suffix+'_calibration.json'),report)
   fig,ax=plt.subplots(figsize=(5.2,5.2));ax.fill_between(qq.expected_log10P,qq.null_95_low,qq.null_95_high,color='#d4d8df',alpha=.65,linewidth=0);ax.scatter(qq.expected_log10P,qq.observed_log10P,s=3,color='#334e68',rasterized=True,linewidths=0)
   xmax=float(qq.expected_log10P.max())*1.025;ymax=max(xmax,float(qq.observed_log10P.max())*1.025);ax.plot([0,xmax],[0,xmax],color='#8c939d',linewidth=.8);ax.set(xlim=(0,xmax),ylim=(0,ymax),xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title='Autosomes: '+label+(' (MHC excluded)' if suffix!='all' else ''))
   ax.text(.04,.96,'N = '+format(len(vals),',')+'\nLambda GC = '+format(lam,'.3f'),transform=ax.transAxes,va='top',fontsize=10);ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/(suffix+'_QQ.png'),dpi=300);fig.savefig(out/(suffix+'_QQ.pdf'));plt.close(fig)
   all_reports[label+'_'+suffix]=report
  save(diag/label/'source_receipts.json',source)
 save(diag/'COMPLETE.json',{'status':'complete','strata':LABELS,'chunks':len(tasks),'reports':all_reports,'script_sha256':sha(Path(__file__)),'utc':stamp()})
 return all_reports

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--pilot-only',action='store_true');args=parser.parse_args();assert 1<=args.workers<=4
 ROOT.mkdir(parents=True,exist_ok=True);SCRATCH.mkdir(parents=True,exist_ok=True);preflight()
 tasks=json.loads((N/'autosomal_tasks.json').read_text());assert len(tasks)==1104
 status=N/'status/autosomal_stratified_status.json'
 for t in tasks[:2]:
  rr=process(t);assert all(m['failed_or_unavailable']==0 for m in rr['strata'].values()),'Pilot has unavailable estimates; requires review'
 save(ROOT/'PILOT_VALIDATED.json',{'status':'validated','chunks':2,'strata':LABELS,'all_available':True,'sex_partition_counts':'matched pooled MAC for every pilot variant','utc':stamp()})
 if args.pilot_only:return
 done=0;it=iter(tasks)
 ctx=multiprocessing.get_context('spawn')
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,mp_context=ctx) as pool:
  active={}
  for _ in range(args.workers):
   t=next(it,None)
   if t:active[pool.submit(process,t)]=t
  while active:
   completed,_=concurrent.futures.wait(active,return_when=concurrent.futures.FIRST_COMPLETED)
   for f in completed:
    active.pop(f);rr=f.result();done+=1
    save(status,{'status':'running','completed_chunks':done,'expected_chunks':len(tasks),'strata':LABELS,'last_completed':rr['task'],'pid':os.getpid(),'utc':stamp()})
    print('AUTOSOMAL_STRATIFIED_CHUNK_COMPLETE',rr['task'],done,flush=True)
    t=next(it,None)
    if t:active[pool.submit(process,t)]=t
 save(ROOT/'COMPLETE.json',{'status':'complete','chunks':len(tasks),'strata':LABELS,'source_scope':'fixed autosomal panel; not fresh source-wide variant discovery','script_sha256':sha(Path(__file__)),'utc':stamp()})
 save(status,{'status':'calculating_full_output_diagnostics','completed_chunks':len(tasks),'pid':os.getpid(),'utc':stamp()})
 calculate_all(tasks);save(status,{'status':'complete','completed_chunks':len(tasks),'strata':LABELS,'utc':stamp()})
if __name__=='__main__':
 try:main()
 except Exception as e:
  save(N/'status/autosomal_stratified_FAILED.json',{'status':'failed','error':str(e),'traceback':traceback.format_exc(),'utc':stamp()});raise
