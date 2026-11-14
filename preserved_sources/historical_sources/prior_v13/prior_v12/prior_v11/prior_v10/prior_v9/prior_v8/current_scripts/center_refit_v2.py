from pathlib import Path
import json,hashlib,sys,os,traceback,subprocess,datetime
import pandas as pd,numpy as np
from threadpoolctl import threadpool_limits
_limit=threadpool_limits(limits=1)
N=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');I=N/'reconciliation_20260912';Q=N.parent/'sample_lead_qc_20260910';O=I/'center_refit_v2';O.mkdir(exist_ok=True)
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(name,v):(O/name).write_text(json.dumps(v,indent=2,default=str)+chr(10))
def run(cmd,prefix):
 save(prefix+'.command.json',cmd)
 with (O/(prefix+'.console.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=ENV)
 assert r.returncode==0,(prefix,r.returncode)
def main():
 c=pd.read_csv(N/'inputs/covariates.tsv',sep=chr(9),dtype={'FID':str,'IID':str});d=pd.read_csv(Q/'sample_audit.private.tsv',sep=chr(9),dtype={'IID':str}).set_index('IID');assert c.IID.is_unique and d.index.is_unique;matched=d.loc[c.IID];assert matched.strict_eligible.all() and len(matched)==int(d.strict_eligible.sum());assert np.array_equal(c.sex,matched.SEX.to_numpy()==1);assert matched.site_id.isin(['bcm','bi','uw']).all()
 for center in ['bi','uw']:c['center_'+center]=matched.site_id.eq(center).to_numpy(int)
 assert np.isfinite(c.drop(columns=['FID','IID']).to_numpy()).all();v=c.drop(columns=['FID','IID']).to_numpy(float);assert np.linalg.matrix_rank(np.column_stack([np.ones(len(v)),v]))==v.shape[1]+1
 cf=O/'center_covariates.private.tsv';c.to_csv(cf,sep=chr(9),index=False)
 cmd=json.loads((N/'step1/command.json').read_text())
 old=N.parent/'stratified_calibration';panel=old/'private_panel';status=json.loads((old/'panel_status.json').read_text());local=O/'training_panel';local.mkdir(exist_ok=True);prefix=local/'calibration_pruned'
 if not Path(str(prefix)+'.bed').exists():
  import shutil
  assert shutil.disk_usage(O).free>100*1024**3
  print('RESTORE_VERIFIED_TRAINING_PANEL_START',flush=True)
  with Path(str(prefix)+'.bed').open('wb') as out,Path(str(prefix)+'.bim').open('w') as bim:
   out.write(bytes([108,27,1]))
   for rec in sorted(status['chromosomes'],key=lambda x:int(x['chromosome'])):
    stem=panel/('chr'+rec['chromosome']);bed=Path(str(stem)+'.bed');assert sha(bed)==rec['bed_sha256']
    with bed.open('rb') as f:
     assert f.read(3)==bytes([108,27,1]);shutil.copyfileobj(f,out,8*1024**2)
    bim.write(Path(str(stem)+'.bim').read_text())
  shutil.copyfile(panel/'chr1.fam',Path(str(prefix)+'.fam'))
 assert sha(Path(str(prefix)+'.bed'))==status['bed_sha256']
 for ext in ['bim','fam']:assert sha(Path(str(prefix)+'.'+ext))==sha(panel/('calibration_pruned.'+ext))
 save('training_panel_verification.private.json',dict(bed_sha256=status['bed_sha256'],per_chromosome_hashes_verified=True,bim_fam_match=True))
 cmd[cmd.index('--bed')+1]=str(prefix)
 for flag,value in [('--covarFile',str(cf)),('--lowmem-prefix',str(O/'step1_lowmem')),('--out',str(O/'step1')),('--threads','4')]:cmd[cmd.index(flag)+1]=value
 save('PROTOCOL.private.json',dict(analysis='Sequencing-center sensitivity in the exact corrected primary cohort, both dosage assumptions; two center indicators with bcm reference added in both REGENIE Step 1 and Step 2.',covariates=c.columns.tolist(),relatedness='REGENIE whole-genome prediction from the same fixed LD-pruned training panel; center-adjusted Step 1 refit.',lead_status='NLGN4X QC passing; the other three historical X leads are QC-excluded diagnostics only.',fine_batch='Flowcell/library/run not present in reviewed released metadata; center is not a substitute.',release='HELD'))
 print('CENTER_STEP1_REFIT_START',flush=True);run(cmd,'step1');print('CENTER_STEP1_REFIT_FINISHED',flush=True)
 log=(O/'step1.log').read_text();assert 'number of individuals used in analysis = 386069' in log and '5705 cases and 380364 controls' in log
 with (O/'step1_1.loco').open() as f:
  h=f.readline().split();assert h[0]=='FID_IID' and set(h[1:])==set(c.FID+'_'+c.IID);rows=[]
  for line in f:
   t=line.split();rows.append(int(t[0]));assert len(t)==len(c)+1 and np.isfinite(np.array(t[1:],float)).all()
 assert rows==list(range(1,24))
 pred=O/'step1_pred.list';assert pred.read_text().split()==['psoriasis_case',str(O/'step1_1.loco')]
 leads=['chrX:5998724:A:G','chrX:22917619:C:T','chrX:67736589:T:C','chrX:67842269:G:T'];(O/'leads.txt').write_text(chr(10).join(leads)+chr(10));pg=N/'X_projection_pilot/chunk_0000/input';pv=pd.read_csv(str(pg)+'.pvar',sep=chr(9)).set_index('ID');assert set(leads)<=set(pv.index)
 exact=pd.read_csv(I/'estimate_audit/exact_association_rows.private.tsv.gz',sep=chr(9));results=[]
 for adjustment in ['matched','center_refit']:
  for model in ['male02','male01']:
   name=adjustment+'_'+model
   cmd=['regenie','--step','2','--bt','--pgen',str(pg),'--keep',str(N/'inputs/keep.tsv'),'--extract',str(O/'leads.txt'),'--phenoFile',str(N/'inputs/phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(cf if adjustment=='center_refit' else N/'inputs/covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(pred if adjustment=='center_refit' else N/'step1/model_pred.list'),'--af-cc','--par-region','hg38','--threads','4','--out',str(O/name)]
   if model=='male01':cmd+=['--skip-dosage-comp']
   run(cmd,name);t=pd.read_csv(O/(name+'_psoriasis_case.regenie'),sep=chr(92)+'s+');assert set(t.ID)==set(leads) and t.TEST.eq('ADD').all();assert np.isfinite(t[['BETA','SE','LOG10P']]).all().all() and t.SE.gt(0).all()
   p=pv.loc[t.ID];assert np.array_equal(t.ALLELE0,p.REF) and np.array_equal(t.ALLELE1,p.ALT)
   if adjustment=='matched':
    e=exact.loc[exact.group.eq('pooled')&exact.domain.eq('X')&exact.model.eq(model)&exact.ID.eq(leads[0])];assert len(e)==1;z=t.loc[t.ID.eq(leads[0])];assert np.allclose(z[['BETA','SE','LOG10P','N']].to_numpy(float),e[['BETA','SE','LOG10P','N']].to_numpy(float),rtol=2e-5,atol=1e-7),'primary NLGN4X reproduction'
   t['adjustment']=adjustment;t['model']=model;t['primary_QC_PASS']=t.ID.eq(leads[0]);t['OR']=np.exp(t.BETA);t['CI95_low']=np.exp(t.BETA-1.959963984540054*t.SE);t['CI95_high']=np.exp(t.BETA+1.959963984540054*t.SE);t['P_from_REGENIE_LOG10P']=10.**(-t.LOG10P);t['source_file']=str(O/(name+'_psoriasis_case.regenie'));t['source_sha256']=sha(O/(name+'_psoriasis_case.regenie'));results.append(t);print('CENTER_MODEL_COMPLETE',name,flush=True)
 out=pd.concat(results,ignore_index=True);out.to_csv(O/'estimates.private.tsv',sep=chr(9),index=False);m=out[out.adjustment.eq('matched')].merge(out[out.adjustment.eq('center_refit')],on=['ID','model'],suffixes=('_matched','_center'));assert np.array_equal(m.N_matched,m.N_center);m['relative_OR_change']=m.OR_center/m.OR_matched-1;m.to_csv(O/'paired_comparison.private.tsv',sep=chr(9),index=False)
 save('COMPLETE.private.json',dict(status='verified',exact_primary_cohort=True,center_refit_both_stages=True,matched_primary_NLGN4X_reproduced=True,models=len(out),loco_rows=rows,loco_sha256=sha(O/'step1_1.loco'),covariates_sha256=sha(cf),script_sha256=sha(Path(__file__)),release='HELD'));print('CENTER_REFIT_COMPLETE',flush=True)
try:main()
except Exception as e:save('FAILED.private.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
