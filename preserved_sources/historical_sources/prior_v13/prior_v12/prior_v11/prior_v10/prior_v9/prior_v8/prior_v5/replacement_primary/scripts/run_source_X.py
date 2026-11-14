from pathlib import Path
import hashlib,json,os,subprocess,time,traceback
import numpy as np,pandas as pd,psutil
from reconcile_X_site_filters import sha,save
N=Path(__file__).resolve().parents[1];ENV=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='4')
def wait_for(path,pidfile,failures):
 while not path.exists():
  for fail in failures:
   if fail.exists():raise RuntimeError(str(fail)+' requires review')
  pid=json.loads(pidfile.read_text())['pid']
  try:alive=psutil.pid_exists(pid) and psutil.Process(pid).status()!=psutil.STATUS_ZOMBIE
  except psutil.NoSuchProcess:alive=False
  if not alive:raise RuntimeError(str(pidfile)+' ended before required completion receipt '+str(path))
  time.sleep(30)
 return json.loads(path.read_text())
def run(cmd,path):
 save(path.with_suffix('.command.json'),cmd)
 with path.open('w') as handle:r=subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT,env=ENV)
 if r.returncode:raise RuntimeError('Command failed; inspect '+str(path))
def hwe_pilot(folder):
 q=pd.read_csv(folder/'variant_QC.private.tsv.gz',sep=chr(9)).set_index('ID')
 for an in ['EUR','AFR','AMR']:
  cmd=['plink2','--pfile',str(folder/'input'),'--split-par','hg38','--keep',str(N/'inputs'/(an+'.female_controls.keep.tsv')),'--hardy','--threads','4','--memory','4096','--out',str(folder/(an+'.female_controls'))]
  run(cmd,folder/(an+'.hwe.console.log'))
  ht=pd.concat([pd.read_csv(p,sep=chr(92)+'s+').set_index('ID') for p in [folder/(an+'.female_controls.hardy'),folder/(an+'.female_controls.hardy.x')] if p.exists()])
  assert ht.index.is_unique and set(ht.index)==set(q.index)
  q[an+'_female_control_HWE_P']=ht.P.reindex(q.index)
 q['HWE_PASS']=q[[an+'_female_control_HWE_P' for an in ['EUR','AFR','AMR']]].ge(1e-15).all(axis=1)
 q['PASS_except_site_FILTER']=q.count_call_maf_PASS&q.HWE_PASS
 q.to_csv(folder/'QC_with_HWE.private.tsv.gz',sep=chr(9));save(folder/'HWE_COMPLETE.json',{'variants':len(q),'passed_except_site_FILTER':int(q.PASS_except_site_FILTER.sum()),'site_FILTER_pending':True})
def process(folder,outroot,sites,pred):
 out=outroot/folder.name;out.mkdir(parents=True,exist_ok=True)
 complete=out/'COMPLETE.json'
 if complete.exists():
  receipt=json.loads(complete.read_text());assert receipt['status']=='validated'
  for s in receipt['models']:assert sha(out/(s+'_psoriasis_case.regenie'))==receipt['models'][s]['results_sha256']
  return receipt
 assert (folder/'HWE_COMPLETE.json').exists()
 pv=pd.read_csv(folder/'input.pvar',sep=chr(9)).set_index('ID');q=pd.read_csv(folder/'QC_with_HWE.private.tsv.gz',sep=chr(9)).set_index('ID');assert pv.index.equals(q.index) and q.index.is_unique
 source=sites.reindex(q.source_index);assert source.site_mapping.eq('matched').all();assert source.site_FILTER_PASS.dtype.kind=='b' and source.ALT_order_differs.dtype.kind=='b' and q.PASS_except_site_FILTER.dtype.kind=='b'
 assert np.array_equal(pv.REF.to_numpy(str),source.REF.to_numpy(str)) and np.array_equal(pv.POS.to_numpy(int),source.POS.to_numpy(int))
 for pr,qr,sr in zip(pv.itertuples(),q.itertuples(),source.itertuples()):
  assert pr.ALT==sr.ALT.split(',')[int(qr.allele_index)-1]
  assert pr.Index=='chrX:'+str(int(pr.POS))+':'+pr.REF+':'+pr.ALT
 q['site_FILTER_PASS']=source.site_FILTER_PASS.to_numpy(bool);q['site_filters']=source.site_filters.to_numpy(str)
 q['site_ALT_order_reconciled']=source.ALT_order_differs.to_numpy(bool);q['PASS']=q.PASS_except_site_FILTER&q.site_FILTER_PASS
 q.to_csv(out/'final_variant_QC.private.tsv.gz',sep=chr(9));selected=pv.loc[q.PASS]
 selected.index.to_series().to_csv(out/'passed_variants.txt',header=False,index=False)
 ps=pd.read_csv(folder/'input.psam',sep=chr(9),dtype={'#FID':str,'IID':str})
 ph=pd.read_csv(N/'inputs/phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str})
 assert len(ps)==len(ph)==386069 and set(ps.IID)==set(ph.IID) and ps.SEX.isin([1,2]).all()
 receipt={'status':'validated','candidate_alleles':len(q),'QC_passed':len(selected),'models':{},'n_total_cohort':386069,'n_cases_total':5705,'source':'unaffected provider PGEN with exact site FILTER and genotype count/call/MAF/female-control HWE QC','PGEN_sha256':sha(folder/'input.pgen'),'PVAR_sha256':sha(folder/'input.pvar'),'PSAM_sha256':sha(folder/'input.psam'),'site_reconciliation_sha256':sha(N/'source_X_sites/RECONCILIATION_COMPLETE.json'),'script_sha256':sha(Path(__file__))}
 for scale in ['02','01']:
  if not len(selected):continue
  prefix='male'+scale
  cmd=['regenie','--step','2','--bt','--pgen',str(folder/'input'),'--keep',str(N/'inputs/keep.tsv'),'--extract',str(out/'passed_variants.txt'),'--phenoFile',str(N/'inputs/phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(N/'inputs/covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(pred),'--af-cc','--par-region','hg38','--threads','4','--out',str(out/prefix)]
  if scale=='01':cmd+=['--skip-dosage-comp']
  run(cmd,out/(prefix+'.console.log'))
  resultpath=out/(prefix+'_psoriasis_case.regenie');res=pd.read_csv(resultpath,sep=chr(92)+'s+')
  assert res.ID.is_unique and set(res.ID)==set(selected.index) and len(res)==len(selected) and res.TEST.eq('ADD').all()
  match=selected.reindex(res.ID)
  assert np.array_equal(res.ALLELE0.to_numpy(str),match.REF.to_numpy(str)) and np.array_equal(res.ALLELE1.to_numpy(str),match.ALT.to_numpy(str)) and np.array_equal(res.GENPOS.to_numpy(int),match.POS.to_numpy(int))
  log=(out/(prefix+'.log')).read_text()
  for term in ['number of individuals used in analysis = 386069','5705 cases and 380364 controls','LOCO predictions']:assert term in log,term
  assert res.N.between(1,386069).all()
  values=res[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');valid=np.isfinite(values).all(axis=1)&values.SE.gt(0)&values.LOG10P.ge(0)
  derived=res.copy();derived['model']='female012_male'+scale;derived['estimate_available']=valid
  derived['OR']=np.exp(values.BETA);derived['CI95_Wald_low']=np.exp(values.BETA-1.959963984540054*values.SE);derived['CI95_Wald_high']=np.exp(values.BETA+1.959963984540054*values.SE)
  derived['P_from_REGENIE_LOG10P']=np.power(10.,-values.LOG10P);derived['N_total_cohort']=386069
  qc=q.reindex(res.ID)
  for name in ['case_called','control_called','case_AC','case_AN','control_AC','control_AN','nonPAR']:
   derived[name]=qc[name].to_numpy()
  derived['N_nonmissing_genotype']=derived.case_called+derived.control_called
  derived.to_csv(out/(prefix+'.derived.private.tsv.gz'),sep=chr(9),index=False)
  receipt['models'][prefix]={'estimated':int(valid.sum()),'failed_or_unavailable':int((~valid).sum()),'results_sha256':sha(resultpath),'effect_allele':'ALLELE1 = selected provider ALT; other alleles collapsed to dosage zero','dosage_units':'female 0/1/2; non-PAR male 0/'+scale[-1]+'; PAR diploid','N_definition':'N is reported by REGENIE; total cohort and observed nonmissing genotype counts retained separately','CI_definition':'Wald 95% from the same REGENIE BETA and SE; P retained from REGENIE LOG10P, not recomputed from BETA/SE','threshold_counts':{str(t):int((valid & values.LOG10P.gt(-np.log10(t))).sum()) for t in [5e-8,3e-8,5e-9,1e-9]}}
 save(complete,receipt);return receipt

def main():
 state=N/'status/X_source_association_status.json'
 save(state,{'status':'waiting_for_site_FILTER_reconciliation_and_strict_Step1','pid':os.getpid()})
 sr=wait_for(N/'source_X_sites/RECONCILIATION_COMPLETE.json',N/'status/X_site_reconciliation_pid.json',[N/'status/X_site_reconciliation_failed.json'])
 sv=wait_for(N/'step1/VALIDATED.json',N/'status/autosomal_pid.json',[N/'status/autosomal_advance_failed.json'])
 assert sv['n']==386069 and sv['cases']==5705 and sv['loco_rows']==list(range(1,24))
 assert sha(N/'step1/model_1.loco')==sv['loco_sha256']
 assert sha(N/'source_X_sites/provider_X_site_QC.tsv.gz')==sr['output_sha256']
 pred=N/'step1/model_pred.list';assert pred.read_text().split()==['psoriasis_case',str(N/'step1/model_1.loco')]
 sites=pd.read_csv(N/'source_X_sites/provider_X_site_QC.tsv.gz',sep=chr(9)).set_index('source_index');assert sites.index.is_unique
 pilot=N/'X_projection_pilot/chunk_0000'
 if not (pilot/'HWE_COMPLETE.json').exists():hwe_pilot(pilot)
 pr=process(pilot,N/'results/X_source_pilot',sites,pred)
 assert pr['QC_passed']>0 and len(pr['models'])==2 and all(r['failed_or_unavailable']==0 for r in pr['models'].values())
 save(N/'status/X_source_association_pilot_validated.json',{'status':'validated','candidate_alleles':pr['candidate_alleles'],'QC_passed':pr['QC_passed'],'male_models':['02','01'],'scope':'source-projection pilot with distributed X positions, PAR/non-PAR and original leads; testing restricted to QC-passing alleles'})
 done={};root=N/'X_projected_source';outroot=N/'results/X_source'
 while True:
  for fail in [root/'FAILED.json',N/'status/X_HWE_FAILED.json']:
   if fail.exists():raise RuntimeError(str(fail)+' requires review')
  for folder in sorted(root.glob('chunk_*')):
   if folder.name in done or not (folder/'HWE_COMPLETE.json').exists():continue
   save(state,{'status':'running','completed_chunks':len(done),'current_chunk':folder.name,'pid':os.getpid()})
   done[folder.name]=process(folder,outroot,sites,pred)
   print('PSORIASIS_SOURCE_X_CHUNK_COMPLETE',folder.name,done[folder.name]['QC_passed'],flush=True)
  if (root/'ALL_HWE_COMPLETE.json').exists():
   expected=json.loads((root/'COMPLETE.json').read_text())['chunks']
   if len(done)==expected:break
  time.sleep(30)
 report={'status':'complete','completed_chunks':len(done),'candidate_alleles':sum(r['candidate_alleles'] for r in done.values()),'QC_passed':sum(r['QC_passed'] for r in done.values()),'failed_or_unavailable':{s:sum(r['models'].get(s,{}).get('failed_or_unavailable',0) for r in done.values()) for s in ['male02','male01']},'scope':'pooled X source-wide candidate scan; whole-output calibration and sex/ancestry stratification still required','pid':os.getpid()}
 save(state,report);save(outroot/'COMPLETE.json',report)
if __name__=='__main__':
 try:main()
 except Exception as e:
  save(N/'status/X_source_association_failed.json',{'error':str(e),'traceback':traceback.format_exc()});raise
