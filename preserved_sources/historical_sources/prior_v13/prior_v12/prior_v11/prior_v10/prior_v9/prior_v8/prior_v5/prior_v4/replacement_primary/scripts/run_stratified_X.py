from pathlib import Path
import concurrent.futures,hashlib,importlib.metadata as im,json,os,threading,time,traceback
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reconcile_X_site_filters import sha,save
from run_source_X import wait_for,run
from finalize_pooled_calibration import quantiles
N=Path(__file__).resolve().parents[1];LABELS=['AMR','AFR','female','male','EUR'];STOP=threading.Event()
def process(folder,outroot,sites,pred,label):
 cfg=json.loads((N/'stratified_step1/strata.json').read_text())[label];ntotal=cfg['n'];ncase=cfg['cases'];ncontrol=cfg['controls'];inp=N/'stratified_step1'/label
 out=outroot/folder.name;out.mkdir(parents=True,exist_ok=True)
 complete=out/'COMPLETE.json'
 if complete.exists():
  receipt=json.loads(complete.read_text());assert receipt['status']=='validated'
  for s in receipt['models']:assert sha(out/(s+'_psoriasis_case.regenie'))==receipt['models'][s]['results_sha256']
  return receipt
 assert (folder/'HWE_COMPLETE.json').exists()
 pv=pd.read_csv(folder/'input.pvar',sep=chr(9)).set_index('ID');q=pd.read_csv(folder/'QC_with_HWE.private.tsv.gz',sep=chr(9)).set_index('ID');assert pv.index.equals(q.index) and q.index.is_unique
 subroot=N/('X_stratified_QC_pilot' if 'X_projection_pilot' in folder.parts else 'X_stratified_QC')
 sub=pd.read_csv(subroot/folder.name/(label+'.count_QC.private.tsv.gz'),sep=chr(9)).set_index('ID');assert sub.index.equals(q.index)
 q['pooled_PASS_except_site_FILTER']=q.PASS_except_site_FILTER
 for name in ['case_AC','case_AN','case_MAC','case_called','case_callrate','control_AC','control_AN','control_MAC','control_called','control_callrate','total_MAC','MAF','count_call_maf_PASS']:q[name]=sub[name]
 q['PASS_except_site_FILTER']=q.pooled_PASS_except_site_FILTER&q.count_call_maf_PASS
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
 ph=pd.read_csv(inp/'phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str})
 assert len(ps)==386069 and len(ph)==ntotal and set(ph.IID).issubset(set(ps.IID)) and ps.SEX.isin([1,2]).all() and int(ph.psoriasis_case.sum())==ncase
 receipt={'status':'validated','candidate_alleles':len(q),'QC_passed':len(selected),'models':{},'n_total_stratum':ntotal,'n_cases_total':ncase,'n_total_primary_cohort':386069,'stratum':label,'source':'unaffected provider PGEN with exact site FILTER and genotype count/call/MAF/female-control HWE QC','PGEN_sha256':sha(folder/'input.pgen'),'PVAR_sha256':sha(folder/'input.pvar'),'PSAM_sha256':sha(folder/'input.psam'),'site_reconciliation_sha256':sha(N/'source_X_sites/RECONCILIATION_COMPLETE.json'),'script_sha256':sha(Path(__file__))}
 for scale in (['02'] if label=='female' else ['02','01']):
  if not len(selected):continue
  prefix='male'+scale
  cmd=['regenie','--step','2','--bt','--pgen',str(folder/'input'),'--keep',str(inp/'keep.tsv'),'--extract',str(out/'passed_variants.txt'),'--phenoFile',str(inp/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(inp/'covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(pred),'--af-cc','--par-region','hg38','--threads','4','--out',str(out/prefix)]
  if scale=='01':cmd+=['--skip-dosage-comp']
  run(cmd,out/(prefix+'.console.log'))
  resultpath=out/(prefix+'_psoriasis_case.regenie');res=pd.read_csv(resultpath,sep=chr(92)+'s+')
  assert res.ID.is_unique and set(res.ID)==set(selected.index) and len(res)==len(selected) and res.TEST.eq('ADD').all()
  match=selected.reindex(res.ID)
  assert np.array_equal(res.ALLELE0.to_numpy(str),match.REF.to_numpy(str)) and np.array_equal(res.ALLELE1.to_numpy(str),match.ALT.to_numpy(str)) and np.array_equal(res.GENPOS.to_numpy(int),match.POS.to_numpy(int))
  log=(out/(prefix+'.log')).read_text()
  for term in ['number of individuals used in analysis = '+str(ntotal),str(ncase)+' cases and '+str(ncontrol)+' controls','LOCO predictions']:assert term in log,term
  assert res.N.between(1,ntotal).all()
  values=res[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');valid=np.isfinite(values).all(axis=1)&values.SE.gt(0)&values.LOG10P.ge(0)
  derived=res.copy();derived['model']='female012_male'+scale;derived['estimate_available']=valid
  derived['OR']=np.exp(values.BETA);derived['CI95_Wald_low']=np.exp(values.BETA-1.959963984540054*values.SE);derived['CI95_Wald_high']=np.exp(values.BETA+1.959963984540054*values.SE)
  derived['P_from_REGENIE_LOG10P']=np.power(10.,-values.LOG10P);derived['N_total_primary_cohort']=386069;derived['N_total_stratum']=ntotal;derived['stratum']=label
  qc=q.reindex(res.ID)
  for name in ['case_called','control_called','case_AC','case_AN','control_AC','control_AN','nonPAR']:
   derived[name]=qc[name].to_numpy()
  derived['N_nonmissing_genotype']=derived.case_called+derived.control_called
  derived.to_csv(out/(prefix+'.derived.private.tsv.gz'),sep=chr(9),index=False)
  receipt['models'][prefix]={'estimated':int(valid.sum()),'failed_or_unavailable':int((~valid).sum()),'results_sha256':sha(resultpath),'effect_allele':'ALLELE1 = selected provider ALT; other alleles collapsed to dosage zero','dosage_units':'female 0/1/2; non-PAR male 0/'+scale[-1]+'; PAR diploid','N_definition':'N is reported by REGENIE; total cohort and observed nonmissing genotype counts retained separately','CI_definition':'Wald 95% from the same REGENIE BETA and SE; P retained from REGENIE LOG10P, not recomputed from BETA/SE','threshold_counts':{str(t):int((valid & values.LOG10P.gt(-np.log10(t))).sum()) for t in [5e-8,3e-8,5e-9,1e-9]}}
 save(complete,receipt);return receipt
def calculate(group,out,stratum):
 seen=set();parts=[];checks=[];bad=0;total=0;perchr=[]
 for p,expected_sha,n,e,u in inputs(group,stratum):
  assert sha(p)==expected_sha,'Result checksum differs: '+str(p)
  d=pd.read_csv(p,sep=chr(92)+'s+',dtype={'ID':str});assert len(d)==n and d.ID.is_unique and d.TEST.eq('ADD').all()
  overlap=seen.intersection(d.ID);assert not overlap,'Duplicate variants across result chunks require reconciliation'
  seen.update(d.ID)
  values=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');ok=np.isfinite(values).all(axis=1)&values.SE.gt(0)&values.LOG10P.ge(0)
  assert int(ok.sum())==e and int((~ok).sum())==u
  parts.append(values.loc[ok,'LOG10P'].to_numpy(float));bad+=u;total+=n
  checks.append({'file':str(p),'sha256':expected_sha,'tested':n,'finite_estimates':e,'unavailable':u})
  for chromosome,indices in d.groupby('CHROM').groups.items():
   perchr.append({'chromosome':str(chromosome),'tested':len(indices),'finite_estimates':int(ok.loc[indices].sum()),'unavailable':int((~ok.loc[indices]).sum())})
 del seen
 allp=np.concatenate(parts);del parts
 lam,qq=quantiles(allp);qq.to_csv(out/(group+'_QQ_full_ranks.tsv.gz'),sep=chr(9),index=False)
 report={'group':group,'lambda_GC_from_full_output_P':lam,'N_tested_rows':total,'N_valid_tests':len(allp),'N_failed_or_unavailable':bad,'source_files':len(checks),'weights':'unit weight for every valid full-output test; no sampling or inclusion weights','lambda_method':'median of one-degree-of-freedom chi-square quantiles from the same REGENIE P values; log-space inverse normal avoids P underflow','QQ_method':'expected quantile = -log10((rank - 0.5)/N_valid_tests); observed values ranked using all valid tests; displayed points are a deterministic subset of these full ranks, never reranked','QQ_display_points':len(qq),'band':'pointwise 95% beta order-statistic null envelope, not a simultaneous confidence band','significance_counts':{str(t):int((allp>-np.log10(t)).sum()) for t in [5e-8,3e-8,5e-9,1e-9]},'stratum':stratum,'scope':'Full chromosome X within-cohort calibration under pooled source/HWE QC and subgroup count/call/MAF QC; not independent replication'}
 save(out/(group+'_source_receipts.json'),checks);save(out/(group+'_calibration.json'),report)
 pd.DataFrame(perchr).groupby('chromosome',sort=False)[['tested','finite_estimates','unavailable']].sum().to_csv(out/(group+'_chromosome_counts.tsv'),sep=chr(9))
 fig,ax=plt.subplots(figsize=(5.2,5.2));ax.fill_between(qq.expected_log10P,qq.null_95_low,qq.null_95_high,color='#d4d8df',alpha=.65,linewidth=0)
 ax.scatter(qq.expected_log10P,qq.observed_log10P,s=3,color={'autosomes':'#334e68','male02':'#007c83','male01':'#94512a'}[group],rasterized=True,linewidths=0)
 xmax=float(qq.expected_log10P.max())*1.025;ymax=max(xmax,float(qq.observed_log10P.max())*1.025);ax.plot([0,xmax],[0,xmax],color='#8c939d',linewidth=.8);ax.set_xlim(0,xmax);ax.set_ylim(0,ymax)
 label='Autosomes' if group=='autosomes' else 'Chromosome X: male 0/'+group[-1]
 ax.set_title('Chromosome X: '+stratum+('' if stratum=='female' else ' (male 0/'+group[-1]+')'));ax.set_xlabel('Expected -log10(P)');ax.set_ylabel('Observed -log10(P)')
 ax.text(.04,.96,'N = '+format(len(allp),',')+chr(10)+'Lambda GC = '+format(lam,'.3f'),transform=ax.transAxes,va='top',fontsize=10)
 ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/(group+'_QQ.png'),dpi=300);fig.savefig(out/(group+'_QQ.pdf'));plt.close(fig)
 return report

def validate_prediction(label):
 inp=N/'stratified_step1'/label
 meta=wait_for(inp/'VALIDATED.json',N/'status/stratified_step1_pid.json',[N/'stratified_step1/FAILED.json'])
 assert meta['status']=='validated' and meta['stratum']==label and meta['loco_rows']==list(range(1,24))
 assert sha(inp/'model_1.loco')==meta['loco_sha256']
 for fname,key in [('phenotype.tsv','phenotype_sha256'),('covariates.tsv','covariate_sha256')]:assert sha(inp/fname)==meta[key]
 pred=inp/'model_pred.list';assert pred.read_text().split()==['psoriasis_case',str(inp/'model_1.loco')]
 return pred,meta

def one_stratum(label,sites):
 try:
  out=N/'results/X_stratified'/label;out.mkdir(parents=True,exist_ok=True);pred,meta=validate_prediction(label)
  pilot=process(N/'X_projection_pilot/chunk_0000',N/'results/X_stratified_pilot'/label,sites,pred,label)
  assert all(x['failed_or_unavailable']==0 for x in pilot['models'].values()),'Subgroup pilot model failure requires review: '+label
  save(out/'PILOT_VALIDATED.json',{'status':'validated','QC_passed':pilot['QC_passed'],'models':list(pilot['models']),'scope':'QC-passing pilot only; zero eligible variants recorded without a significance claim'})
  root=N/'X_projected_source';expected=json.loads((root/'COMPLETE.json').read_text())['chunks'];folders=sorted(root.glob('chunk_*'));assert len(folders)==expected;results=[]
  for folder in folders:
   if STOP.is_set():raise RuntimeError('Another subgroup failed; stop before the next chunk')
   assert (folder/'HWE_COMPLETE.json').exists() and (N/'X_stratified_QC'/folder.name/'COMPLETE.json').exists()
   save(out/'status.json',{'status':'running','stratum':label,'completed_chunks':len(results),'expected_chunks':expected,'current':folder.name,'pid':os.getpid()})
   results.append(process(folder,out,sites,pred,label));print('STRATIFIED_X_ASSOCIATION_COMPLETE',label,folder.name,flush=True)
  record={'status':'complete','stratum':label,'chunks':len(results),'QC_passed':sum(x['QC_passed'] for x in results),'candidate_alleles':sum(x['candidate_alleles'] for x in results),'n_total_stratum':meta['n'],'n_cases_stratum':meta['cases'],'covariates':meta['covariates'],'LOCO_row':23,'LOCO_sha256':meta['loco_sha256'],'scales':['02'] if label=='female' else ['02','01'],'interpretation':'Within-cohort exploratory association and calibration; no independent replication. Subgroup eligibility does not establish power.'}
  save(out/'COMPLETE.json',record);return record
 except Exception:
  STOP.set();raise

def inputs(group,stratum):
 root=N/'results/X_stratified'/stratum;r=json.loads((root/'COMPLETE.json').read_text());assert r['status']=='complete'
 paths=sorted(root.glob('chunk_*/COMPLETE.json'));assert len(paths)==r['chunks']
 for path in paths:
  meta=json.loads(path.read_text());assert meta['status']=='validated'
  if meta['QC_passed']:
   model=meta['models'][group];yield path.parent/(group+'_psoriasis_case.regenie'),model['results_sha256'],meta['QC_passed'],model['estimated'],model['failed_or_unavailable']

def main():
 save(N/'status/X_stratified_association_status.json',{'status':'waiting_for_source_QC_and_predictions','pid':os.getpid()})
 proof=json.loads((N/'X_stratified_QC_pilot/INDEPENDENT_PLINK_VALIDATION.json').read_text());assert proof['status']=='passed' and proof['comparisons']==830
 sr=wait_for(N/'source_X_sites/RECONCILIATION_COMPLETE.json',N/'status/X_site_reconciliation_pid.json',[N/'status/X_site_reconciliation_failed.json'])
 wait_for(N/'status/X_source_association_pilot_validated.json',N/'status/X_source_association_pid.json',[N/'status/X_source_association_failed.json'])
 wait_for(N/'X_projected_source/ALL_HWE_COMPLETE.json',N/'status/X_HWE_pid.json',[N/'status/X_HWE_FAILED.json',N/'X_projected_source/FAILED.json'])
 cq=wait_for(N/'X_stratified_QC/COMPLETE.json',N/'status/X_stratified_QC_pid.json',[N/'status/X_stratified_QC_failed.json']);assert cq['status']=='complete'
 assert sha(N/'source_X_sites/provider_X_site_QC.tsv.gz')==sr['output_sha256'];sites=pd.read_csv(N/'source_X_sites/provider_X_site_QC.tsv.gz',sep=chr(9)).set_index('source_index');assert sites.index.is_unique;records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for future in concurrent.futures.as_completed([pool.submit(one_stratum,label,sites) for label in LABELS]):records.append(future.result())
 for r in records:
  out=N/'diagnostics/X_stratified_complete'/r['stratum'];out.mkdir(parents=True,exist_ok=True)
  if not r['QC_passed']:save(out/'COMPLETE.json',{'status':'no_eligible_tests','stratum':r['stratum']});continue
  reports=[calculate('male'+scale,out,r['stratum']) for scale in r['scales']]
  save(out/'COMPLETE.json',{'status':'complete','stratum':r['stratum'],'diagnostics':reports,'script_sha256':sha(Path(__file__))})
 save(N/'results/X_stratified/COMPLETE.json',{'status':'complete','strata':records,'diagnostics':'diagnostics/X_stratified_complete','script_sha256':sha(Path(__file__))})
 save(N/'status/X_stratified_association_status.json',{'status':'complete','strata':LABELS,'pid':os.getpid()});print('FULL_STRATIFIED_X_ASSOCIATION_AND_CALIBRATION_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(N/'status/X_stratified_association_failed.json',{'error':str(e),'traceback':traceback.format_exc()});raise
