from pathlib import Path
import json,sys,subprocess,hashlib,traceback,re
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Y=Path(__file__).resolve().parents[1];N=Y.parent;I=N/'stratified_step1/male';G=Y/'projected';O=Y/'results';O.mkdir(exist_ok=True)
sys.path.insert(0,str(N/'scripts'))
from finalize_pooled_calibration import quantiles
from reconcile_X_site_filters import sha,save
def main():
 meta=json.loads((G/'COMPLETE.json').read_text());assert meta['status']=='validated' and not meta['pilot']
 for f,h in meta['inputs'].items():assert sha(G/f)==h
 validation=json.loads((I/'VALIDATED.json').read_text());assert validation['status']=='validated' and validation['n']==151282 and validation['cases']==2259
 for f,key in [('phenotype.tsv','phenotype_sha256'),('covariates.tsv','covariate_sha256'),('model_1.loco','loco_sha256')]:assert sha(I/f)==validation[key]
 rows=[]
 with (I/'model_1.loco').open() as f:
  next(f)
  for line in f:rows.append(int(line.split(maxsplit=1)[0]))
 assert rows==list(range(1,24))
 predline=(I/'model_pred.list').read_text().strip().split();assert len(predline)==2 and Path(predline[1]).resolve()==(I/'model_1.loco').resolve()
 cov=pd.read_csv(I/'covariates.tsv',sep=chr(9),nrows=0).columns.tolist();assert cov==['FID','IID','age','age2']+['PC'+str(i) for i in range(1,11)]
 q=pd.read_csv(G/'final_variant_QC.private.tsv.gz',sep=chr(9)).set_index('ID');pv=pd.read_csv(G/'input.pvar',sep=chr(9)).set_index('ID');assert q.index.equals(pv.index)
 selected=q.index[q.PASS];assert len(selected)>0
 reports={};models={}
 for scale in ['02','01']:
  name='male'+scale;prefix=O/name
  cmd=['regenie','--step','2','--bt','--pgen',str(G/'input'),'--extract',str(G/'passed_variants.txt'),'--keep',str(I/'keep.tsv'),'--phenoFile',str(I/'phenotype.tsv'),'--phenoCol','psoriasis_case','--covarFile',str(I/'covariates.tsv'),'--bsize','200','--firth','--approx','--pThresh','.01','--pred',str(I/'model_pred.list'),'--sex-specific','male','--af-cc','--threads','4','--out',str(prefix)]
  if scale=='01':cmd+=['--skip-dosage-comp']
  save(O/(name+'.command.json'),cmd)
  with (O/(name+'.console.log')).open('w') as f:run=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
  assert run.returncode==0,(name,run.returncode)
  log=Path(str(prefix)+'.log').read_text()
  for text in ['End time:','number of individuals used in analysis = 151282','2259 cases and 149023 controls','Chromosome 23','LOCO predictions']:assert text in log,(name,text)
  fail=re.findall(r'Number of failed tests\s*:\s*\((\d+)/(\d+)\)',log);assert fail and all(int(a)==0 for a,b in fail)
  assert re.search(r'Number of ignored tests due to low MAC\s*:\s*0',log)
  rp=O/(name+'_psoriasis_case.regenie');d=pd.read_csv(rp,sep=r'\s+').set_index('ID');assert d.index.is_unique and set(d.index)==set(selected) and d.TEST.eq('ADD').all()
  d=d.reindex(selected);v=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');assert np.isfinite(v).all().all() and v.SE.gt(0).all() and v.LOG10P.ge(0).all()
  source=pv.reindex(d.index);assert d.ALLELE0.eq(source.REF).all() and d.ALLELE1.eq(source.ALT).all() and d.GENPOS.eq(source.POS).all()
  assert d.CHROM.eq(23).all();d['reported_CHROM']='Y';d['engine_CHROM']=23
  d['P_from_REGENIE_LOG10P']=10.**(-v.LOG10P);d['OR_per_coded_unit']=np.exp(v.BETA);d['CI95_Wald_low']=np.exp(v.BETA-1.959963984540054*v.SE);d['CI95_Wald_high']=np.exp(v.BETA+1.959963984540054*v.SE)
  d['carrier_vs_noncarrier_logOR']=v.BETA*(2 if scale=='02' else 1);d['carrier_vs_noncarrier_OR']=np.exp(d.carrier_vs_noncarrier_logOR)
  d['N_total_male_cohort']=151282;d['case_total_male_cohort']=2259;d['control_total_male_cohort']=149023
  for c in ['case_AC','case_AN','control_AC','control_AN','case_called','control_called','case_callrate','control_callrate','MAF']:d[c]=q.loc[d.index,c]
  d['N_nonmissing_genotypes']=d.case_called+d.control_called
  d['Y_family_BH_q']=stats.false_discovery_control(d.P_from_REGENIE_LOG10P.to_numpy(),method='bh')
  d.to_csv(O/(name+'.derived.private.tsv.gz'),sep=chr(9))
  lam,qq=quantiles(v.LOG10P.to_numpy());qq.to_csv(O/(name+'_QQ_full_ranks.private.tsv.gz'),sep=chr(9),index=False)
  report={'model':name,'valid_tests':len(d),'failed_or_unavailable':0,'lambda_GC':lam,'weights':'unit weights for lambda and complete-data QQ ranks; no significance-enriched reranking','threshold_counts':{str(t):int(v.LOG10P.gt(-np.log10(t)).sum()) for t in [5e-8,3e-8,5e-9,1e-9]},'Y_family_Bonferroni_threshold':.05/len(d),'Y_family_Bonferroni_count':int(d.P_from_REGENIE_LOG10P.lt(.05/len(d)).sum()),'Y_family_BH_q_lt_05':int(d.Y_family_BH_q.lt(.05).sum()),'dosage_units':'male 0/'+scale[-1]+'; selected ALT vs other alleles','prediction_row':23,'source_output_sha256':sha(rp),'model_log_sha256':sha(Path(str(prefix)+'.log')),'cohort_denominator':'eligible total distinct from variant-specific REGENIE N and nonmissing genotype count','release':'HELD'}
  save(O/(name+'_calibration.json'),report);reports[name]=report;models[name]=d
  fig,ax=plt.subplots(1,2,figsize=(7.2,3.6));ax[0].scatter(d.GENPOS/1e6,v.LOG10P,s=3,color='#247b87',linewidths=0,rasterized=True)
  ax[0].axhline(-np.log10(5e-8),ls='--',color='#a95142',lw=.8);ax[0].axhline(-np.log10(5e-9),ls=':',color='#737373',lw=.8)
  ax[0].set(xlabel='Chromosome Y position (Mb; GRCh38)',ylabel='-log10(P)',title='A  Non-PAR Y: male 0/'+scale[-1],ylim=(0,max(10,float(v.LOG10P.max())*1.04)))
  ax[1].fill_between(qq.expected_log10P,qq.null_95_low,qq.null_95_high,color='#dce3e6');ax[1].scatter(qq.expected_log10P,qq.observed_log10P,s=3,color='#247b87',linewidths=0,rasterized=True)
  xmax=float(qq.expected_log10P.max())*1.05;ymax=max(xmax,float(qq.observed_log10P.max())*1.05);ax[1].plot([0,xmax],[0,xmax],color='#87959a',lw=.7)
  ax[1].set(xlim=(0,xmax),ylim=(0,ymax),xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title='B  Full-output calibration')
  ax[1].text(.03,.97,'Tests = '+format(len(d),',')+chr(10)+'Lambda GC = '+format(lam,'.3f'),transform=ax[1].transAxes,va='top',fontsize=8)
  for a in ax:a.spines[['top','right']].set_visible(False)
  fig.subplots_adjust(left=.1,right=.98,bottom=.19,top=.90,wspace=.4)
  for ext in ['pdf','svg']:fig.savefig(O/('Chromosome_Y_'+name+'.'+ext))
  fig.savefig(O/('Chromosome_Y_'+name+'_600dpi.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'});fig.savefig(O/('Chromosome_Y_'+name+'_preview.png'),dpi=150);plt.close(fig)
  print('Y_MODEL_COMPLETE',name,len(d),flush=True)
 a=models['male02'];b=models['male01'];assert a.index.equals(b.index) and a.N.equals(b.N)
 checks={'max_abs_beta01_minus_twice_beta02':float(np.max(np.abs(b.BETA-2*a.BETA))),'max_abs_se01_minus_twice_se02':float(np.max(np.abs(b.SE-2*a.SE))),'max_abs_LOG10P_difference':float(np.max(np.abs(a.LOG10P-b.LOG10P)))}
 assert np.allclose(b.BETA,2*a.BETA,rtol=5e-5,atol=1e-6) and np.allclose(b.SE,2*a.SE,rtol=5e-5,atol=1e-6) and np.allclose(a.LOG10P,b.LOG10P,rtol=5e-5,atol=1e-6),checks
 save(O/'DOSAGE_COMPARISON.json',dict(status='validated',**checks,interpretation='male-only scale reparameterization; not independent replication'))
 save(O/'COMPLETE.json',{'status':'validated','reports':reports,'dosage_comparison':checks,'QC_sha256':sha(G/'COMPLETE.json'),'prediction_validation_sha256':sha(I/'VALIDATED.json'),'script_sha256':sha(Path(__file__)),'scope':'Source-corrected non-PAR Y follow-up; figure QA and manuscript/supplement integration pending','release':'HELD'})
 print('Y_SOURCE_CORRECTED_MODELS_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:save(O/'FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
