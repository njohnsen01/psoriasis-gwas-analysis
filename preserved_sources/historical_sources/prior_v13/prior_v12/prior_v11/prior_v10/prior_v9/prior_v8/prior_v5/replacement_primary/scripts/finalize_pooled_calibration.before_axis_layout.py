from pathlib import Path
import argparse,importlib.metadata as im,json,os,time,traceback
import numpy as np,pandas as pd
from scipy import stats,special
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reconcile_X_site_filters import sha,save
N=Path(__file__).resolve().parents[1]
def chi1(log10p):
 # Two-sided normal tail inversion in log space avoids P underflow.
 return special.ndtri_exp(-np.asarray(log10p,dtype=float)*np.log(10.)-np.log(2.))**2

def quantiles(log10p,max_points=60000):
 v=np.sort(np.asarray(log10p,dtype=float))[::-1];n=len(v);assert n>0 and np.isfinite(v).all() and (v>=0).all()
 middle=v[(n-1)//2:(n//2)+1];lam=float(chi1(middle).mean()/stats.chi2.ppf(.5,1))
 ranks=np.unique(np.concatenate([np.arange(1,min(n,1000)+1),np.linspace(1,n,min(n,max_points),dtype=int),np.arange(max(1,n-999),n+1)]))
 expected=-np.log10((ranks-.5)/n);observed=v[ranks-1]
 lower=-np.log10(stats.beta.ppf(.975,ranks,n-ranks+1));upper=-np.log10(stats.beta.ppf(.025,ranks,n-ranks+1))
 return lam,pd.DataFrame({'full_output_rank':ranks,'N_valid_tests':n,'expected_log10P':expected,'observed_log10P':observed,'null_95_low':lower,'null_95_high':upper})

def inputs(group):
 if group=='autosomes':
  status=json.loads((N/'status/autosomal_status.json').read_text());assert status['status']=='complete' and status['completed']==status['expected']==1104
  receipts=sorted((N/'results/autosomes').glob('*/COMPLETE.json'));assert len(receipts)==1104
  for path in receipts:
   r=json.loads(path.read_text());assert r['status']=='validated'
   if r['passed']:yield path.parent/'assoc_psoriasis_case.regenie',r['results_sha256'],r['passed'],r['estimated'],r['failed_or_unavailable']
 else:
  status=json.loads((N/'results/X_source/COMPLETE.json').read_text());assert status['status']=='complete'
  receipts=sorted((N/'results/X_source').glob('*/COMPLETE.json'));assert len(receipts)==status['completed_chunks']
  for path in receipts:
   r=json.loads(path.read_text());assert r['status']=='validated'
   if r['QC_passed']:
    model=r['models'][group];yield path.parent/(group+'_psoriasis_case.regenie'),model['results_sha256'],r['QC_passed'],model['estimated'],model['failed_or_unavailable']

def calculate(group,out):
 seen=set();parts=[];checks=[];bad=0;total=0;perchr=[]
 for p,expected_sha,n,e,u in inputs(group):
  assert sha(p)==expected_sha,'Result checksum differs: '+str(p)
  d=pd.read_csv(p,sep=chr(92)+'s+',dtype={'ID':str});assert len(d)==n and d.ID.is_unique and d.TEST.eq('ADD').all()
  overlap=seen.intersection(d.ID);assert not overlap,'Duplicate variants across result chunks require reconciliation'
  seen.update(d.ID)
  values=d[['BETA','SE','LOG10P']].apply(pd.to_numeric,errors='coerce');ok=np.isfinite(values).all(axis=1)&values.SE.gt(0)&values.LOG10P.ge(0)
  assert int(ok.sum())==e and int((~ok).sum())==u
  parts.append(values.loc[ok,'LOG10P'].to_numpy(float));bad+=u;total+=n
  checks.append({'file':str(p),'sha256':expected_sha,'tested':n,'finite_estimates':e,'unavailable':u})
  perchr.append({'chromosome':str(d.CHROM.iloc[0]),'tested':n,'finite_estimates':e,'unavailable':u})
 del seen
 allp=np.concatenate(parts);del parts
 lam,qq=quantiles(allp);qq.to_csv(out/(group+'_QQ_full_ranks.tsv.gz'),sep=chr(9),index=False)
 report={'group':group,'lambda_GC_from_full_output_P':lam,'N_tested_rows':total,'N_valid_tests':len(allp),'N_failed_or_unavailable':bad,'source_files':len(checks),'weights':'unit weight for every valid full-output test; no sampling or inclusion weights','lambda_method':'median of one-degree-of-freedom chi-square quantiles from the same REGENIE P values; log-space inverse normal avoids P underflow','QQ_method':'expected quantile = -log10((rank - 0.5)/N_valid_tests); observed values ranked using all valid tests; displayed points are a deterministic subset of these full ranks, never reranked','QQ_display_points':len(qq),'band':'pointwise 95% beta order-statistic null envelope, not a simultaneous confidence band','significance_counts':{str(t):int((allp>-np.log10(t)).sum()) for t in [5e-8,3e-8,5e-9,1e-9]},'scope':'pooled analysis; sex/ancestry-specific calibration remains separate'}
 save(out/(group+'_source_receipts.json'),checks);save(out/(group+'_calibration.json'),report)
 pd.DataFrame(perchr).groupby('chromosome',sort=False)[['tested','finite_estimates','unavailable']].sum().to_csv(out/(group+'_chromosome_counts.tsv'),sep=chr(9))
 fig,ax=plt.subplots(figsize=(5.2,5.2));ax.fill_between(qq.expected_log10P,qq.null_95_low,qq.null_95_high,color='#d4d8df',alpha=.65,linewidth=0)
 ax.scatter(qq.expected_log10P,qq.observed_log10P,s=3,color={'autosomes':'#334e68','male02':'#007c83','male01':'#94512a'}[group],rasterized=True,linewidths=0)
 limit=max(float(qq.expected_log10P.max()),float(qq.observed_log10P.max()))*1.025;ax.plot([0,limit],[0,limit],color='#8c939d',linewidth=.8);ax.set_xlim(0,limit);ax.set_ylim(0,limit)
 label='Autosomes' if group=='autosomes' else 'Chromosome X: male 0/'+group[-1]
 ax.set_title(label);ax.set_xlabel('Expected -log10(P)');ax.set_ylabel('Observed -log10(P)')
 ax.text(.04,.96,'N = '+format(len(allp),',')+chr(10)+'Lambda GC = '+format(lam,'.3f'),transform=ax.transAxes,va='top',fontsize=10)
 ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/(group+'_QQ.png'),dpi=300);fig.savefig(out/(group+'_QQ.pdf'));plt.close(fig)
 return report

def main(wait):
 while True:
  for p in [N/'status/autosomal_advance_failed.json',N/'status/X_source_association_failed.json']:
   if p.exists():raise RuntimeError(str(p)+' requires review')
  auto=N/'status/autosomal_status.json';x=N/'results/X_source/COMPLETE.json'
  a=json.loads(auto.read_text()) if auto.exists() else {}
  if a.get('status')=='complete' and a.get('completed')==a.get('expected')==1104 and x.exists():break
  if not wait:raise RuntimeError('Full pooled association output is not complete')
  time.sleep(30)
 O=N/'diagnostics/pooled_complete';O.mkdir(parents=True,exist_ok=True)
 assert not (O/'COMPLETE.json').exists(),'Preserve prior diagnostic output; inspect before rerunning'
 reports=[calculate(g,O) for g in ['autosomes','male02','male01']]
 save(O/'COMPLETE.json',{'status':'complete','groups':reports,'code_sha256':sha(Path(__file__)),'versions':{p:im.version(p) for p in ['numpy','pandas','scipy','matplotlib']},'scope':'full pooled autosomal and X diagnostics; final scientific review, sex/ancestry diagnostics and public-release checks are separate'})
 print('FULL_POOLED_CALIBRATION_COMPLETE',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--wait',action='store_true');a=p.parse_args()
 try:main(a.wait)
 except Exception as e:save(N/'status/pooled_calibration_failed.json',{'error':str(e),'traceback':traceback.format_exc()});raise
