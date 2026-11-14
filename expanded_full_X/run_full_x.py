"""Resumable full nonPAR X joint-Wald scan; starts only after independent validation."""
from pathlib import Path
import json,time,subprocess,concurrent.futures,datetime,hashlib,os
import pandas as pd,numpy as np
from scipy.stats import chi2
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x';Q=R.parent/'chen_scan_benchmark'
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911');P='/opt/workbench-tools/binaries/bin/plink2'
def one(item):
 chunk,ids=item;D=O/chunk;D.mkdir(exist_ok=True)
 if (D/'DONE.json').exists():return json.loads((D/'DONE.json').read_text())
 src=B/'X_projected_source'/chunk/'input'
 for ext in ['.pgen','.psam']:
  dst=D/('proxy'+ext)
  if not dst.exists():dst.symlink_to(str(src)+ext)
 with Path(str(src)+'.pvar').open() as f,(D/'proxy.pvar').open('w') as w:
  for l in f:
   if l.startswith('##'):continue
   if l.startswith('#'):w.write(l);continue
   a=l.rstrip().split('\t');assert a[0] in ['X','chrX','23'];a[0]='1';w.write('\t'.join(a)+'\n')
 (D/'extract.txt').write_text('\n'.join(ids)+'\n')
 cmd=[P,'--pfile',str(D/'proxy'),'--keep',str(Q/'keep.private.tsv'),'--covar',str(Q/'covariates.private.tsv'),'--pheno',str(Q/'phenotype.private.tsv'),'--extract',str(D/'extract.txt'),'--glm','genotypic','interaction','no-x-sex','no-firth','hide-covar','--parameters','1-28','--tests','1,2,28','--threads','4','--memory','10000','--out',str(D/'fit')]
 (D/'COMMAND.json').write_text(json.dumps(cmd));t=time.time()
 with (D/'console.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(chunk,p.returncode)
 df=pd.read_csv(D/'fit.CASE.glm.logistic',sep=r'\s+');joint=df[df.TEST=='USER_3DF'].copy();assert len(joint)==len(ids) and set(joint.ID)==set(ids)
 joint['#CHROM']='X';joint['WALD_CHISQ_3DF']=3*joint.Z_OR_F_STAT;joint['P_CHISQ_3DF']=chi2.sf(joint.WALD_CHISQ_3DF,3)
 joint.to_csv(D/'joint.tsv',sep='\t',index=False)
 rec={'chunk':chunk,'variants':len(joint),'valid_p':int(joint.P.notna().sum()),'errors':joint.ERRCODE.fillna('OK').value_counts().to_dict(),'seconds':time.time()-t}
 (D/'DONE.json').write_text(json.dumps(rec,indent=2)+'\n');return rec

def main():
 assert json.loads((R/'proxy_benchmark/VALIDATED.json').read_text())['status']=='passed','Independent implementation validation required'
 O.mkdir(exist_ok=True);jobs=[];total=0
 for p in sorted((B/'X_projected_source').glob('chunk_*/input.pvar')):
  passed=set((B/'results/X_source'/p.parent.name/'passed_variants.txt').read_text().split())
  ids=[]
  for l in p.open():
   if l.startswith('#'):continue
   a=l.rstrip().split('\t')
   if a[2] in passed and 2781479<int(a[1])<155701383:ids.append(a[2])
  if ids:jobs.append((p.parent.name,ids));total+=len(ids)
 assert total==334926,total
 protocol={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Full previously QC-passed common nonPAR X; not PAR or rare-variant scan','variants':total,'chunks':len(jobs),'n_before_variant_missingness':357530,'cases_before_variant_missingness':5324,'method':'PLINK 3df joint Wald ADD + female DOMDEV + ADDxsex, expanded nuisance model; temporary chr1 wrapper, X restored in outputs','limitation':'Ordinary logistic; provider-pruned cohort but no mixed model/LOCO; sparse-cell/convergence failures reported. Implementation validation is not tail calibration.','primary_reference':'df3 chi-square; also retain native PLINK reported P','followup':'Top/threshold-crossing signals require exact/Firth verification','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 if not (O/'PROTOCOL.json').exists():(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2)+'\n')
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
  for rec in ex.map(one,jobs):print(json.dumps(rec),flush=True)
 allj=pd.concat([pd.read_csv(O/c/'joint.tsv',sep='\t') for c,_ in jobs],ignore_index=True);assert len(allj)==total and not allj.ID.duplicated().any()
 allj.to_csv(O/'full_X_joint_Wald.tsv.gz',sep='\t',index=False,compression='gzip')
 valid=allj.P_CHISQ_3DF.notna() & np.isfinite(allj.WALD_CHISQ_3DF) & (allj.WALD_CHISQ_3DF>=0)
 z=allj[valid].copy();z.sort_values('P_CHISQ_3DF').head(100).to_csv(O/'top100.tsv',sep='\t',index=False)
 n=len(z);q=z.P_CHISQ_3DF.to_numpy();qq=pd.DataFrame({'expected':-np.log10((np.arange(1,n+1)-.5)/n),'observed':-np.log10(np.maximum(np.sort(q),np.nextafter(0.,1.)))})
 qq.iloc[np.unique(np.r_[np.arange(min(1000,n)),np.linspace(0,n-1,min(n,5000)).astype(int)])].to_csv(O/'qq_points.tsv',sep='\t',index=False)
 bins=[]
 for lo,hi in [(0,.01),(.01,.05),(.05,.5)]:
  a=np.minimum(z.A1_FREQ,1-z.A1_FREQ);sel=(a>lo)&(a<=hi);s=z[sel]
  bins.append({'maf_interval':[lo,hi],'n':len(s),'lambda_df3':float(s.WALD_CHISQ_3DF.median()/chi2.ppf(.5,3)) if len(s) else None})
 rec={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'scan_completed_pending_scientific_review','attempted':total,'estimated':n,'unestimated':total-n,'lambda_df3':float(z.WALD_CHISQ_3DF.median()/chi2.ppf(.5,3)),'minimum_p':float(q.min()),'below_5e8':int((q<5e-8).sum()),'errors':allj.ERRCODE.fillna('OK').value_counts().to_dict(),'maf_diagnostics':bins}
 (O/'SCAN_COMPLETE.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
if __name__=='__main__':main()
