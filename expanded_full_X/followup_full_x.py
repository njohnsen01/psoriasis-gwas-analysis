"""Independent ML/LRT and Firth checks of predeclared leading scan signals."""
import os
from pathlib import Path
import pandas as pd,numpy as np,json,sys,hashlib,subprocess,datetime
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');O=R/'full_x_followup';Q=R.parent/'chen_scan_benchmark';B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/X_projected_source')
def main():
 assert (R/'full_x/SCAN_COMPLETE.json').exists();O.mkdir(exist_ok=True)
 x=pd.read_csv(R/'full_x/full_X_joint_Wald.tsv.gz',sep='\t');x=x[(x.ERRCODE=='.')&np.isfinite(x.P_CHISQ_3DF)&np.isfinite(x.WALD_CHISQ_3DF)]
 selected=[]
 for row in x.sort_values('P_CHISQ_3DF').itertuples():
  if all(abs(row.POS-p)>1000000 for p in selected):selected.append(row.POS)
  if len(selected)==10:break
 leads=x[x.POS.isin(selected)].sort_values('P_CHISQ_3DF').drop_duplicates('POS')
 ids=list(leads.ID)
 ids+=['chrX:5995921:G:A','chrX:5996494:GCT:G','chrX:5998724:A:G'];ids=list(dict.fromkeys(ids))
 sources={}
 for pv in sorted(B.glob('chunk_*/input.pvar')):
  v=pd.read_csv(pv,sep='\t',comment=None,skiprows=sum(1 for l in pv.open() if l.startswith('##')))
  for i,row in v[v.ID.isin(ids)].iterrows():sources[row.ID]=(pv.parent,int(i))
 assert set(sources)==set(ids)
 sys.path.insert(0,f'{os.environ["PGENLIB_SITE_DIR"]}');import pgenlib
 cov=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});ph=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str});d=cov.merge(ph,on=['FID','IID'],validate='one_to_one');assert len(d)==357530
 audit=[]
 for j,vid in enumerate(ids):
  folder,index=sources[vid];s=pd.read_csv(folder/'input.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str}).rename(columns={'#FID':'FID'});s['sample_index']=np.arange(len(s));aligned=d[['FID','IID','S']].merge(s[['FID','IID','SEX','sample_index']],on=['FID','IID'],validate='one_to_one');assert len(aligned)==len(d) and np.array_equal(aligned.S,(aligned.SEX==1).astype(int))
  arr=np.empty(len(s),np.int8)
  with pgenlib.PgenReader(str(folder/'input.pgen').encode(),raw_sample_ct=len(s)) as reader:reader.read(index,arr)
  g=arr[aligned.sample_index.to_numpy()].astype(float);g[g<0]=np.nan;assert set(g[np.isfinite(g)&(d.S==1)])<=set([0.,2.]);d[f'g{j}']=g
  z=x[x.ID==vid];audit.append({'key':f'g{j}','ID':vid,'n_complete':int(np.isfinite(g).sum()),'scan_p':float(z.P_CHISQ_3DF.iloc[0]) if len(z) else None})
 d[['CASE','S']+[f'Q{i}' for i in range(1,25)]+[f'g{i}' for i in range(len(ids))]].to_csv(O/'selected.private.csv',index=False)
 (O/'selected_variants.json').write_text(json.dumps(audit,indent=2)+'\n')
 # Reuse the already installed research R library; locate the package directory.
 lib=[]
 for candidate in [R.parent/'r_library',R.parent/'Rlib',R.parent/'R_library']:
  if (candidate/'logistf/DESCRIPTION').exists():lib.append(candidate)
 if not lib:
  for candidate in R.parent.glob('**/logistf/DESCRIPTION'):lib.append(candidate.parent.parent)
 assert lib,'Previously installed logistf library not found'
 cmd=['/opt/conda/envs/jupyter/bin/Rscript',str(R/'followup_full_x.R'),str(O),str(lib[0]),str(len(ids))]
 with (O/'console.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,'Follow-up R script failed'
 results=pd.read_csv(O/'followup_summary.csv').merge(pd.DataFrame(audit),on='key',validate='many_to_one');results.to_csv(O/'followup_named.csv',index=False)
 (O/'FOLLOWUP_COMPLETE.json').write_text(json.dumps({'status':'independent_lead_checks_completed_pending_review','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selected':len(ids),'failed':int((results.status!='estimated').sum()),'limitation':'Selection on the scan is follow-up, not independent replication; 1 Mb spacing does not establish LD independence'},indent=2)+'\n')
if __name__=='__main__':main()
