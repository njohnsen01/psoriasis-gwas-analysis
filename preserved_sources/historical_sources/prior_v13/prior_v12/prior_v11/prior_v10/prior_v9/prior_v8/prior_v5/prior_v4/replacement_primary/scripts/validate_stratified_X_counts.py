from pathlib import Path
import json,os,subprocess,traceback
import numpy as np,pandas as pd
from reconcile_X_site_filters import save,sha
N=Path(__file__).resolve().parents[1];root=N/'X_stratified_QC_pilot';folder=N/'X_projection_pilot/chunk_0000';out=root/'independent_PLINK_check';out.mkdir(exist_ok=True)
try:
 checks=[]
 for label in ['EUR','AFR','AMR','female','male']:
  ph=pd.read_csv(N/'stratified_step1'/label/'phenotype.tsv',sep=chr(9),dtype={'FID':str,'IID':str})
  qc=pd.read_csv(root/'chunk_0000'/(label+'.count_QC.private.tsv.gz'),sep=chr(9)).set_index('ID')
  for name,value in [('case',1),('control',0)]:
   keep=out/(label+'.'+name+'.keep.private.tsv');ph.loc[ph.psoriasis_case.eq(value),['FID','IID']].to_csv(keep,sep=chr(9),index=False,header=False)
   prefix=out/(label+'.'+name)
   cmd=['plink2','--pfile',str(folder/'input'),'--split-par','hg38','--keep',str(keep),'--freq','counts','--missing','variant-only','--threads','2','--memory','4096','--out',str(prefix)]
   save(prefix.with_suffix(prefix.suffix+'.command.json'),cmd)
   with prefix.with_suffix(prefix.suffix+'.console.log').open('w') as handle:r=subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT)
   assert r.returncode==0,str(prefix)
   ac=pd.read_csv(str(prefix)+'.acount',sep=chr(92)+'s+').set_index('ID');miss=pd.read_csv(str(prefix)+'.vmiss',sep=chr(92)+'s+').set_index('ID')
   assert ac.index.is_unique and set(ac.index)==set(qc.index) and set(miss.index)==set(qc.index)
   ac=ac.reindex(qc.index);miss=miss.reindex(qc.index)
   assert np.array_equal(ac.ALT_CTS.to_numpy(float),qc[name+'_AC'].to_numpy(float)),(label,name,'AC')
   assert np.array_equal(ac.OBS_CT.to_numpy(float),qc[name+'_AN'].to_numpy(float)),(label,name,'AN')
   assert np.array_equal((miss.OBS_CT-miss.MISSING_CT).to_numpy(int),qc[name+'_called'].to_numpy(int)),(label,name,'called')
   assert np.allclose(1-miss.F_MISS.to_numpy(float),qc[name+'_callrate'].to_numpy(float),rtol=0,atol=5e-7),(label,name,'callrate precision')
   checks.append({'stratum':label,'case_control':name,'variants':len(qc),'AC_AN_called_exact_match':True,'callrate_atol_5e-7':True})
 save(root/'INDEPENDENT_PLINK_VALIDATION.json',{'status':'passed','checks':checks,'comparisons':len(checks)*83,'scope':'All 83 pilot allele projections, five subgroups and case/control separately; biological chromosome counts and called samples independently verified with PLINK2','script_sha256':sha(Path(__file__))})
 print('INDEPENDENT_STRATIFIED_COUNTS_VALIDATED',len(checks)*83,flush=True)
except Exception as e:save(root/'INDEPENDENT_PLINK_FAILED.json',{'error':str(e),'traceback':traceback.format_exc()});raise
