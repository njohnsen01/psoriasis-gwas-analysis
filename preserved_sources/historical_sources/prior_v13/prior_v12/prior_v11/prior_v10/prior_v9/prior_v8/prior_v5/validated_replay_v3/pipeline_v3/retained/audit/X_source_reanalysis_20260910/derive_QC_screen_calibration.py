from pathlib import Path
import json,time,hashlib
import numpy as np,pandas as pd
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;Q=R/'private_panel_QC';D=R/'corrected_cohort_diagnostics';O=R/'QC_screen_calibration';O.mkdir(exist_ok=True)
results=[]
for an in ['EUR','AFR','AMR','pooled']:
 q=pd.read_csv(Q/(an+'.variant_QC_mask.private.tsv'),sep='\t',index_col='ID')
 for model in ['male02','male01']:
  receipt=D/(an+'_'+model+'_panel.validated.json')
  while not receipt.exists():
   if (D/'pipeline_error.json').exists():raise RuntimeError((D/'pipeline_error.json').read_text())
   time.sleep(10)
  src=D/(an+'_'+model+'_panel_psoriasis_case.regenie');d=pd.read_csv(src,sep=r'\s+').set_index('ID');assert d.index.is_unique and len(d)==376171
  valid=np.isfinite(d.LOG10P)&(d.LOG10P>=0)
  for mask in ['mac_call_only_pass','count_call_pass','plus_two_sided_HWE_pass']:
   sel=q[mask].reindex(d.index);assert sel.notna().all();v=d[valid&sel]
   rec={'stratum':an,'model':model,'screen':mask,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'fixed original X test panel; same association outputs, variant screen only; no final genotype remasking/sample-QC refit'}
   for reg,t in {'X_all':v,'X_nonPAR':v[(v.GENPOS>2781479)&(v.GENPOS<155701383)],'X_PAR':v[(v.GENPOS<=2781479)|(v.GENPOS>=155701383)]}.items():
    lp=t.LOG10P.to_numpy();p=np.maximum(10.**(-lp),np.finfo(float).tiny);lam=float(np.median(chi2.isf(p,1))/chi2.ppf(.5,1));rec[reg]={'n':len(p),'lambda_gc':lam}
    if reg=='X_nonPAR':
     obs=np.sort(lp);expected=-np.log10((np.arange(len(lp),0,-1)-.5)/len(lp));fig,ax=plt.subplots(figsize=(5,5));ax.scatter(expected,obs,s=1,rasterized=True);ax.plot([0,expected.max()],[0,expected.max()],color='gray',lw=.8);ax.set(xlabel='Expected -log10(P)',ylabel='Observed -log10(P)',title=f'{an} {model}: {mask}\nnon-PAR X; n={len(lp)}; lambda={lam:.4f}');fig.tight_layout();fig.savefig(O/(an+'_'+model+'_'+mask+'_nonPAR_QQ.pdf'));plt.close(fig)
   results.append(rec)
   (O/'progress.json').write_text(json.dumps(results,indent=2))
  print('COMPLETE',an,model,flush=True)
(O/'complete.json').write_text(json.dumps({'status':'complete','qq_and_lambda_use':'all finite output P values passing the specified variant screen; unit weights; no significant-variant enrichment','results':results},indent=2))
