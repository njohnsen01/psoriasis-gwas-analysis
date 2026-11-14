from pathlib import Path
import time,json,pandas as pd,numpy as np
from scipy.stats import chi2
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/joint_G_sex_20260914/full_X')
while not (O/'DRIVER_FINISHED.json').exists():time.sleep(30)
frames=[];audit=[]
for p in sorted(O.glob('chunk_*')):
 f=p/'joint_psoriasis_case.regenie'
 if not f.exists():audit.append({'chunk':p.name,'missing_output':True});continue
 d=pd.read_csv(f,sep=r'\s+');j=d[d.TEST.eq('ADD-INT_2DF')].copy();j['chunk']=p.name;frames.append(j)
 c=json.loads((p/'command.json').read_text());expected=set(Path(c[c.index('--extract')+1]).read_text().split());audit.append({'chunk':p.name,'expected':len(expected),'joint_rows':len(j),'missing_ids':sorted(expected-set(j.ID)),'log_tail':(p/'console.log').read_text().splitlines()[-10:]})
d=pd.concat(frames,ignore_index=True);assert not d.ID.duplicated().any();d.to_csv(O/'all_joint_results.private.tsv.gz',sep='\t',index=False,compression='gzip')
valid=np.isfinite(d.CHISQ)&np.isfinite(d.LOG10P)&d.CHISQ.ge(0);v=d.loc[valid].copy();v['MAF']=np.minimum(v.A1FREQ,1-v.A1FREQ)
groups={'all':v,'excluding_NLGN4X_1Mb':v.loc[~v.GENPOS.between(5498724,6498724)],'MAF_below_0.01':v.loc[v.MAF.lt(.01)],'MAF_0.01_to_0.05':v.loc[v.MAF.ge(.01)&v.MAF.lt(.05)],'MAF_at_least_0.05':v.loc[v.MAF.ge(.05)]}
r=[];fig,axes=plt.subplots(2,3,figsize=(12,8))
for (name,z),ax in zip(groups.items(),axes.flat):
 n=len(z)
 if not n:continue
 lam=float(np.median(z.CHISQ)/chi2.ppf(.5,2));r.append({'group':name,'variants':n,'lambda_2df':lam,'p_lt_5e8':int(z.LOG10P.gt(-np.log10(5e-8)).sum()),'p_lt_5e9':int(z.LOG10P.gt(-np.log10(5e-9)).sum())})
 y=np.sort(z.LOG10P.to_numpy())[::-1];x=-np.log10((np.arange(n)+.5)/n);ax.plot(x,y,'.',ms=1);m=max(x.max(),y.max());ax.plot([0,m],[0,m],color='gray');ax.set(title=name+'; lambda='+format(lam,'.3f'),xlabel='Expected -log10(P)',ylabel='Observed -log10(P)')
axes.flat[-1].axis('off');fig.tight_layout();fig.savefig(O/'joint_X_QQ.pdf');fig.savefig(O/'joint_X_QQ.png',dpi=200)
(O/'calibration_review.json').write_text(json.dumps({'diagnostics':r,'weights':'none; all valid joint output rows; exact ranks before plotting','df':2,'chunk_audit':audit,'invalid_rows':int((~valid).sum()),'limitations':['QQ and lambda diagnostics do not prove calibration or replication']},indent=2));print('CALIBRATION_COMPLETE',flush=True)
