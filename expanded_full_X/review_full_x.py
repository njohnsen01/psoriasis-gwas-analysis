"""Strict completed-scan review and disclosure-limited common-variant outputs."""
from pathlib import Path
import pandas as pd,numpy as np,json,datetime,hashlib
from scipy.stats import chi2
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');S=R/'full_x';O=R/'aggregate_review'
def main():
 assert (S/'SCAN_COMPLETE.json').exists() and (R/'full_x_followup/FOLLOWUP_COMPLETE.json').exists();assert json.loads((R/'full_x_ploidy_audit/GLOBAL_PLOIDY_AUDIT.json').read_text())['status']=='passed';O.mkdir(exist_ok=True)
 x=pd.read_csv(S/'full_X_joint_Wald.tsv.gz',sep='\t');assert len(x)==334926 and x.ID.is_unique
 good=(x.ERRCODE=='.')&np.isfinite(x.WALD_CHISQ_3DF)&(x.WALD_CHISQ_3DF>=0)&np.isfinite(x.P_CHISQ_3DF);z=x[good].copy();assert len(z)>0
 assert np.allclose(chi2.sf(z.WALD_CHISQ_3DF,3),z.P_CHISQ_3DF,rtol=1e-10,atol=1e-15)
 fields=[c for c in ['#CHROM','POS','ID','REF','ALT','ALT1','A1','TEST','WALD_CHISQ_3DF','P_CHISQ_3DF','ERRCODE'] if c in x]
 safe=x[fields].copy();safe['estimated']=good
 # Exact per-variant N can reveal a small missing-call count when subtracted from cohort N.
 # Keep exact N only in protected outputs; exported bin has width 100 and is explicitly coarsened.
 safe['N_rounded_to_100']=(np.round(pd.to_numeric(x.OBS_CT,errors='coerce')/100)*100).astype('Int64')
 safe.to_csv(O/'full_X_joint_Wald.aggregate.tsv.gz',sep='\t',index=False,compression='gzip');safe[good].sort_values('P_CHISQ_3DF').head(100).to_csv(O/'full_X_top100.aggregate.tsv',sep='\t',index=False)
 f=pd.read_csv(R/'full_x_followup/followup_named.csv');f=f.merge(x[['ID','WALD_CHISQ_3DF']],on='ID',how='left',validate='many_to_one');f['independent_wald_absolute_difference']=abs(f.wald-f.WALD_CHISQ_3DF)
 assert (f.loc[f.status=='estimated','independent_wald_absolute_difference']<.005).all(),'Independent fit differs from full scan'
 fc=[c for c in ['ID','status','warnings','wald','ml_lrt','ml_lrt_p','firth_stat','firth_p','scan_p','independent_wald_absolute_difference'] if c in f];f[fc].to_csv(O/'full_X_lead_verification.aggregate.tsv',sep='\t',index=False)
 thresholds={str(t):int((z.P_CHISQ_3DF<t).sum()) for t in [.05,.01,1e-3,1e-4,1e-5,5e-8]};quantiles=[]
 for q in [.5,.9,.95,.99,.999]:quantiles.append({'quantile':q,'observed_W':float(z.WALD_CHISQ_3DF.quantile(q)),'null_W':float(chi2.ppf(q,3))})
 receipt={'status':'full_X_scan_and_lead_checks_completed_pending_interpretation','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'334926 previously QC-passed non-PAR X variants; PAR regions excluded','participants_before_variant_missingness':357530,'cases_before_variant_missingness':5324,'attempted':len(x),'estimated_strict':len(z),'unestimated':len(x)-len(z),'errors':x.ERRCODE.fillna('NA').value_counts().to_dict(),'lambda_median_3df':float(z.WALD_CHISQ_3DF.median()/chi2.ppf(.5,3)),'minimum_p':float(z.P_CHISQ_3DF.min()),'threshold_counts':thresholds,'distribution_quantiles':quantiles,'independent_followup_variants':len(f),'followup_failed':int((f.status!='estimated').sum()),'followup_Firth_below_5e8':int((f.firth_p<5e-8).sum()),'limitations':['Ordinary logistic joint Wald scan on provider-pruned cohort; no LOCO or mixed model in this full-X scan','Implementation validation and median lambda do not establish extreme-tail calibration','Selected lead follow-up is not independent replication','Per-variant sample size is coarsened to nearest 100 in exported results; exact counts remain protected']}
 (O/'FULL_X_REVIEW.json').write_text(json.dumps(receipt,indent=2)+'\n')
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,2,figsize=(12,4.7));p=np.sort(z.P_CHISQ_3DF.to_numpy());obs=-np.log10(np.maximum(p,np.nextafter(0.,1.)));exp=-np.log10((np.arange(1,len(p)+1)-.5)/len(p));sel=np.unique(np.r_[np.arange(min(2000,len(p))),np.linspace(0,len(p)-1,min(6000,len(p))).astype(int)])
 axs[0].scatter(exp[sel],obs[sel],s=5,alpha=.6,color='#345d85');mx=max(exp.max(),obs.max());axs[0].plot([0,mx],[0,mx],color='#888888',lw=1);axs[0].set(xlabel='Expected −log10(P)',ylabel='Observed −log10(P)',title=f'Non-PAR X: joint 3-df Wald\nMedian ratio = {receipt["lambda_median_3df"]:.3f}')
 axs[1].scatter(z.POS/1e6,-np.log10(np.maximum(z.P_CHISQ_3DF,np.nextafter(0.,1.))),s=2,alpha=.45,color='#345d85');axs[1].axhline(-np.log10(5e-8),ls='--',color='#ab403b',lw=1);axs[1].axvspan(5.7,6.3,color='#dfb651',alpha=.2);axs[1].set(xlabel='Chromosome X position (Mb; GRCh38)',ylabel='−log10(P)',title='Chromosome-wide association signal')
 fig.text(.01,.005,'Exploratory analysis; dashed line P = 5 × 10⁻⁸. Shading marks the NLGN4X region. Lead verification is reported separately.',fontsize=8);fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(O/'full_X_diagnostics.png',dpi=180);plt.close(fig)
 print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
