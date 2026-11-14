import os
from pathlib import Path
import sys,pandas as pd,numpy as np,subprocess,json,hashlib,datetime
sys.path.insert(0,f'{os.environ["PGENLIB_SITE_DIR"]}')
import pgenlib
R=Path('/home/jupyter/x_followup_20261013/continuation_20261013');D=R/'proxy_benchmark';Q=R.parent/'chen_scan_benchmark'
B=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/X_projected_source/chunk_0009/input')
s=pd.read_csv(str(B)+'.psam',sep=r'\s+',dtype={'#FID':str,'FID':str,'IID':str});s.rename(columns={'#FID':'FID'},inplace=True)
c=pd.read_csv(Q/'covariates.private.tsv',sep='\t',dtype={'FID':str,'IID':str});y=pd.read_csv(Q/'phenotype.private.tsv',sep='\t',dtype={'FID':str,'IID':str})
d=s[['FID','IID','SEX']].assign(sample_index=np.arange(len(s))).merge(c,on=['FID','IID'],validate='one_to_one').merge(y,on=['FID','IID'],validate='one_to_one');d=d.sort_values('sample_index')
assert len(d)==357530 and np.array_equal(d.S,(d.SEX==1).astype(int))
pv=pd.read_csv(str(B)+'.pvar',sep='\t',skiprows=sum(1 for l in Path(str(B)+'.pvar').open() if l.startswith('##')));assert 'ID' in pv
ex=(Q/'extract100.txt').read_text().split();targets=['chrX:5995921:G:A','chrX:5996494:GCT:G','chrX:5998724:A:G'];extra=sorted(set(ex)-set(targets),key=lambda x:hashlib.sha256(('validation|'+x).encode()).digest())[:3];chosen=targets+extra
reader=pgenlib.PgenReader(str(B.with_suffix('.pgen')).encode(),raw_sample_ct=len(s));arr=np.empty(len(s),np.int8);audit=[]
for j,vid in enumerate(chosen):
 ix=int(pv.index[pv.ID==vid][0]);reader.read(ix,arr);g=arr[d.sample_index.to_numpy()].astype(float);g[g<0]=np.nan
 assert set(np.unique(g[np.isfinite(g)&(d.S.to_numpy()==1)]))<=set([0.,2.])
 d[f'g{j}']=g;audit.append({'key':f'g{j}','ID':vid,'n':int(np.isfinite(g).sum()),'male_nonmissing_encoding':'0/2_only'})
reader.close();d[['CASE','S']+[f'Q{i}' for i in range(1,25)]+[f'g{i}' for i in range(6)]].to_csv(D/'validation.private.csv',index=False)
(D/'validation_ids.json').write_text(json.dumps(audit,indent=2))
r='''d<-read.csv("validation.private.csv"); C<-as.matrix(cbind(intercept=1,d[,c("S",paste0("Q",1:24))])); y<-d$CASE-1; out<-list(); for(j in 0:5) for(coding in c("02","01")) for(reverse in c(FALSE,TRUE)){ g<-d[[paste0("g",j)]];k<-is.finite(g);g<-g[k];sex<-d$S[k];if(reverse)g<-2-g; dom<-as.numeric(g==1 & sex==0);if(coding=="01")g[sex==1]<-g[sex==1]/2;X<-cbind(C[k,],ADD=g,DOMDEV=dom,ADDS=g*sex); w<-character(); f<-withCallingHandlers(glm.fit(X,y[k],family=binomial(),control=glm.control(epsilon=1e-10,maxit=100)),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")});stopifnot(f$converged,length(w)==0,f$rank==ncol(X)); V<-chol2inv(qr.R(f$qr)); b<-f$coefficients; invpivot<-order(f$qr$pivot);V<-V[invpivot,invpivot];ii<-(ncol(X)-2):ncol(X);W<-as.numeric(t(b[ii])%*%solve(V[ii,ii],b[ii]));score<-max(abs(crossprod(X,y[k]-f$fitted.values)));stopifnot(score<1e-5);out[[length(out)+1]]<-data.frame(key=paste0("g",j),coding=coding,reverse=reverse,n=sum(k),wald_chisq=W,p_chisq=pchisq(W,3,lower.tail=FALSE),p_F=pf(W/3,3,sum(k)-ncol(X),lower.tail=FALSE),score=score)};write.csv(do.call(rbind,out),"independent_wald_checks.csv",row.names=FALSE);print(do.call(rbind,out)[c(1,5,9,13,17,21),],row.names=FALSE)'''
(D/'validate_independent.R').write_text(r)
with (D/'validation.console.log').open('w') as f:p=subprocess.run(['/opt/conda/envs/jupyter/bin/Rscript','validate_independent.R'],cwd=D,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,'Independent R check failed'
v=pd.read_csv(D/'independent_wald_checks.csv');z=pd.read_csv(D/'fit.CASE.glm.logistic',sep=r'\s+');checks=[]
for a in audit:
 vs=v[v.key==a['key']];p0=z[(z.ID==a['ID'])&(z.TEST=='USER_3DF')].iloc[0]
 delta=abs(float(p0.Z_OR_F_STAT)*3-vs.wald_chisq.iloc[0]);spread=vs.wald_chisq.max()-vs.wald_chisq.min()
 assert int(p0.OBS_CT)==a['n']==int(vs.n.iloc[0]) and delta<1e-3 and spread<1e-5,(a,delta,spread)
 checks.append({**a,'statistic_absolute_difference':delta,'coding_reversal_spread':spread,'plink_p':float(p0.P),'R_chisq_p':float(vs.p_chisq.iloc[0]),'R_F_p':float(vs.p_F.iloc[0])})
receipt={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'passed','source':'Original X pgen read with pgenlib; independent R glm on aligned participants','checks':checks,'limitations':'No mixed model or sparse-cell calibration established by implementation checks'}
(D/'VALIDATED.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
