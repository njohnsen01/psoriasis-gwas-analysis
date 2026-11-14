from pathlib import Path
b=Path('/home/jupyter/psoriasis_manuscript_cohort_calibration_2026-09-09/replacement_primary_20260911/reconciliation_20260912/joint_G_sex_20260914/independent_validation_20260916')
o=b/'verified_irls_20260918';o.mkdir(exist_ok=False)
for n in ['Rlib','model_data.private.csv']:(o/n).symlink_to(b/n,target_is_directory=n=='Rlib')
s=(b/'validate.R').read_text()
s=s.replace('maxit=100,maxhs=50,lconv','maxit=200,maxhs=50,fit="IRLS",lconv')
old='test<-logistftest(b,test=idx,values=c(0,0));list('
new='''test<-logistftest(b,test=idx,values=c(0,0));
X<-model.matrix(f1,z); yy<-z$psoriasis_case
scorecheck<-function(bb,free){pp<-plogis(drop(X%*%bb));ww<-pp*(1-pp);info<-crossprod(X,X*ww);hh<-rowSums((X%*%solve(info))*X)*ww;max(abs(drop(crossprod(X,yy-pp+hh*(0.5-pp)))[free]))}
fullscore<-scorecheck(coef(b),seq_len(ncol(X)))
start<-coef(b);start[idx]<-0;free<-setdiff(seq_len(ncol(X)),idx)
nfit<-logistf:::logistf.fit(x=X,y=yy,firth=TRUE,init=start,control=b$control,modcontrol=logistf.mod.control(terms.fit=free))
nullscore<-scorecheck(nfit$beta,free)
checkp<-pchisq(2*(unname(b$loglik[1])-nfit$loglik),2,lower.tail=FALSE)
stopifnot(length(w)==0,all(is.finite(b$conv)),all(abs(b$conv)<=c(1e-8,1e-6,1e-8)),fullscore<=1e-6,nullscore<=1e-6,nfit$iter<200,all(abs(nfit$conv[c(1,3)])<=c(1e-8,1e-8)),is.finite(checkp),abs(log(checkp)-log(test$prob))<1e-4)
list(full_adjusted_score=fullscore,null_adjusted_score=nullscore,verified_joint_p=checkp,'''
assert old in s
s=s.replace(old,new)
s=s.replace('writeLines("complete",file.path(root,"TARGETS_FINISHED"))','''if(any(vapply(details,function(x)!is.null(x$error)||length(x$warnings)>0,logical(1))))stop("Target convergence gate failed; simulations not started")
writeLines("verified complete",file.path(root,"TARGETS_FINISHED"))''')
s=s.replace('convergence=a$convergence)','convergence=a$convergence,full_adjusted_score=a$full_adjusted_score,null_adjusted_score=a$null_adjusted_score)')
s=s.replace('writeLines("100 replicates; inspect failures and warnings before inference",file.path(root,"SIMULATION_FINISHED"))','''if(length(sim)!=200 || any(vapply(sim,function(x)!is.null(x$error)||length(x$warnings)>0||!is.finite(x$p),logical(1))))stop("Simulation verification failed; inspect saved records")
writeLines("100 replicates verified; not genome-wide tail validation",file.path(root,"SIMULATION_FINISHED"))''')
(o/'validate.R').write_text(s)
print('Created verified_irls_20260918 with explicit full and constrained score gates')
