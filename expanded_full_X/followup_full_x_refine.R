a<-commandArgs(TRUE);out<-normalizePath(a[1]);.libPaths(c(normalizePath(a[2]),.libPaths()));library(logistf)
d<-read.csv(file.path(out,'selected.private.csv'));K<-as.integer(a[3]);stopifnot(all(d$CASE %in% 1:2),all(d$S %in% 0:1))
for(q in paste0('Q',1:24))d[[q]]<-d[[q]]/sd(d[[q]])
covs<-c('S',paste0('Q',1:24));d$y<-d$CASE-1
runone<-function(j){
 old<-file.path(a[4],paste0('g',j,'.csv'))
 if(file.exists(old)){r<-read.csv(old);if(identical(as.character(r$status),'estimated')){file.copy(old,file.path(out,basename(old)));file.copy(file.path(a[4],paste0('g',j,'.private.rds')),file.path(out,paste0('g',j,'.private.rds')));return(r)}}
 w<-character();ans<-withCallingHandlers(tryCatch({
 z<-d[is.finite(d[[paste0('g',j)]]),];z$g<-z[[paste0('g',j)]];z$gd<-as.numeric(z$S==0 & z$g==1);z$gs<-z$g*z$S
 form<-reformulate(c(covs,'g','gd','gs'),'y');X<-model.matrix(form,z);y<-z$y;ii<-(ncol(X)-2):ncol(X);stopifnot(qr(X)$rank==ncol(X))
 fit<-glm.fit(X,y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100));nullml<-glm.fit(X[,-ii],y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
 stopifnot(fit$converged,nullml$converged,max(abs(crossprod(X,y-fit$fitted.values)))<1e-5)
 V<-chol2inv(qr.R(fit$qr));pv<-order(fit$qr$pivot);V<-V[pv,pv];wald<-drop(t(fit$coefficients[ii])%*%solve(V[ii,ii],fit$coefficients[ii]));lrt<-nullml$deviance-fit$deviance
 f<-logistf(form,data=z,pl=FALSE,control=logistf.control(maxit=500,maxhs=100,fit='IRLS',lconv=1e-12,gconv=1e-10,xconv=1e-12))
 test<-logistftest(f,test=ii,values=c(0,0,0));start<-coef(f);start[ii]<-0;free<-setdiff(seq_len(ncol(X)),ii)
 nf<-logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=start,control=f$control,modcontrol=logistf.mod.control(terms.fit=free))
 score<-function(b,ix){p<-plogis(drop(X%*%b));ww<-p*(1-p);I<-crossprod(X,X*ww);h<-rowSums((X%*%solve(I))*X)*ww;max(abs(crossprod(X,y-p+h*(.5-p))[ix]))}
 fs<-score(coef(f),seq_len(ncol(X)));ns<-score(nf$beta,free);stat<-2*(unname(f$loglik[1])-nf$loglik);p<-pchisq(stat,3,lower.tail=FALSE)
 write.csv(data.frame(full_score=fs,null_score=ns,full_iter=max(f$iter),null_iter=max(nf$iter),f_lconv=f$conv[1],f_gconv=f$conv[2],f_xconv=f$conv[3]),file.path(out,paste0('g',j,'.numerics.private.csv')),row.names=FALSE)
 stopifnot(length(w)==0,fs<=1e-6,ns<=1e-6,stat>=0,f$iter<500,nf$iter<500,all(abs(f$conv)<=c(1e-8,1e-6,1e-8)),abs(log(p)-log(test$prob))<1e-4)
 list(wald=wald,ml_lrt=lrt,firth_stat=stat,firth_p=p,full_score=fs,null_score=ns,n=nrow(z),fit=f,private_cells=table(z$S,z$g,z$y))
 },error=function(e)list(error=conditionMessage(e))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 saveRDS(ans,file.path(out,paste0('g',j,'.private.rds')))
 r<-data.frame(key=paste0('g',j),status=if(is.null(ans$error))'estimated' else ans$error,warnings=paste(w,collapse='; '),wald=if(is.null(ans$wald))NA_real_ else ans$wald,ml_lrt=if(is.null(ans$ml_lrt))NA_real_ else ans$ml_lrt,ml_lrt_p=if(is.null(ans$ml_lrt))NA_real_ else pchisq(ans$ml_lrt,3,lower.tail=FALSE),firth_stat=if(is.null(ans$firth_stat))NA_real_ else ans$firth_stat,firth_p=if(is.null(ans$firth_p))NA_real_ else ans$firth_p)
 write.csv(r,file.path(out,paste0('g',j,'.csv')),row.names=FALSE);r
}
r<-do.call(rbind,parallel::mclapply(0:(K-1),runone,mc.cores=8));write.csv(r,file.path(out,'followup_summary.csv'),row.names=FALSE);print(r)
