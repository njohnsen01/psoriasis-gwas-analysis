a<-commandArgs(TRUE);out<-normalizePath(a[1]);.libPaths(c(normalizePath(a[2]),.libPaths()));library(logistf)
d<-read.csv(file.path(out,'selected.private.csv'));sel<-read.csv(file.path(out,'selection.private.csv'));for(q in paste0('Q',1:24))d[[q]]<-d[[q]]/sd(d[[q]]);d$y<-d$CASE-1;covs<-c('S',paste0('Q',1:24))
runone<-function(j){
 key<-paste0('g',j);dest<-file.path(out,paste0(key,'.csv'));if(file.exists(dest))return(read.csv(dest));w<-character()
 ans<-withCallingHandlers(tryCatch({
 z<-d[is.finite(d[[key]]),];z$g<-z[[key]];z$gd<-as.numeric(z$S==0 & z$g==1);z$gs<-z$S*z$g;form<-reformulate(c(covs,'g','gd','gs'),'y');X<-model.matrix(form,z);y<-z$y;ii<-(ncol(X)-2):ncol(X);free<-setdiff(seq_len(ncol(X)),ii);stopifnot(qr(X)$rank==ncol(X));wald<-NA_real_
 if(sel$ml_eligible[j+1]){
  ml<-glm.fit(X,y,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200));stopifnot(ml$converged,max(abs(crossprod(X,y-ml$fitted.values)))<=1e-6);ww<-ml$fitted.values*(1-ml$fitted.values);V<-chol2inv(chol(crossprod(X,X*ww)));wald<-drop(t(ml$coefficients[ii])%*%solve(V[ii,ii],ml$coefficients[ii]))
 }
 f<-logistf(form,data=z,pl=FALSE,control=logistf.control(maxit=500,maxhs=100,fit='IRLS',lconv=1e-12,gconv=1e-10,xconv=1e-12));test<-logistftest(f,test=ii,values=c(0,0,0));start<-coef(f);start[ii]<-0
 nf<-logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=start,control=f$control,modcontrol=logistf.mod.control(terms.fit=free))
 score<-function(b,ix){p<-plogis(drop(X%*%b));ww<-p*(1-p);I<-crossprod(X,X*ww);h<-rowSums((X%*%solve(I))*X)*ww;max(abs(crossprod(X,y-p+h*(.5-p))[ix]))}
 fs<-score(coef(f),seq_len(ncol(X)));ns<-score(nf$beta,free);stat<-2*(unname(f$loglik[1])-nf$loglik);p<-pchisq(stat,3,lower.tail=FALSE)
 write.csv(data.frame(full_score=fs,null_score=ns,full_iterations=max(f$iter),null_iterations=max(nf$iter)),file.path(out,paste0(key,'.numerics.private.csv')),row.names=FALSE)
 stopifnot(length(w)==0,fs<=1e-6,ns<=1e-6,stat>=0,f$iter<500,nf$iter<500,all(abs(f$conv)<=c(1e-8,1e-6,1e-8)),abs(log(p)-log(test$prob))<1e-4)
 list(wald=wald,firth_stat=stat,firth_p=p,full_score=fs,null_score=ns,fit=f,private_cells=table(z$S,z$g,z$y))
 },error=function(e)list(error=conditionMessage(e))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 saveRDS(ans,file.path(out,paste0(key,'.private.rds')))
 r<-data.frame(key=key,ID=sel$ID[j+1],status=if(is.null(ans$error))'estimated' else ans$error,warnings=paste(w,collapse='; '),wald=if(is.null(ans$wald))NA_real_ else ans$wald,firth_stat=if(is.null(ans$firth_stat))NA_real_ else ans$firth_stat,firth_p=if(is.null(ans$firth_p))NA_real_ else ans$firth_p)
 write.csv(r,dest,row.names=FALSE);r
}
r<-do.call(rbind,parallel::mclapply(0:(nrow(sel)-1),runone,mc.cores=8));write.csv(r,file.path(out,'native_summary.csv'),row.names=FALSE)
