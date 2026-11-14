a<-commandArgs(TRUE);out<-normalizePath(a[1]);d<-read.csv(file.path(out,'selected.private.csv'));sel<-read.csv(file.path(out,'selection.private.csv'));d$y<-d$CASE-1;for(q in paste0('Q',1:24))d[[q]]<-d[[q]]/sd(d[[q]]);covs<-c('S',paste0('Q',1:24))
runone<-function(j){
 key<-paste0('g',j);dest<-file.path(out,paste0(key,'.csv'));if(file.exists(dest))return(read.csv(dest));w<-character();res<-withCallingHandlers(tryCatch({
 z<-d[is.finite(d[[key]]),];z$g<-z[[key]];z$gd<-as.numeric(z$S==0 & z$g==1);z$gs<-z$g*z$S;X<-model.matrix(reformulate(c(covs,'g','gd','gs')),z);ii<-(ncol(X)-2):ncol(X);fit<-glm.fit(X,z$y,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200));score<-max(abs(crossprod(X,z$y-fit$fitted.values)));stopifnot(length(w)==0,fit$converged,fit$rank==ncol(X),score<=1e-6);V<-chol2inv(qr.R(fit$qr));pv<-order(fit$qr$pivot);V<-V[pv,pv];wald<-drop(t(fit$coefficients[ii])%*%solve(V[ii,ii],fit$coefficients[ii]));list(wald=wald,score=score)
 },error=function(e)list(error=conditionMessage(e))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 r<-data.frame(ID=sel$ID[j+1],status=if(is.null(res$error))'estimated' else res$error,warnings=paste(w,collapse=';'),wald=if(is.null(res$wald))NA_real_ else res$wald,score=if(is.null(res$score))NA_real_ else res$score);write.csv(r,dest,row.names=FALSE);r
}
r<-do.call(rbind,parallel::mclapply(0:(nrow(sel)-1),runone,mc.cores=8));write.csv(r,file.path(out,'R_summary.csv'),row.names=FALSE)
