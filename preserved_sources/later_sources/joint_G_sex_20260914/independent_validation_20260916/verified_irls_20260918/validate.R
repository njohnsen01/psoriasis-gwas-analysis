args<-commandArgs(TRUE);root<-args[1];.libPaths(c(file.path(root,"Rlib"),.libPaths()))
if(file.exists(file.path(root,"ANALYSIS_STARTED")))stop("Already started")
writeLines(as.character(Sys.time()),file.path(root,"ANALYSIS_STARTED"))
for(i in seq_len(120)){if(requireNamespace("logistf",quietly=TRUE))break;Sys.sleep(15)}
if(!requireNamespace("logistf",quietly=TRUE))stop("logistf installation unavailable; no fits performed")
library(logistf);writeLines(capture.output(sessionInfo()),file.path(root,"SESSION_INFO.txt"))
d<-read.csv(file.path(root,"model_data.private.csv"));covs<-setdiff(names(d),c("psoriasis_case","g0","g1","g2"));f0<-reformulate(covs,response="psoriasis_case");f1<-reformulate(c(covs,"g","gs"),response="psoriasis_case")
records<-list();details<-list()
fitone<-function(z,method){
 w<-character();value<-withCallingHandlers(tryCatch({
 if(method=="ML"){
 a<-glm(f0,data=z,family=binomial(),control=glm.control(maxit=100,epsilon=1e-9));b<-glm(f1,data=z,family=binomial(),control=glm.control(maxit=100,epsilon=1e-9));stopifnot(a$converged,b$converged,b$rank-a$rank==2);stat<-2*(as.numeric(logLik(b))-as.numeric(logLik(a)));stopifnot(stat>= -1e-7);list(p=pchisq(max(0,stat),2,lower.tail=FALSE),stat=stat,coef=coef(b),var=vcov(b),convergence=c(a$converged,b$converged))
 }else{
 b<-logistf(f1,data=z,pl=FALSE,control=logistf.control(maxit=200,maxhs=50,fit="IRLS",lconv=1e-8,gconv=1e-6,xconv=1e-8));idx<-match(c("g","gs"),names(coef(b)));stopifnot(!anyNA(idx));test<-logistftest(b,test=idx,values=c(0,0));
X<-model.matrix(f1,z); yy<-z$psoriasis_case
scorecheck<-function(bb,free){pp<-plogis(drop(X%*%bb));ww<-pp*(1-pp);info<-crossprod(X,X*ww);hh<-rowSums((X%*%solve(info))*X)*ww;max(abs(drop(crossprod(X,yy-pp+hh*(0.5-pp)))[free]))}
fullscore<-scorecheck(coef(b),seq_len(ncol(X)))
start<-coef(b);start[idx]<-0;free<-setdiff(seq_len(ncol(X)),idx)
nfit<-logistf:::logistf.fit(x=X,y=yy,firth=TRUE,init=start,control=b$control,modcontrol=logistf.mod.control(terms.fit=free))
nullscore<-scorecheck(nfit$beta,free)
checkp<-pchisq(2*(unname(b$loglik[1])-nfit$loglik),2,lower.tail=FALSE)
stopifnot(length(w)==0,all(is.finite(b$conv)),all(abs(b$conv)<=c(1e-8,1e-6,1e-8)),fullscore<=1e-6,nullscore<=1e-6,nfit$iter<200,all(abs(nfit$conv[c(1,3)])<=c(1e-8,1e-8)),is.finite(checkp),abs(log(checkp)-log(test$prob))<1e-4)
list(full_adjusted_score=fullscore,null_adjusted_score=nullscore,verified_joint_p=checkp,p=as.numeric(test$prob),coef=coef(b),var=b$var,convergence=b$conv,test=test)
 }
 },error=function(e)list(error=conditionMessage(e))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")});value$warnings<-w;value
}
for(v in c("g0","g1","g2")){
 z<-d[!is.na(d[[v]]),];for(coding in c("01","02")){
 z$g<-z[[v]];if(coding=="01")z$g[z$sex==1]<-z$g[z$sex==1]/2;z$gs<-z$g*z$sex
 for(method in c("ML","Firth")){
 key<-paste(v,coding,method,sep="_");a<-fitone(z,method);details[[key]]<-a;records[[key]]<-data.frame(target=v,coding=coding,method=method,p=if(is.null(a$p))NA_real_ else a$p,error=if(is.null(a$error))"" else a$error,warnings=paste(a$warnings,collapse="; "));saveRDS(details,file.path(root,"target_fit_details.private.rds"));write.csv(do.call(rbind,records),file.path(root,"target_results.csv"),row.names=FALSE);cat("TARGET_DONE",key,"\n");flush.console()
 }
 }}
if(any(vapply(details,function(x)!is.null(x$error)||length(x$warnings)>0,logical(1))))stop("Target convergence gate failed; simulations not started")
writeLines("verified complete",file.path(root,"TARGETS_FINISHED"))
# Fixed representative target chosen in advance; parametric-null pilot, not extreme-tail validation.
z<-d[!is.na(d$g2),];z$g<-z$g2;if(any(z$sex==1))z$g[z$sex==1]<-z$g[z$sex==1]/2;z$gs<-z$g*z$sex
null<-glm(f0,data=z,family=binomial(),control=glm.control(maxit=100,epsilon=1e-9));stopifnot(null$converged);prob<-fitted(null);set.seed(20260916);sim<-list()
for(i in seq_len(100)){
 z$psoriasis_case<-rbinom(nrow(z),1,prob)
 for(method in c("ML","Firth")){
 a<-fitone(z,method);key<-paste(i,method,sep="_");sim[[key]]<-list(replicate=i,method=method,p=a$p,error=a$error,warnings=a$warnings,convergence=a$convergence,full_adjusted_score=a$full_adjusted_score,null_adjusted_score=a$null_adjusted_score)
 }
 saveRDS(sim,file.path(root,"null_simulation.private.rds"));cat("NULL_REPLICATE_DONE",i,"\n");flush.console()
}
if(length(sim)!=200 || any(vapply(sim,function(x)!is.null(x$error)||length(x$warnings)>0||!is.finite(x$p),logical(1))))stop("Simulation verification failed; inspect saved records")
writeLines("100 replicates verified; not genome-wide tail validation",file.path(root,"SIMULATION_FINISHED"));writeLines("complete",file.path(root,"FINISHED"))
