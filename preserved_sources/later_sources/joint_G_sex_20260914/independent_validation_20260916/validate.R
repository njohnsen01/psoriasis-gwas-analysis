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
 b<-logistf(f1,data=z,pl=FALSE,control=logistf.control(maxit=100,maxhs=50,lconv=1e-8,gconv=1e-6,xconv=1e-8));idx<-match(c("g","gs"),names(coef(b)));stopifnot(!anyNA(idx));test<-logistftest(b,test=idx,values=c(0,0));list(p=as.numeric(test$prob),coef=coef(b),var=b$var,convergence=b$conv,test=test)
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
writeLines("complete",file.path(root,"TARGETS_FINISHED"))
# Fixed representative target chosen in advance; parametric-null pilot, not extreme-tail validation.
z<-d[!is.na(d$g2),];z$g<-z$g2;if(any(z$sex==1))z$g[z$sex==1]<-z$g[z$sex==1]/2;z$gs<-z$g*z$sex
null<-glm(f0,data=z,family=binomial(),control=glm.control(maxit=100,epsilon=1e-9));stopifnot(null$converged);prob<-fitted(null);set.seed(20260916);sim<-list()
for(i in seq_len(100)){
 z$psoriasis_case<-rbinom(nrow(z),1,prob)
 for(method in c("ML","Firth")){
 a<-fitone(z,method);key<-paste(i,method,sep="_");sim[[key]]<-list(replicate=i,method=method,p=a$p,error=a$error,warnings=a$warnings,convergence=a$convergence)
 }
 saveRDS(sim,file.path(root,"null_simulation.private.rds"));cat("NULL_REPLICATE_DONE",i,"\n");flush.console()
}
writeLines("100 replicates; inspect failures and warnings before inference",file.path(root,"SIMULATION_FINISHED"));writeLines("complete",file.path(root,"FINISHED"))
