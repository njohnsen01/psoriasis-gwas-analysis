args<-commandArgs(TRUE); root<-args[1]; .libPaths(c(file.path(root,"Rlib"),.libPaths()));library(logistf)
out<-file.path(root,"convergence_pilot_20260918");if(!dir.create(out))stop("Pilot directory already exists")
writeLines(capture.output(sessionInfo()),file.path(out,"SESSION_INFO.txt"))
d<-read.csv(file.path(root,"model_data.private.csv")); covs<-setdiff(names(d),c("psoriasis_case","g0","g1","g2"));z<-d[!is.na(d$g0),];z$g<-z$g0;z$g[z$sex==1]<-z$g[z$sex==1]/2;z$gs<-z$g*z$sex
f<-reformulate(c(covs,"g","gs"),response="psoriasis_case")
old<-readRDS(file.path(root,"target_fit_details.private.rds"))[["g0_01_Firth"]]
ctrl<-logistf.control(maxit=500,maxhs=50,lconv=1e-8,gconv=1e-6,xconv=1e-8)
w<-character();cat("PILOT_STARTED\n");flush.console()
b<-withCallingHandlers(logistf(f,data=z,pl=FALSE,init=old$coef,control=ctrl),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")})
saveRDS(b,file.path(out,"fit.private.rds"));cat("FIT_CONVERGENCE\n");print(b$conv);cat("ITER\n");print(b$iter);cat("WARNINGS\n");print(w);flush.console()
ok<-length(w)==0&&all(is.finite(b$conv))&&all(abs(b$conv)<=c(1e-8,1e-6,1e-8))
writeLines(paste("fit_pass",ok),file.path(out,"STATUS.txt"))
if(ok){tw<-character();test<-withCallingHandlers(logistftest(b,test=match(c("g","gs"),names(coef(b))),values=c(0,0),control=ctrl),warning=function(e){tw<<-c(tw,conditionMessage(e));invokeRestart("muffleWarning")});saveRDS(list(test=test,warnings=tw),file.path(out,"test.private.rds"));cat("JOINT_P",test$prob,"\n");print(tw)}
writeLines("Pilot finished; inspect convergence, not completion alone",file.path(out,"FINISHED"))
