a<-commandArgs(TRUE)[1];t<-readRDS(file.path(a,"target_fit_details.private.rds"));s<-readRDS(file.path(a,"null_simulation.private.rds"));bad<-function(x)!is.null(x$error)||length(x$warnings)>0||length(x$p)!=1||!is.finite(x$p)
r<-capture.output({cat("target_records",length(t),"target_failures",sum(vapply(t,bad,logical(1))),"simulation_records",length(s),"simulation_failures",sum(vapply(s,bad,logical(1))),"
");cat("target_max_score",max(unlist(lapply(t,function(x)c(x$full_adjusted_score,x$null_adjusted_score)))),"simulation_max_score",max(unlist(lapply(s,function(x)c(x$full_adjusted_score,x$null_adjusted_score)))),"
");for(m in c("ML","Firth")){v<-Filter(function(x)identical(x$method,m),s);p<-vapply(v,function(x)x$p,numeric(1));k<-sum(p<.05);ci<-binom.test(k,length(p))$conf.int;cat(m,"replicates",length(p),"unique_replicates",length(unique(vapply(v,function(x)x$replicate,numeric(1)))),"nominal_05_rate",mean(p<.05),"CI",ci,"lambda_2df",median(qchisq(p,2,lower.tail=FALSE))/qchisq(.5,2),"
")};cat("markers",file.exists(file.path(a,c("TARGETS_FINISHED","SIMULATION_FINISHED","FINISHED"))),"
")});writeLines(r,file.path(a,"verified_audit.txt"));cat(r,sep="
")
