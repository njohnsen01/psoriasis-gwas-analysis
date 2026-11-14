args <- commandArgs(TRUE); out <- normalizePath(args[1]); v <- normalizePath(file.path(out,'../joint_G_sex_20260914/independent_validation_20260916'))
.libPaths(c(file.path(v,'Rlib'),.libPaths())); library(logistf)
writeLines(capture.output(sessionInfo()),file.path(out,'SESSION_INFO.txt'))
d <- read.csv(file.path(v,'verified_irls_20260918/model_data.private.csv'))
stopifnot(all(d$sex %in% c(0,1))); covs <- c('age','age2',paste0('PC',1:10))
f0 <- reformulate(covs,response='psoriasis_case'); f1 <- reformulate(c(covs,'g'),response='psoriasis_case')
ctrl <- logistf.control(maxit=300,maxhs=50,fit='IRLS',lconv=1e-8,gconv=1e-6,xconv=1e-8)
records <- list(); fits <- list()
for (target in c('g0','g1','g2')) {
 z <- d[d$sex==0 & !is.na(d[[target]]),]; z$g <- z[[target]]
 stopifnot(all(complete.cases(z[,c('psoriasis_case',covs,'g')]))); X <- model.matrix(f1,z); stopifnot(qr(X)$rank==ncol(X))
 for (method in c('ML','Firth')) {
  warnings <- character(); key <- paste(target,method,sep='_')
  result <- withCallingHandlers(tryCatch({
   if (method=='ML') {
    fit <- glm(f1,z,family=binomial(),control=glm.control(maxit=150,epsilon=1e-10)); null <- glm(f0,z,family=binomial(),control=glm.control(maxit=150,epsilon=1e-10)); stopifnot(fit$converged,null$converged)
    beta <- coef(fit); covariance <- vcov(fit); stat <- 2*(as.numeric(logLik(fit))-as.numeric(logLik(null))); stopifnot(stat>=-1e-8); pvalue <- pchisq(max(stat,0),1,lower.tail=FALSE); testconv <- TRUE
   } else {
    fit <- logistf(f1,data=z,pl=FALSE,control=ctrl); beta <- coef(fit); covariance <- fit$var
    test <- logistftest(fit,test=match('g',names(beta)),values=0,control=ctrl); pvalue <- as.numeric(test$prob); testconv <- all(is.finite(fit$conv)) && all(abs(fit$conv)<=c(1e-8,1e-6,1e-8))
   }
   fittedp <- plogis(drop(X %*% beta)); w <- fittedp*(1-fittedp); adjustment <- rep(0,nrow(X))
   if(method=='Firth') { information <- crossprod(X,X*w); hat <- rowSums((X %*% solve(information))*X)*w; adjustment <- hat*(.5-fittedp) }
   scoremax <- max(abs(crossprod(X,z$psoriasis_case-fittedp+adjustment)))
   stopifnot(is.finite(pvalue),pvalue>=0,pvalue<=1,scoremax<1e-5,testconv)
   j <- match('g',names(beta)); se <- sqrt(covariance[j,j]); fits[[key]] <- fit
   data.frame(target=target,method=method,beta=beta[j],se=se,OR=exp(beta[j]),wald_lower=exp(beta[j]-qnorm(.975)*se),wald_upper=exp(beta[j]+qnorm(.975)*se),p=pvalue,max_adjusted_score=scoremax,convergence_pass=TRUE,error='')
  },error=function(e)data.frame(target=target,method=method,beta=NA,se=NA,OR=NA,wald_lower=NA,wald_upper=NA,p=NA,max_adjusted_score=NA,convergence_pass=FALSE,error=conditionMessage(e))),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  result$warnings <- paste(warnings,collapse='; '); records[[key]] <- result
  write.csv(do.call(rbind,records),file.path(out,'RESULTS.csv'),row.names=FALSE);saveRDS(fits,file.path(out,'fits.private.rds'));cat('COMPLETED',key,fill=TRUE);flush.console()
 }
}
allresults <- do.call(rbind,records);stopifnot(all(allresults$convergence_pass),all(allresults$warnings==''));writeLines('All six bounded fits completed; source and input verification still required',file.path(out,'FITS_COMPLETE'))
