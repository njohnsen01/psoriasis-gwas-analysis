args <- commandArgs(TRUE)
stopifnot(length(args)==3)
root <- normalizePath(args[1]); input <- normalizePath(args[2]); lib <- normalizePath(args[3])
.libPaths(c(lib,.libPaths()))
stopifnot(requireNamespace("logistf",quietly=TRUE))
library(logistf)
out <- file.path(root,"chen_firth_run")
if(dir.exists(out)) stop("Output already exists")
dir.create(out)
writeLines(c(as.character(Sys.time()),capture.output(sessionInfo())),file.path(out,"SESSION.txt"))
d <- read.csv(input,check.names=FALSE)
covs <- c("sex","age","age2","age_sex",paste0("PC",1:10))
stopifnot(all(complete.cases(d[,c("psoriasis_case",covs)])),all(d$sex %in% c(0,1)))
for(v in c(paste0("PC",1:10),"age2")) d[[paste0("sex_",v)]] <- d$sex*d[[v]]
primary_covs <- c(covs,paste0("sex_",c(paste0("PC",1:10),"age2")))
work <- expand.grid(target=c("g0","g1","g2"),coding=c("02","01"),stringsAsFactors=FALSE)
runone <- function(k) {
  target <- work$target[k]; coding <- work$coding[k]; warnings <- character()
  ans <- withCallingHandlers(tryCatch({
    z <- d[is.finite(d[[target]]),]
    stopifnot(all(z[[target]] %in% 0:2),all(z[[target]][z$sex==1] %in% c(0,2)))
    z$g <- z[[target]]
    z$gd <- as.numeric(z$sex==0 & z$g==1)
    if(coding=="01") z$g[z$sex==1] <- z$g[z$sex==1]/2
    z$gs <- z$g*z$sex
    form <- reformulate(c(primary_covs,"g","gd","gs"),response="psoriasis_case")
    X <- model.matrix(form,z); y <- z$psoriasis_case
    stopifnot(qr(X)$rank==ncol(X))
    b <- logistf(form,data=z,pl=FALSE,control=logistf.control(maxit=200,maxhs=50,
      fit="IRLS",lconv=1e-8,gconv=1e-6,xconv=1e-8))
    idx <- match(c("g","gd","gs"),names(coef(b))); stopifnot(!anyNA(idx))
    test <- logistftest(b,test=idx,values=c(0,0,0))
    scorecheck <- function(beta,free) {
      p <- plogis(drop(X%*%beta)); w <- p*(1-p)
      info <- crossprod(X,X*w)
      h <- rowSums((X%*%solve(info))*X)*w
      max(abs(drop(crossprod(X,y-p+h*(0.5-p)))[free]))
    }
    fullscore <- scorecheck(coef(b),seq_len(ncol(X)))
    start <- coef(b);start[idx] <- 0;free <- setdiff(seq_len(ncol(X)),idx)
    # Same full-model Jeffreys penalty under both hypotheses; constrain only genetic terms.
    null <- logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=start,control=b$control,
      modcontrol=logistf.mod.control(terms.fit=free))
    nullscore <- scorecheck(null$beta,free)
    stat <- 2*(unname(b$loglik[1])-null$loglik)
    p <- pchisq(stat,3,lower.tail=FALSE)
    stopifnot(length(warnings)==0,stat>=0,is.finite(p),p>0,
      all(is.finite(b$conv)),all(abs(b$conv)<=c(1e-8,1e-6,1e-8)),
      fullscore<=1e-6,nullscore<=1e-6,null$iter<200,
      all(abs(null$conv[c(1,3)])<=c(1e-8,1e-8)),abs(log(p)-log(test$prob))<1e-4)
    list(stat=stat,p=p,full_adjusted_score=fullscore,null_adjusted_score=nullscore,
      iterations=b$iter,coefficients=coef(b),cells=table(z$sex,z[[target]],y))
  },error=function(e)list(error=conditionMessage(e))),warning=function(w){
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  saveRDS(ans,file.path(out,paste0(target,"_",coding,".private.rds")))
  r <- data.frame(target=target,coding=coding,method="Firth_3df_sex_specific_covariates",
    stat=if(is.null(ans$stat))NA_real_ else ans$stat,df=3,
    p=if(is.null(ans$p))NA_real_ else ans$p,
    status=if(is.null(ans$error))"estimated" else ans$error,
    warnings=paste(warnings,collapse="; "))
  write.csv(r,file.path(out,paste0(target,"_",coding,".csv")),row.names=FALSE)
  r
}
r <- do.call(rbind,parallel::mclapply(seq_len(nrow(work)),runone,mc.cores=6))
write.csv(r,file.path(out,"all_firth_checks.csv"),row.names=FALSE)
if(any(r$status!="estimated") || any(nzchar(r$warnings)))stop("Firth fit gate failed")
for(t in unique(r$target))stopifnot(diff(range(r$stat[r$target==t]))<1e-5)
s <- r[r$coding=="02",];s$p_holm_three_sensitivity_tests <- p.adjust(s$p,"holm",n=3)
write.csv(s,file.path(out,"firth_summary.csv"),row.names=FALSE)
print(s,row.names=FALSE)
writeLines("Completed Firth sensitivity on the three previously selected targets; no chromosome-wide calibration claim.",file.path(out,"FINISHED"))
