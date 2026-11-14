args <- commandArgs(TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1]); input <- normalizePath(args[2])
source(file.path(root,"chen_model.R"))
out <- file.path(root,"chen_target_run")
if(dir.exists(out)) stop("Output exists; do not overwrite an earlier run")
dir.create(out)
writeLines(c(as.character(Sys.time()),capture.output(sessionInfo())),file.path(out,"SESSION.txt"))
d <- read.csv(input,check.names=FALSE)
targets <- c("g0","g1","g2")
covs <- c("sex","age","age2","age_sex",paste0("PC",1:10))
stopifnot(identical(sort(names(d)),sort(c("psoriasis_case",covs,targets))),
          all(complete.cases(d[,c("psoriasis_case",covs)])),
          all(d$sex %in% c(0,1)))
Cbase <- cbind(intercept=1,as.matrix(d[,covs]))
Cprimary <- cbind(Cbase,sweep(as.matrix(d[,c(paste0("PC",1:10),"age2")]),1,d$sex,`*`))
colnames(Cprimary) <- c(colnames(Cbase),paste0("sex_",c(paste0("PC",1:10),"age2")))
matrices <- list(sex_specific_covariates=Cprimary,historical_covariates=Cbase)
writeLines(c("Target follow-up only; no genome-wide discovery inference.",
 "Primary: original nuisance covariates plus sex-by-PC1..10 and sex-by-age2.",
 "Historical covariate model is a sensitivity. Six distinct target/model tests reported.",
 "Input is the retained provider-pruned cohort. Residual relatedness remains possible.",
 "ML 3-df likelihood ratio. Exact formula and coding invariance checked on synthetic data.",
 "Nominal null simulation is not genome-wide-tail validation."),file.path(out,"MODEL.txt"))
work <- expand.grid(target=targets,model=names(matrices),stringsAsFactors=FALSE)
runone <- function(k) {
 target <- work$target[k]; model <- work$model[k]; C <- matrices[[model]]
 fits <- list(); records <- list()
 for(coding in c("02","01")) for(reverse in c(FALSE,TRUE)) {
   label <- paste(target,model,coding,reverse,sep="_")
   f <- tryCatch(chen_test(d$psoriasis_case,C,d[[target]],d$sex,coding,reverse),
                 error=function(e) list(error=conditionMessage(e)))
   fits[[label]] <- f
   records[[label]] <- data.frame(target=target,model=model,coding=coding,reverse=reverse,
    stat=if(is.null(f$stat)) NA_real_ else f$stat,
    df=if(is.null(f$df)) NA_integer_ else f$df,
    p=if(is.null(f$p)) NA_real_ else f$p,
    status=if(is.null(f$error)) "estimated" else f$error)
 }
 saveRDS(fits,file.path(out,paste0(target,"_",model,".private.rds")))
 r <- do.call(rbind,records)
 write.csv(r,file.path(out,paste0(target,"_",model,".csv")),row.names=FALSE)
 r
}
res <- do.call(rbind,parallel::mclapply(seq_len(nrow(work)),runone,mc.cores=6))
stopifnot(nrow(res)==24)
write.csv(res,file.path(out,"all_encoding_checks.csv"),row.names=FALSE)
valid <- all(res$status=="estimated")
if(valid) for(t in targets) for(m in names(matrices)) {
  z <- res[res$target==t & res$model==m,]
  if(max(z$stat)-min(z$stat)>1e-5) valid <- FALSE
}
if(!valid) {writeLines("Target fit or coding-invariance gate failed",file.path(out,"FAILED"));stop("Inspect recorded target failures")}
primary <- res[res$coding=="02" & !res$reverse,]
primary$p_holm_six_tests <- p.adjust(primary$p,method="holm",n=6)
write.csv(primary,file.path(out,"target_summary.csv"),row.names=FALSE)
print(primary,row.names=FALSE)
writeLines("Complete target pilot; NOT a chromosome-wide analysis",file.path(out,"TARGETS_FINISHED"))
