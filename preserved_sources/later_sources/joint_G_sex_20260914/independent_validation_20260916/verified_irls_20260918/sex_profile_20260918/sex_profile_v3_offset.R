args <- commandArgs(TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'STARTED'))) stop('Existing run; will not duplicate')
writeLines(as.character(Sys.time()),file.path(out,'STARTED'))
.libPaths(c(file.path(input,'Rlib'),file.path(dirname(input),'Rlib'),.libPaths()))
library(logistf)
options(warn=2)
writeLines(capture.output(sessionInfo()),file.path(out,'SESSION_INFO.txt'))
d <- read.csv(file.path(input,'model_data.private.csv'))
old <- readRDS(file.path(input,'target_fit_details.private.rds'))
covs <- setdiff(names(d),c('psoriasis_case','g0','g1','g2'))
stopifnot('sex'%in%covs,all(d$sex%in%c(0,1)))
formula <- reformulate(c(covs,'female_g','male_carrier'),response='psoriasis_case')
ctrl <- logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-8,gconv=1e-6,xconv=1e-8)
strictctrl <- logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-10,gconv=1e-8,xconv=1e-10)
records <- list(); diagnostics <- list()
run <- function(){
for(v in c('g0','g1','g2')) {
  z <- d[!is.na(d[[v]]),]
  z$female_g <- z[[v]]*(1-z$sex)
  z$male_carrier <- (z[[v]]/2)*z$sex
  stopifnot(all(z$male_carrier%in%c(0,1)))
  # Native-ploidy and compensated inputs map to exactly the same carrier design.
  g01 <- z[[v]];g01[z$sex==1] <- g01[z$sex==1]/2
  stopifnot(identical(z$female_g,g01*(1-z$sex)),identical(z$male_carrier,g01*z$sex))
  X <- model.matrix(formula,z); y <- z$psoriasis_case
  prev <- old[[paste(v,'01','Firth',sep='_')]]
  stopifnot(is.null(prev$error),length(prev$warnings)==0)
  init <- prev$coef[match(colnames(X),names(prev$coef))];names(init)<-colnames(X)
  init['female_g'] <- prev$coef['g'];init['male_carrier']<-prev$coef['g']+prev$coef['gs']
  stopifnot(all(is.finite(init)))
  score <- function(beta,free){
    p<-plogis(drop(X%*%beta));w<-p*(1-p);I<-crossprod(X,X*w)
    h<-rowSums((X%*%solve(I))*X)*w
    max(abs(drop(crossprod(X,y-p+h*(.5-p)))[free]))
  }
  fit <- logistf(formula,data=z,pl=FALSE,init=init,control=ctrl)
  fs <- score(coef(fit),seq_len(ncol(X)))
  stopifnot(fs<=1e-6,all(abs(fit$conv)<=c(1e-8,1e-6,1e-8)))
  stopifnot(max(abs(coef(fit)-init))<1e-5)
  fullLL <- unname(fit$loglik[1])
  cache <- list()
  constrained <- function(idx,value,start=coef(fit)){
    start[idx]<-0;free<-setdiff(seq_len(ncol(X)),idx)
    fixed_offset<-drop(X[,idx,drop=FALSE]%*%value)
    a<-logistf:::logistf.fit(x=X,y=y,offset=fixed_offset,firth=TRUE,init=start,control=ctrl,
      modcontrol=logistf.mod.control(terms=free))
    a$beta[idx]<-value
    sc<-score(a$beta,free)
    if(sc>1e-6){
      saveRDS(list(target=v,idx=idx,value=value,fit=a,score=sc),file.path(out,'initial_profile_score_failure.private.rds'))
      cat('PROFILE_POLISH',v,'initial_score',sc,'\n');flush.console()
      restart<-a$beta;restart[idx]<-0
      a<-logistf:::logistf.fit(x=X,y=y,offset=fixed_offset,firth=TRUE,init=restart,control=strictctrl,
        modcontrol=logistf.mod.control(terms=free))
      a$beta[idx]<-value
      sc<-score(a$beta,free)
    }
    eta_check<-drop(X%*%a$beta);p_check<-plogis(eta_check)
    ll_check<-sum(y*eta_check-log1p(exp(eta_check)))+0.5*as.numeric(determinant(crossprod(X,X*(p_check*(1-p_check))),logarithm=TRUE)$modulus)
    stopifnot(abs(ll_check-a$loglik)<1e-6)
    saveRDS(list(target=v,idx=idx,value=value,fit=a,score=sc),file.path(out,'last_constrained_fit.private.rds'))
    stopifnot(all(is.finite(a$beta)),is.finite(a$loglik),a$iter<500,
      sc<=1e-6,all(abs(a$conv[c(1,3)])<=c(1e-8,1e-8)),
      max(abs(a$beta[idx]-value))<1e-12,a$loglik<=fullLL+1e-6)
    list(beta=a$beta,loglik=a$loglik,score=sc,iter=a$iter,conv=a$conv)
  }
  idxboth <- match(c('female_g','male_carrier'),colnames(X))
  jointnull <- constrained(idxboth,c(0,0))
  jp<-pchisq(max(0,2*(fullLL-jointnull$loglik)),2,lower.tail=FALSE)
  stopifnot(abs(log(jp)-log(prev$p))<1e-4)
  diagnostics[[v]]<-list(full_score=fs,full_conv=fit$conv,joint_p=jp,
    canonical_design_identical=TRUE,full_coefficients=coef(fit),joint_null=jointnull)
  for(term in c('female_g','male_carrier')) {
    idx<-match(term,colnames(X));est<-coef(fit)[idx];se<-sqrt(vcov(fit)[idx,idx])
    null<-constrained(idx,0)
    stat<-max(0,2*(fullLL-null$loglik));p<-pchisq(stat,1,lower.tail=FALSE)
    lt<-logistftest(fit,test=idx,values=0)
    stopifnot(is.finite(lt$prob),abs(log(p)-log(lt$prob))<1e-4)
    cache<-list(list(value=est,beta=coef(fit),loglik=fullLL,score=fs,iter=0))
    evaluate<-function(value){
      near<-which.min(vapply(cache,function(a)abs(a$value-value),numeric(1)))
      a<-constrained(idx,value,start=cache[[near]]$beta)
      a$value<-value;cache[[length(cache)+1]]<<-a
      2*(fullLL-a$loglik)-qchisq(.95,1)
    }
    bound<-function(direction){
      width<-2*se;endpoint<-est+direction*width;f<-evaluate(endpoint)
      expansions<-0L
      while(f<0 && expansions<6){width<-width*2;endpoint<-est+direction*width;f<-evaluate(endpoint);expansions<-expansions+1L}
      stopifnot(f>=0)
      interval<-sort(c(est,endpoint))
      r<-uniroot(evaluate,interval=interval,tol=1e-6,maxiter=40)$root
      stopifnot(abs(evaluate(r))<1e-4)
      r
    }
    lower<-bound(-1);upper<-bound(1)
    key<-paste(v,term,sep='_')
    records[[key]]<-data.frame(target=v,sex=if(term=='female_g')'female' else 'male',
      effect_unit=if(term=='female_g')'per effect allele' else 'carrier versus noncarrier',
      beta=unname(est),OR=exp(unname(est)),profile_CI_low=exp(lower),profile_CI_high=exp(upper),
      PLRT_statistic=stat,PLRT_p=p,p_bonferroni_6=min(1,6*p),joint_p=jp,
      full_score=fs,constrained_score=null$score,
      profile_max_score=max(vapply(cache,function(a)a$score,numeric(1))),
      profile_evaluations=length(cache),status='verified')
    diagnostics[[key]]<-list(null=null,profile=cache,logistftest=lt)
    saveRDS(diagnostics,file.path(out,'profile_diagnostics.private.rds'))
    write.csv(do.call(rbind,records),file.path(out,'sex_profile_results.csv'),row.names=FALSE)
    cat('VERIFIED',v,term,'P',format(p,digits=10),'CI',exp(lower),exp(upper),'\n');flush.console()
  }
}
stopifnot(length(records)==6)
writeLines('Six sex-specific tests and profile intervals verified; canonical designs agree across encodings.',file.path(out,'FINISHED'))
}
tryCatch(run(),error=function(e){writeLines(conditionMessage(e),file.path(out,'FAILED'));stop(e)})

