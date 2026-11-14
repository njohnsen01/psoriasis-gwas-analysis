args <- commandArgs(TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[1]); out <- normalizePath(args[2])
stopifnot(file.exists(file.path(out,'PROTOCOL.json')))
if(file.exists(file.path(out,'MODEL_STARTED'))) stop('Run already started; preserve existing results')
writeLines(as.character(Sys.time()),file.path(out,'MODEL_STARTED'))
.libPaths(c(file.path(dirname(input),'Rlib'),.libPaths()));library(logistf)
writeLines(capture.output(sessionInfo()),file.path(out,'MODEL_SESSION.txt'))
d<-read.csv(file.path(input,'model_data.private.csv'))
previous<-readRDS(file.path(input,'target_fit_details.private.rds'))
female<-read.csv(file.path(dirname(dirname(dirname(input))),'targeted_model_comparison_20260919/RESULTS.csv'))
covs<-setdiff(names(d),c('psoriasis_case','g0','g1','g2'))
stopifnot(identical(covs,c('sex','age','age2','age_sex',paste0('PC',1:10))))
ctrl<-logistf.control(maxit=300,maxhs=50,fit='IRLS',lconv=1e-8,gconv=1e-6,xconv=1e-8)
strict<-logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-10,gconv=1e-8,xconv=1e-10)
models<-list(A_original=character(),B_age2='sex_age2',C_PCs=paste0('sex_PC',1:10),D_all=c('sex_age2',paste0('sex_PC',1:10)))
records<-list();diagnostics<-list();objects<-list();nested<-list()
for(t in c('g0','g1','g2')){
 z<-d[!is.na(d[[t]]),];stopifnot(all(z$sex %in% c(0,1)))
 z$gf<-z[[t]]*(1-z$sex);z$gm<-z[[t]]*z$sex/2
 z$sex_age2<-z$sex*z$age2
 for(k in 1:10)z[[paste0('sex_PC',k)]]<-z$sex*z[[paste0('PC',k)]]
 ml_ll<-numeric();ml_df<-integer()
 for(label in names(models)){
  f<-reformulate(c(covs,models[[label]],'gf','gm'),response='psoriasis_case')
  X<-model.matrix(f,z);y<-z$psoriasis_case;stopifnot(all(is.finite(X)),qr(X)$rank==ncol(X))
  ij<-match(c('gf','gm'),colnames(X))
  evaluate<-function(beta,firth){
   eta<-drop(X%*%beta);p<-plogis(eta);w<-p*(1-p);I<-crossprod(X,X*w)
   ll<-sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))
   adj<-0
   if(firth){ h<-rowSums((X%*%solve(I))*X)*w;adj<-h*(.5-p);ll<-ll+.5*as.numeric(determinant(I,logarithm=TRUE)$modulus) }
   list(ll=ll,score=drop(crossprod(X,y-p+adj)))
  }
  ml<-glm(f,data=z,family=binomial(),control=glm.control(maxit=150,epsilon=1e-10))
  stopifnot(ml$converged);ml_ll[label]<-as.numeric(logLik(ml));ml_df[label]<-length(coef(ml))
  for(method in c('ML','Firth')){
   key<-paste(t,label,method,sep='_');warn<-character();cat('START',key,fill=TRUE);flush.console()
   result<-withCallingHandlers(tryCatch({
    isf<-method=='Firth'
    fit<-if(isf)logistf(f,data=z,pl=FALSE,init=coef(ml),control=ctrl) else ml
    beta<-coef(fit);V<-vcov(fit);full<-evaluate(beta,isf)
    stopifnot(max(abs(full$score))<1e-6)
    if(isf)stopifnot(all(abs(fit$conv)<=c(1e-8,1e-6,1e-8)),abs(full$ll-fit$loglik[1])<1e-6)
    else stopifnot(abs(full$ll-as.numeric(logLik(fit)))<1e-6)
    nullfit<-function(idx){
     free<-setdiff(seq_len(ncol(X)),idx)
     if(isf){
      init<-beta;init[idx]<-0
      null<-logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=init,control=strict,modcontrol=logistf.mod.control(terms=free))
      b<-null$beta;check<-evaluate(b,TRUE)
      stopifnot(null$iter<500,all(abs(null$conv[c(1,3)])<=c(1e-8,1e-8)),abs(check$ll-null$loglik)<1e-6)
     }else{
      fn<-reformulate(setdiff(attr(terms(f),'term.labels'),colnames(X)[idx]),response='psoriasis_case')
      null<-glm(fn,data=z,family=binomial(),control=glm.control(maxit=150,epsilon=1e-10))
      stopifnot(null$converged);b<-setNames(rep(0,ncol(X)),colnames(X));b[names(coef(null))]<-coef(null);check<-evaluate(b,FALSE)
     }
     stopifnot(max(abs(b[idx]))<1e-12,max(abs(check$score[free]))<1e-6,check$ll<=full$ll+1e-6)
     list(p=pchisq(max(0,2*(full$ll-check$ll)),length(idx),lower.tail=FALSE),score=max(abs(check$score[free])))
    }
    nf<-nullfit(ij[1]);nj<-nullfit(ij)
    se<-sqrt(diag(V));reference<-previous[[paste(t,'01',method,sep='_')]]
    if(label=='A_original')stopifnot(abs(beta['gf']-reference$coef['g'])<1e-5,abs(log(nj$p)-log(reference$p))<1e-4)
    if(label=='D_all'){
     ff<-female[female$target==t & female$method==method,];stopifnot(nrow(ff)==1,abs(beta['gf']-ff$beta)<1e-5,abs(se['gf']-ff$se)<1e-5,abs(log(nf$p)-log(ff$p))<1e-4)
    }
    difference<-beta['gm']-beta['gf'];difference_se<-sqrt(V['gm','gm']+V['gf','gf']-2*V['gm','gf'])
    objects[[key]]<-list(coef=beta,var=V,ll=full$ll,columns=colnames(X))
    data.frame(target=t,model=label,method=method,female_OR=exp(beta['gf']),female_lower=exp(beta['gf']-qnorm(.975)*se['gf']),female_upper=exp(beta['gf']+qnorm(.975)*se['gf']),female_p=nf$p,male_OR=exp(beta['gm']),male_lower=exp(beta['gm']-qnorm(.975)*se['gm']),male_upper=exp(beta['gm']+qnorm(.975)*se['gm']),joint_p=nj$p,equal_native_effect_Wald_p=2*pnorm(-abs(difference/difference_se)),full_score=max(abs(full$score)),null_score=max(nf$score,nj$score),verified=TRUE,error='')
   },error=function(e)data.frame(target=t,model=label,method=method,verified=FALSE,error=conditionMessage(e))),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
   result$warnings<-paste(warn,collapse='; ');records[[key]]<-result
   saveRDS(list(records=records,models=objects),file.path(out,'model_progress.private.rds'))
   columns<-unique(unlist(lapply(records,names)));aligned<-lapply(records,function(r){for(n in setdiff(columns,names(r)))r[[n]]<-NA;r[,columns]})
   write.csv(do.call(rbind,aligned),file.path(out,'MODEL_RESULTS.csv'),row.names=FALSE)
   cat('DONE',key,'verified',result$verified,fill=TRUE);flush.console()
  }
 }
 for(pair in list(c('A_original','B_age2'),c('A_original','C_PCs'),c('B_age2','D_all'),c('C_PCs','D_all'),c('A_original','D_all'))){
  stat<-2*(ml_ll[pair[2]]-ml_ll[pair[1]]);df<-ml_df[pair[2]]-ml_df[pair[1]];stopifnot(stat>=-1e-6,df>0)
  nested[[paste(t,paste(pair,collapse='_'),sep='_')]]<-data.frame(target=t,reduced=pair[1],full=pair[2],df=df,ML_LR=max(stat,0),p=pchisq(max(stat,0),df,lower.tail=FALSE))
 }
 write.csv(do.call(rbind,nested),file.path(out,'NUISANCE_BLOCK_TESTS.csv'),row.names=FALSE)
}
final<-read.csv(file.path(out,'MODEL_RESULTS.csv'));stopifnot(nrow(final)==24,all(final$verified),all(is.na(final$warnings)))
writeLines('All 24 model fits and constrained likelihood tests verified; baseline and separate-female equivalence passed.',file.path(out,'MODEL_VERIFIED'))
