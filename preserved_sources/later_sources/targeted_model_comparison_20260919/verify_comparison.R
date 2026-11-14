args<-commandArgs(TRUE);o<-normalizePath(args[1]);v<-normalizePath(file.path(o,'../joint_G_sex_20260914/independent_validation_20260916'))
stopifnot(file.exists(file.path(o,'FITS_COMPLETE')));.libPaths(c(file.path(v,'Rlib'),.libPaths()));library(logistf)
d<-read.csv(file.path(v,'verified_irls_20260918/model_data.private.csv'));fits<-readRDS(file.path(o,'fits.private.rds'));results<-read.csv(file.path(o,'RESULTS.csv'));old<-readRDS(file.path(v,'verified_irls_20260918/target_fit_details.private.rds'))
covs<-c('age','age2',paste0('PC',1:10));f<-reformulate(c(covs,'g'),response='psoriasis_case');ctrl<-logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-10,gconv=1e-8,xconv=1e-10);allchecks<-list()
for(t in c('g0','g1','g2')){
 z<-d[d$sex==0 & !is.na(d[[t]]),];z$g<-z[[t]];X<-model.matrix(f,z);y<-z$psoriasis_case;j<-match('g',colnames(X));fit<-fits[[paste(t,'Firth',sep='_')]]
 stopifnot(identical(colnames(X),names(coef(fit))))
 evalfit<-function(beta){ eta<-drop(X%*%beta);p<-plogis(eta);w<-p*(1-p);info<-crossprod(X,X*w);h<-rowSums((X%*%solve(info))*X)*w;score<-drop(crossprod(X,y-p+h*(.5-p)));ll<-sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))+.5*as.numeric(determinant(info,logarithm=TRUE)$modulus);list(ll=ll,score=score) }
 full<-evalfit(coef(fit));stopifnot(abs(full$ll-fit$loglik[1])<1e-6,max(abs(full$score))<1e-6)
 init<-coef(fit);init[j]<-0;free<-setdiff(seq_len(ncol(X)),j)
 null<-logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=init,control=ctrl,modcontrol=logistf.mod.control(terms=free));ncheck<-evalfit(null$beta)
 stopifnot(abs(null$beta[j])<1e-12,abs(ncheck$ll-null$loglik)<1e-6,max(abs(ncheck$score[free]))<1e-6,null$iter<500,all(abs(null$conv[c(1,3)])<=c(1e-8,1e-8)),ncheck$ll<=full$ll+1e-6)
 p<-pchisq(2*(full$ll-ncheck$ll),1,lower.tail=FALSE);previous<-results[results$target==t & results$method=='Firth','p'];stopifnot(abs(log(p)-log(previous))<1e-4)
 prior<-old[[paste(t,'01','Firth',sep='_')]];stopifnot(is.null(prior$error),length(prior$warnings)==0)
 allchecks[[t]]<-data.frame(target=t,full_score=max(abs(full$score)),null_score=max(abs(ncheck$score[free])),recomputed_p=p,pooled_female_OR=exp(prior$coef['g']),separate_female_OR=exp(coef(fit)['g']),verified=TRUE)
 write.csv(do.call(rbind,allchecks),file.path(o,'INDEPENDENT_NUMERICAL_CHECK.csv'),row.names=FALSE);cat('VERIFIED',t,fill=TRUE);flush.console()
}
writeLines('All three full and constrained Firth adjusted scores and penalized likelihoods independently verified.',file.path(o,'NUMERICAL_VERIFIED'))
