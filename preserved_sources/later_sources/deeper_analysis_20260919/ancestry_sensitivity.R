args<-commandArgs(TRUE);stopifnot(length(args)==1);out<-normalizePath(args[1])
stopifnot(file.exists(file.path(out,'ANCESTRY_INPUT_VERIFIED.json')))
if(file.exists(file.path(out,'ANCESTRY_STARTED')))stop('Already started')
writeLines('Female ancestry sensitivity: all case/carrier cells >=21 for each reported ancestry; fixed age, age2, PC1-PC10; ML and Firth; no LOCO; same-cohort exploratory analysis.',file.path(out,'ANCESTRY_STARTED'))
.libPaths(c(file.path(dirname(out),'joint_G_sex_20260914/independent_validation_20260916/Rlib'),.libPaths()));library(logistf)
d<-read.csv(file.path(out,'annotated_model_data.private.csv'));d<-d[d$sex==0,]
ctrl<-logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-10,gconv=1e-8,xconv=1e-10)
rows<-list();objects<-list();gates<-list();covs<-c('age','age2',paste0('PC',1:10))
for(t in c('g0','g1','g2'))for(group in c('ALL_ancestry_adjusted','AFR','AMR','EUR')){
 z<-d[!is.na(d[[t]]),];if(group!='ALL_ancestry_adjusted')z<-z[z$ancestry4==group,]
 z$g<-z[[t]];cells<-table(factor(z$psoriasis_case,levels=0:1),factor(z$g>0,levels=c(FALSE,TRUE)))
 gate<-all(cells>=21);gates[[paste(t,group)]]<-list(cells=cells,passed=gate)
 if(!gate){rows[[paste(t,group)]]<-data.frame(target=t,group=group,method='both',status='support_gate_not_met');next}
 z$ancestry4<-droplevels(factor(z$ancestry4))
 f<-reformulate(c(covs,if(group=='ALL_ancestry_adjusted' && nlevels(z$ancestry4)>1)'ancestry4','g'),response='psoriasis_case')
 X<-model.matrix(f,z);y<-z$psoriasis_case
 if(qr(X)$rank!=ncol(X)){rows[[paste(t,group)]]<-data.frame(target=t,group=group,method='both',status='rank_deficient');next}
 idx<-match('g',colnames(X));free<-setdiff(seq_len(ncol(X)),idx)
 evaluate<-function(b,firth){eta<-drop(X%*%b);p<-plogis(eta);w<-p*(1-p);I<-crossprod(X,X*w);ll<-sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))));a<-0;if(firth){h<-rowSums((X%*%solve(I))*X)*w;a<-h*(.5-p);ll<-ll+.5*as.numeric(determinant(I,logarithm=TRUE)$modulus)};list(ll=ll,score=drop(crossprod(X,y-p+a)))}
 ml<-glm(f,data=z,family=binomial(),control=glm.control(maxit=200,epsilon=1e-10));stopifnot(ml$converged)
 for(method in c('ML','Firth')){
  key<-paste(t,group,method,sep='_');cat('START',key,fill=TRUE);flush.console();warnings<-character()
  r<-withCallingHandlers(tryCatch({
   firth<-method=='Firth';fit<-if(firth)logistf(f,data=z,pl=FALSE,init=coef(ml),control=ctrl) else ml
   b<-coef(fit);v<-vcov(fit);ev<-evaluate(b,firth);stopifnot(max(abs(ev$score))<1e-6)
   if(firth){stopifnot(all(abs(fit$conv)<=c(1e-8,1e-6,1e-8)),abs(ev$ll-fit$loglik[1])<1e-6);init<-b;init[idx]<-0;nf<-logistf:::logistf.fit(x=X,y=y,firth=TRUE,init=init,control=ctrl,modcontrol=logistf.mod.control(terms=free));nb<-nf$beta;stopifnot(nf$iter<500,all(abs(nf$conv[c(1,3)])<=c(1e-8,1e-8)))}else{nf<-glm(update(f,.~.-g),data=z,family=binomial(),control=glm.control(maxit=200,epsilon=1e-10));stopifnot(nf$converged,abs(ev$ll-as.numeric(logLik(fit)))<1e-6);nb<-setNames(rep(0,ncol(X)),colnames(X));nb[names(coef(nf))]<-coef(nf)}
   en<-evaluate(nb,firth);stopifnot(max(abs(en$score[free]))<1e-6,abs(nb[idx])<1e-12,en$ll<=ev$ll+1e-6)
   if(firth)stopifnot(abs(en$ll-nf$loglik)<1e-6)
   se<-sqrt(v[idx,idx]);objects[[key]]<-list(coef=b,var=v)
   data.frame(target=t,group=group,method=method,status='verified',OR=exp(b[idx]),lower=exp(b[idx]-qnorm(.975)*se),upper=exp(b[idx]+qnorm(.975)*se),p=pchisq(max(0,2*(ev$ll-en$ll)),1,lower.tail=FALSE),full_score=max(abs(ev$score)),null_score=max(abs(en$score[free])))
  },error=function(e)data.frame(target=t,group=group,method=method,status='failed',error=conditionMessage(e))),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  r$warnings<-paste(warnings,collapse='; ');rows[[key]]<-r
  cols<-unique(unlist(lapply(rows,names)));aligned<-lapply(rows,function(a){for(n in setdiff(cols,names(a)))a[[n]]<-NA;a[,cols]});write.csv(do.call(rbind,aligned),file.path(out,'ANCESTRY_RESULTS.csv'),row.names=FALSE)
  saveRDS(list(gates=gates,objects=objects),file.path(out,'ancestry_details.private.rds'));cat('DONE',key,r$status,fill=TRUE);flush.console()
 }
}
cols<-unique(unlist(lapply(rows,names)));aligned<-lapply(rows,function(a){for(n in setdiff(cols,names(a)))a[[n]]<-NA;a[,cols]});r<-do.call(rbind,aligned);write.csv(r,file.path(out,'ANCESTRY_RESULTS.csv'),row.names=FALSE)
saveRDS(list(gates=gates,objects=objects),file.path(out,'ancestry_details.private.rds'))
stopifnot(!any(r$status=='failed'),all(is.na(r$warnings)|r$warnings==''))
writeLines('All estimable ancestry sensitivities numerically verified; unsupported strata retained as unavailable.',file.path(out,'ANCESTRY_VERIFIED'))
