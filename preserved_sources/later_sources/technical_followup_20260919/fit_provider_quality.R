args<-commandArgs(TRUE);stopifnot(length(args)==1);O<-normalizePath(args[1]);B<-dirname(O)
stopifnot(file.exists(file.path(O,'MODEL_INPUT_VERIFIED.json')),!file.exists(file.path(O,'QUALITY_MODELS_STARTED')))
writeLines(as.character(Sys.time()),file.path(O,'QUALITY_MODELS_STARTED'))
.libPaths(c(file.path(B,'joint_G_sex_20260914/independent_validation_20260916/Rlib'),.libPaths()));library(logistf)
d<-read.csv(file.path(O,'model_data.private.csv'));d<-d[d$sex==0,];covs<-c('age','age2',paste0('PC',1:10))
old<-read.csv(file.path(B,'targeted_model_comparison_20260919/RESULTS.csv'));ancestry_old<-read.csv(file.path(B,'deeper_analysis_20260919/ANCESTRY_RESULTS.csv'))
ctrl<-logistf.control(maxit=500,maxhs=50,fit='IRLS',lconv=1e-10,gconv=1e-8,xconv=1e-10)
rows<-list();objects<-list();plan<-list();gates<-list()
support<-function(z)all(table(factor(z$psoriasis_case,levels=0:1),factor(z$g>0,levels=c(FALSE,TRUE)))>=21)
add<-function(t,name,z,extra=character(),off=FALSE,baseline=FALSE){
 key<-paste(t,name,sep='_');plan[[key]]<<-list(target=t,name=name,data=z,extra=extra,offset=off,baseline=baseline)
}
for(t in c('g0','g1','g2')){
 z<-d[!is.na(d[[t]]),];z$g<-z[[t]]
 add(t,'base_no_offset',z,baseline=TRUE);add(t,'base_LOCO',z,off=TRUE)
 add(t,'ancestry_baseline',z,'ancestry4');add(t,'ancestry_site',z,c('ancestry4','site_cat'))
 c<-z[is.finite(z$mean_coverage),];c$coverage_z<-(c$mean_coverage-mean(c$mean_coverage))/sd(c$mean_coverage)
 add(t,'coverage_matched_base',c,c('ancestry4','site_cat'));add(t,'site_plus_coverage',c,c('ancestry4','site_cat','coverage_z'))
 q<-z[tolower(as.character(z[[paste0(t,'_quality_pass')]])) %in% c('true','1'),];add(t,'fixed_provider_quality',q,c('ancestry4','site_cat'))
 for(a in c('AFR','AMR')){
  s<-z[z$ancestry4==a,];tab<-table(s$site_cat,factor(s$psoriasis_case,levels=0:1));good<-rownames(tab)[apply(tab>=21,1,all)];s$site_cat<-ifelse(s$site_cat %in% good,s$site_cat,'OTHER')
  add(t,paste0(a,'_site'),s,'site_cat')
 }
 for(site in sort(setdiff(unique(z$site_cat),'OTHER'))){
  omit<-z[z$site_cat==site,];keep<-z[z$site_cat!=site,];eligible<-support(omit)&&support(keep)
  gates[[paste(t,site)]]<-list(eligible=eligible,excluded_cells=table(factor(omit$psoriasis_case,levels=0:1),factor(omit$g>0,levels=c(FALSE,TRUE))))
  if(eligible)add(t,paste0('leave_',site,'_out'),keep,c('ancestry4','site_cat'))
 }
}
saveRDS(list(plan=lapply(plan,function(z)list(target=z$target,name=z$name,ids=z$data$IID,extra=z$extra,offset=z$offset)),gates=gates),file.path(O,'QUALITY_MODEL_PLAN.private.rds'))
flush_results<-function(){cols<-unique(unlist(lapply(rows,names)));aligned<-lapply(rows,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[,cols]});write.csv(do.call(rbind,aligned),file.path(O,'QUALITY_SENSITIVITY_RESULTS.csv'),row.names=FALSE);saveRDS(objects,file.path(O,'quality_sensitivity_details.private.rds'))}
for(key in names(plan)[grepl('_fixed_provider_quality$',names(plan))]){
 spec<-plan[[key]];z<-spec$data;t<-spec$target;name<-spec$name
 if(!support(z)){rows[[key]]<-data.frame(target=t,model=name,method='both',status='support_gate_not_met');flush_results();next}
 extra<-spec$extra
 for(c in intersect(extra,c('ancestry4','site_cat'))){z[[c]]<-droplevels(factor(z[[c]]));if(nlevels(z[[c]])<2)extra<-setdiff(extra,c)}
 z$off<-if(spec$offset)z$loco else 0;stopifnot(all(is.finite(z$off)))
 f<-reformulate(c(covs,extra,'g','offset(off)'),response='psoriasis_case');X<-model.matrix(f,z);y<-z$psoriasis_case;o<-z$off
 if(qr(X)$rank!=ncol(X)){rows[[key]]<-data.frame(target=t,model=name,method='both',status='rank_deficient');flush_results();next}
 idx<-match('g',colnames(X));free<-setdiff(seq_len(ncol(X)),idx)
 evaluate<-function(b,firth){eta<-drop(X%*%b)+o;p<-plogis(eta);w<-p*(1-p);I<-crossprod(X,X*w);ll<-sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))));adj<-0;if(firth){h<-rowSums((X%*%solve(I))*X)*w;adj<-h*(.5-p);ll<-ll+.5*as.numeric(determinant(I,logarithm=TRUE)$modulus)};list(ll=ll,score=drop(crossprod(X,y-p+adj)))}
 for(method in c('ML','Firth')){
  k<-paste(key,method,sep='_');cat('START',k,fill=TRUE);flush.console();warn<-character()
  result<-withCallingHandlers(tryCatch({
   isf<-method=='Firth';fit<-if(isf)logistf(f,data=z,pl=FALSE,control=ctrl) else glm(f,data=z,family=binomial(),control=glm.control(maxit=200,epsilon=1e-10))
   b<-coef(fit);v<-vcov(fit);ev<-evaluate(b,isf);stopifnot(all(is.finite(b)),max(abs(ev$score))<1e-6)
   if(isf){stopifnot(all(abs(fit$conv)<=c(1e-8,1e-6,1e-8)),abs(ev$ll-fit$loglik[1])<1e-6);init<-b;init[idx]<-0;nf<-logistf:::logistf.fit(x=X,y=y,offset=o,firth=TRUE,init=init,control=ctrl,modcontrol=logistf.mod.control(terms=free));nb<-nf$beta;stopifnot(nf$iter<500,all(abs(nf$conv[c(1,3)])<=c(1e-8,1e-8)))}else{stopifnot(fit$converged,abs(ev$ll-as.numeric(logLik(fit)))<1e-6);nf<-glm(update(f,.~.-g),data=z,family=binomial(),control=glm.control(maxit=200,epsilon=1e-10));stopifnot(nf$converged);nb<-setNames(rep(0,ncol(X)),colnames(X));nb[names(coef(nf))]<-coef(nf)}
   en<-evaluate(nb,isf);stopifnot(abs(nb[idx])<1e-12,max(abs(en$score[free]))<1e-6,en$ll<=ev$ll+1e-6);if(isf)stopifnot(abs(en$ll-nf$loglik)<1e-6)
   se<-sqrt(v[idx,idx]);p<-pchisq(max(0,2*(ev$ll-en$ll)),1,lower.tail=FALSE)
   if(spec$baseline){r<-old[old$target==t & old$method==method,];stopifnot(nrow(r)==1,abs(b[idx]-r$beta)<1e-5,abs(se-r$se)<1e-5,abs(log(p)-log(r$p))<1e-4)}
   if(name=='ancestry_baseline'){r<-ancestry_old[ancestry_old$target==t & ancestry_old$method==method & ancestry_old$group=='ALL_ancestry_adjusted',];stopifnot(nrow(r)==1,abs(exp(b[idx])-r$OR)<1e-5,abs(log(p)-log(r$p))<1e-4)}
   objects[[k]]<-list(coef=b,var=v,columns=colnames(X))
   data.frame(target=t,model=name,method=method,status='verified',OR=exp(b[idx]),lower=exp(b[idx]-qnorm(.975)*se),upper=exp(b[idx]+qnorm(.975)*se),p=p,full_score=max(abs(ev$score)),null_score=max(abs(en$score[free])))
  },error=function(e)data.frame(target=t,model=name,method=method,status='failed',error=conditionMessage(e))),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
  result$warnings<-paste(warn,collapse='; ');rows[[k]]<-result;flush_results();cat('DONE',k,result$status,fill=TRUE);flush.console()
 }
}
r<-read.csv(file.path(O,'QUALITY_SENSITIVITY_RESULTS.csv'));writeLines(if(all(r$status %in% c('verified','support_gate_not_met')) && all(is.na(r$warnings)))'All estimable fits verified with no retained warnings.' else 'Completed with exceptions; inspect all results and warnings.',file.path(O,'QUALITY_MODELS_COMPLETED'))
