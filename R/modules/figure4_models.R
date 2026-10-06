# Figure 4 / S8-S9: calculations from the common input data.
# Uses the shared data-building and model helpers in R/core.
figure4_one <- function(d,...){z<-list(...);for(k in names(z))d<-d[d[[k]]==z[[k]],,drop=FALSE];if(nrow(d)!=1L)stop('Expected one Figure4 row: ',paste(z,collapse=' / '));d}
figure4_plain <- function(d){d<-as.data.frame(d);for(k in names(d))if(length(dim(d[[k]]))==1L)d[[k]]<-as.vector(d[[k]]);d}
figure4_write_workbook <- function(tables,out) {
  tables<-lapply(tables,figure4_plain)
  # Full tables retain descriptive CSV filenames; worksheets use <=31 characters.
  sheet<-setNames(names(tables),names(tables))
  short<-c(S9_CTAG_component_decomposition='S9_components',S9_clinical_normal_reference_sensitivity='S9_clinical_normal',
    S9_crossfitted_clinical_sensitivity='S9_crossfit_clinical',S9_primary_candidate_features='S9_candidates',
    S9_nested_selection_frequency='S9_nested_selection',S9_candidate_OOF_predictions='S9_candidate_OOF',
    S9_fixed_OOF_predictions='S9_fixed_OOF',S9_nested_OOF_predictions='S9_nested_OOF')
  sheet[names(short)]<-short
  stopifnot(all(nchar(sheet)<=31L),!anyDuplicated(sheet))
  wb<-openxlsx::createWorkbook(creator='KTB140 R pipeline')
  hdr<-openxlsx::createStyle(fontName='Arial',fontSize=10,textDecoration='bold',fgFill='#E8F0F0',wrapText=TRUE)
  for(n in names(tables)){
    write_csv(tables[[n]],file.path(out,'results',paste0(n,'.csv')))
    sn<-sheet[[n]];openxlsx::addWorksheet(wb,sn);openxlsx::writeData(wb,sn,tables[[n]],headerStyle=hdr,keepNA=FALSE)
    openxlsx::freezePane(wb,sn,firstRow=TRUE);openxlsx::setColWidths(wb,sn,seq_along(tables[[n]]),19)
  }
  file<-file.path(out,'results/Figure4_S8-S9_R_results.xlsx');openxlsx::saveWorkbook(wb,file,overwrite=FALSE)
  checks<-list()
  for(n in names(tables)){
    a<-tables[[n]];b<-openxlsx::read.xlsx(file,sheet=sheet[[n]],check.names=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE)
    if(!identical(names(a),names(b))||nrow(a)!=nrow(b))stop('Workbook schema mismatch: ',n)
    for(k in names(a)){
      x<-a[[k]];y<-b[[k]];if(is.factor(x))x<-as.character(x)
      if(!identical(is.na(x),is.na(y)))stop('Workbook missingness mismatch: ',n,'/',k)
      ok<-!is.na(x)
      eq<-if(is.numeric(x))all(abs(as.numeric(y[ok])-x[ok])<=1e-12+1e-12*abs(x[ok])) else
        if(is.logical(x))identical(flag(x),flag(y)) else
        if(is.character(x)&&is.logical(y)&&all(tolower(x[ok])%in%c('true','false')))identical(tolower(x[ok]),tolower(as.character(y[ok]))) else identical(as.character(x[ok]),as.character(y[ok]))
      if(!eq)stop('Workbook contents mismatch: ',n,'/',k)
    }
    checks[[n]]<-data.frame(table=n,worksheet=sheet[[n]],rows=nrow(a),columns=ncol(a),pass=TRUE)
  }
  write_csv(do.call(rbind,checks),file.path(out,'checks/workbook_readback.csv'))
  write_csv(data.frame(table=names(sheet),worksheet=unname(sheet)),file.path(out,'results/WORKBOOK_SHEET_MAP.csv'))
  invisible(tables)
}
figure4_read_plans <- function(root,ids) {
  files<-c(comparison='CV_comparison_folds.csv',fixed='CV_fixed_folds.csv',outer='CV_nested_outer_folds.csv',inner='CV_nested_inner_folds.csv')
  p<-lapply(files,function(n)read_csv(file.path(root,figure4_config$folds_dir,n)))
  reps<-c(comparison=50L,fixed=100L,outer=50L)
  for(role in names(reps)){
    q<-p[[role]];assert_cols(q,c('Barcode','repeat_id','fold'))
    stopifnot(setequal(q$Barcode,ids),!anyDuplicated(q[,c('Barcode','repeat_id')]),
      identical(sort(unique(q$repeat_id)),seq_len(reps[[role]])),all(q$fold%in%1:10),nrow(q)==length(ids)*reps[[role]])
    for(r in seq_len(reps[[role]])){
      z<-q[q$repeat_id==r,];stopifnot(setequal(z$Barcode,ids),setequal(z$fold,1:10))
    }
  }
  q<-p$inner;assert_cols(q,c('Barcode','repeat_id','fold','outer_repeat','outer_fold'))
  stopifnot(all(q$repeat_id==1L),all(q$fold%in%1:5),!anyDuplicated(q[,c('Barcode','outer_repeat','outer_fold')]))
  for(r in 1:50)for(f in 1:10){
    train<-p$outer$Barcode[p$outer$repeat_id==r&p$outer$fold!=f]
    z<-q[q$outer_repeat==r&q$outer_fold==f,];stopifnot(setequal(z$Barcode,train),length(train)==nrow(z),setequal(z$fold,1:5))
  }
  p
}
# Training-only FWL nuisance adjustment. Pseudocounts, state definition and initial
# feature nomination are cohort-defined, not claimed as fully nested discovery.
figure4_predict_split <- function(train,test,features) {
  cv<-figure4_config$serum_covariates
  stopifnot(length(intersect(train$Barcode,test$Barcode))==0L,!anyDuplicated(train$Barcode),!anyDuplicated(test$Barcode))
  ym<-mean(train$.target);ys<-sd(train$.target);if(!is.finite(ys)||ys<=0)stop('Invalid training-target scale')
  yt<-(train$.target-ym)/ys;ye<-(test$.target-ym)/ys
  Ctr<-cbind(1,as.matrix(train[,cv,drop=FALSE]));Cte<-cbind(1,as.matrix(test[,cv,drop=FALSE]))
  cy<-stats::lm.fit(Ctr,yt);if(cy$rank<ncol(Ctr))stop('Rank-deficient training nuisance fit')
  yr<-yt-as.numeric(Ctr%*%cy$coefficients);ey<-ye-as.numeric(Cte%*%cy$coefficients)
  Xt<-as.matrix(train[,features,drop=FALSE]);Xe<-as.matrix(test[,features,drop=FALSE]);storage.mode(Xt)<-'double';storage.mode(Xe)<-'double'
  mu<-colMeans(Xt);sc<-apply(Xt,2,sd);sc[!is.finite(sc)|sc==0]<-1
  Xt<-sweep(sweep(Xt,2,mu,'-'),2,sc,'/');Xe<-sweep(sweep(Xe,2,mu,'-'),2,sc,'/')
  cx<-stats::lm.fit(Ctr,Xt);XR<-Xt-Ctr%*%cx$coefficients;ER<-Xe-Cte%*%cx$coefficients
  fit<-stats::lm.fit(XR,yr);if(fit$rank<ncol(XR))stop('Rank-deficient residual predictor')
  b<-as.numeric(fit$coefficients);pred<-as.numeric(ER%*%b)
  pr<-data.frame(Barcode=test$Barcode,target=ey,pred=pred,serum_only_score=as.numeric(Xe%*%b),
    baseline_error=ey^2,squared_error=(ey-pred)^2,stringsAsFactors=FALSE)
  stopifnot(all(is.finite(as.matrix(pr[,-1]))))
  # Independent full-covariate formulation: residual prediction = full prediction
  # minus the covariate-only fitted contribution, all fitted on training donors.
  full<-stats::lm.fit(cbind(Ctr,Xt),yt)
  expected<-as.numeric(cbind(Cte,Xe)%*%full$coefficients-Cte%*%cy$coefficients)
  delta<-max(abs(pred-expected));if(delta>config$deterministic_atol+config$deterministic_rtol*max(abs(expected)))stop('FWL/full-covariate identity failed')
  audit<-data.frame(n_train=nrow(train),n_test=nrow(test),n_features=length(features),nuisance_rank=cy$rank,predictor_rank=fit$rank,
    target_mean=ym,target_sd=ys,FWL_max_difference=delta,train_test_overlap=0L,stringsAsFactors=FALSE)
  list(predictions=pr,audit=audit,parameters=data.frame(feature=features,mean=mu,sd=sc,beta=b,row.names=NULL))
}
figure4_cv_summary <- function(P,features) {
  avg<-stats::aggregate(P[,c('target','pred','serum_only_score')],list(Barcode=P$Barcode),mean)
  rm<-tapply(P$squared_error,P$repeat_id,mean);r<-stats::cor(avg$target,avg$pred)
  data.frame(subset=paste(features,collapse=' + '),n_features=length(features),mean_cv_mse=mean(P$squared_error),
    repeat_mse_sd=sd(rm),repeat_mse_q025=unname(quantile(rm,.025,type=7)),repeat_mse_q975=unname(quantile(rm,.975,type=7)),
    cv_r=r,cv_R2_vs_covariates=1-sum(P$squared_error)/sum(P$baseline_error),donor_averaged_r_squared=r*r,n=nrow(avg),
    note='Restricted-candidate internal evaluation. Repeat MSE percentiles describe split variability, not an independent-sample CI.',stringsAsFactors=FALSE)
}
figure4_cv <- function(d,subsets,folds,role) {
  all<-list();perf<-list();aud<-list();aa<-0L
  for(i in seq_along(subsets)){
    features<-subsets[[i]];label<-paste(features,collapse=' + ');pr<-list();k<-0L
    for(r in sort(unique(folds$repeat_id)))for(f in 1:10){
      ids<-folds$Barcode[folds$repeat_id==r&folds$fold==f];test<-d[d$Barcode%in%ids,,drop=FALSE];train<-d[!d$Barcode%in%ids,,drop=FALSE]
      z<-figure4_predict_split(train,test,features);q<-z$predictions;q$repeat_id<-r;q$fold<-f;q$subset<-label;k<-k+1L;pr[[k]]<-q
      aa<-aa+1L;z$audit$role<-role;z$audit$repeat_id<-r;z$audit$fold<-f;z$audit$subset<-label;aud[[aa]]<-z$audit
    }
    P<-do.call(rbind,pr);perf[[i]]<-figure4_cv_summary(P,features);all[[label]]<-P
    cat('Completed',role,label,'with',length(unique(folds$repeat_id)),'repeats\n')
  }
  perf<-do.call(rbind,perf);perf<-perf[order(perf$mean_cv_mse,perf$n_features,perf$subset),,drop=FALSE]
  list(performance=perf,predictions=do.call(rbind,all),winner=perf$subset[1],audit=do.call(rbind,aud))
}
figure4_nested <- function(d,subsets,plans) {
  predictions<-list();sel<-list();audit<-list();errors<-list();k<-0L;ia<-0L
  for(r in 1:50)for(f in 1:10){
    k<-k+1L;ids<-plans$outer$Barcode[plans$outer$repeat_id==r&plans$outer$fold==f]
    tr<-d[!d$Barcode%in%ids,,drop=FALSE];te<-d[d$Barcode%in%ids,,drop=FALSE]
    inner<-plans$inner[plans$inner$outer_repeat==r&plans$inner$outer_fold==f,,drop=FALSE]
    stopifnot(setequal(inner$Barcode,tr$Barcode))
    mse<-numeric(length(subsets))
    for(s in seq_along(subsets)){
      ee<-numeric()
      for(j in 1:5){
        va<-inner$Barcode[inner$fold==j];q<-figure4_predict_split(tr[!tr$Barcode%in%va,,drop=FALSE],tr[tr$Barcode%in%va,,drop=FALSE],subsets[[s]])
        ee<-c(ee,q$predictions$squared_error);ia<-ia+1L;q$audit$role<-'nested_inner';q$audit$repeat_id<-r;q$audit$fold<-f;q$audit$inner_fold<-j;q$audit$subset<-paste(subsets[[s]],collapse=' + ');audit[[ia]]<-q$audit
      }
      mse[s]<-mean(ee)
    }
    # Same deterministic tie rule as the source: fewer predictors; then fixed subset order.
    pick<-order(mse,lengths(subsets),seq_along(subsets))[1];features<-subsets[[pick]]
    q<-figure4_predict_split(tr,te,features);P<-q$predictions;P$repeat_id<-r;P$fold<-f;P$subset<-paste(features,collapse=' + ');predictions[[k]]<-P
    ia<-ia+1L;q$audit$role<-'nested_outer';q$audit$repeat_id<-r;q$audit$fold<-f;q$audit$inner_fold<-NA_integer_;q$audit$subset<-P$subset[1];audit[[ia]]<-q$audit
    errors[[k]]<-data.frame(repeat_id=r,fold=f,subset=vapply(subsets,paste,character(1),collapse=' + '),inner_mse=mse,selected=seq_along(subsets)==pick)
    sel[[k]]<-data.frame(repeat_id=r,fold=f,subset=P$subset[1])
    if(f==10)cat('Nested selection completed outer repeat',r,'of 50\n')
  }
  P<-do.call(rbind,predictions);selection<-do.call(rbind,sel);avg<-aggregate(P[,c('target','pred')],list(Barcode=P$Barcode),mean)
  metrics<-data.frame(cv_r=cor(avg$target,avg$pred),cv_R2_vs_covariates=1-sum(P$squared_error)/sum(P$baseline_error),n=nrow(avg),
    scope='Subset selection nested within the fixed nominated candidate family; initial screen and tissue construction are not nested.')
  count<-table(selection$subset);freq<-data.frame(subset=names(count),outer_folds=as.integer(count),frequency=as.numeric(count)/nrow(selection))
  list(predictions=P,metrics=metrics,frequency=freq,selection=selection,audit=dplyr::bind_rows(audit),inner_errors=do.call(rbind,errors))
}
figure4_lobule_screen <- function(dat,record) {
  mis<-dat$mispa$log;lob<-dat$lobules[flag(dat$lobules$primary_serum),,drop=FALSE];ep<-list()
  for(y in c('p63_fraction','er_fraction','erp63_enrichment_called')){
    d<-if(y=='erp63_enrichment_called')lob[flag(lob$eligible_erp63_100),,drop=FALSE] else lob
    nm<-switch(y,p63_fraction='p63_fraction',er_fraction='ER_fraction',erp63_enrichment_called='ERp63_enrichment')
    ep[[nm]]<-list(data=d,y=y,weight=NULL,interaction=FALSE)
  }
  for(sub in c('all_lobules','count_rich'))for(n in c('p63_high','ER_low','ERp63_high','composite')){
    d<-dat$burdens[[sub]];d<-d[flag(d$primary_serum),,drop=FALSE]
    ep[[paste(sub,n,sep='_')]]<-list(data=d,y=paste0(n,'_burden'),weight=paste0(n,'_eligible_n'),interaction=FALSE)
  }
  d<-lob[flag(lob$eligible_erp63_100),,drop=FALSE]
  d$mean_enrichment<-ave(d$erp63_enrichment_called,d$Barcode,FUN=mean);d$within_enrichment<-d$erp63_enrichment_called-d$mean_enrichment
  for(n in c('er','p63'))ep[[paste0('coupling_',n)]]<-list(data=d,y=paste0(n,'_fraction'),weight=NULL,interaction=TRUE)
  stopifnot(identical(names(ep),figure4_config$endpoints));rows<-list();membership<-list();scales<-list()
  for(n in names(ep)){
    o<-ep[[n]];d<-o$data;ids<-sort(unique(d$Barcode));idx<-match(ids,mis$Barcode);stopifnot(!anyNA(idx))
    # Preserve the original standardization frame: one value per donor in the
    # endpoint input table, before outcome/weight complete-case fitting. This is
    # 119 for the burden input tables, even when a spatial burden fits 113 donors.
    # Repeated lobule rows never increase a donor's weight in serum standardization.
    for(f in dat$mispa$features){
      v<-mis[[f]][idx];d$.feature<-stdz(v)[match(d$Barcode,ids)]
      terms<-c('.feature','Age','BMI');term<-'.feature'
      if(o$interaction){d$.interaction<-d$.feature*d$within_enrichment;terms<-c('.feature','mean_enrichment','within_enrichment','.interaction','Age','BMI');term<-'.interaction'}
      z<-robust_fit(d,o$y,terms,term,cluster=is.null(o$weight),weight=o$weight);r<-z$row;r$endpoint<-n;r$feature<-f
      r$p_cluster_t<-if(is.null(o$weight))2*pt(-abs(r$beta/r$se),r$donors-1) else NA_real_
      rows[[length(rows)+1L]]<-r;record(paste('S10',n,f,sep='|'),z,save=f%in%figure4_config$targeted_features)
      scales[[length(scales)+1L]]<-data.frame(endpoint=n,feature=f,mean=mean(v),sd=sd(v),standardization_donors=length(ids),model_donors=r$donors)
      if(f==dat$mispa$features[1])membership[[n]]<-data.frame(endpoint=n,Barcode=z$data$Barcode,
        lobule_uid=if(is.null(o$weight))z$data$lobule_uid else NA_character_,weight=if(is.null(o$weight))1 else z$data$.w)
    }
    cat('Completed exploratory serum endpoint',n,'\n')
  }
  screen<-adjust_family(do.call(rbind,rows),'endpoint','BH',name='q_BH209');screen<-adjust_family(screen,method='BH',name='q_BH2717')
  list(screen=screen,membership=do.call(rbind,membership),scales=do.call(rbind,scales))
}
recompute_figure4 <- function(root,out) {
  for(n in c('checks','results','plot_data','models','resampling'))dir.create(file.path(out,n),recursive=TRUE,showWarnings=FALSE)
  dat<-read_common(root);dsp<-build_dsp(root,dat$clinical);states<-build_states(dat,dsp)
  write_csv(dsp$donor,file.path(out,'checks/rebuilt_DSP_donors.csv'));write_csv(states,file.path(out,'checks/rebuilt_Tissue_states.csv'))
  d<-states[flag(states$serum_common),,drop=FALSE];mis<-dat$mispa$log[match(d$Barcode,dat$mispa$log$Barcode),,drop=FALSE]
  stopifnot(nrow(d)==109L,!anyNA(mis$Barcode),identical(d$Barcode,mis$Barcode),length(dat$mispa$features)==209L)
  fits<-list();coeff<-list();register<-list()
  record<-function(id,obj,save=FALSE){
    m<-obj$fit;cf<-coef(m);coeff[[length(coeff)+1L]]<<-data.frame(model=id,term=names(cf),beta=unname(cf),se=sqrt(diag(obj$V)),n=nrow(obj$data),df=df.residual(m))
    register[[length(register)+1L]]<<-data.frame(model=id,formula=paste(deparse(formula(m)),collapse=' '),n=nrow(obj$data),donors=length(unique(obj$data$Barcode)),rank=m$rank,df=df.residual(m),covariance=obj$row$covariance,p_reference=obj$row$p_reference,ci_reference=obj$row$ci_reference)
    if(save)fits[[id]]<<-list(fit=m,V=obj$V,inference=obj$row)
  }
  rows<-list();spec<-list(primary=figure4_config$serum_covariates,sensitivity=figure4_config$sensitivity_covariates)
  for(a in names(spec))for(s in figure4_config$states){
    q<-d;q$.state<-stdz(q[[s]])
    for(f in dat$mispa$features){
      q$.feature<-stdz(mis[[f]]);z<-robust_fit(q,'.state',c('.feature',spec[[a]]),'.feature',serum=TRUE)
      r<-z$row;r$state<-s;r$feature<-f;r$adjustment<-a;rows[[length(rows)+1L]]<-r
      record(paste('S8',a,s,f,sep='|'),z,save=f%in%figure4_config$targeted_features)
    }
  }
  screen<-adjust_family(do.call(rbind,rows),c('state','adjustment'),'BH',name='q_BH209')
  screen<-adjust_family(screen,'adjustment','BH',name='q_BH836')
  hits<-aggregate(as.integer(screen$q_BH209<.05),screen[,c('state','adjustment')],sum);names(hits)[3]<-'hits'
  candidates<-sort(screen$feature[screen$state=='myeloid'&screen$adjustment=='primary'&screen$q_BH209<.05])
  # Candidates are taken from the screen results; they are not an input to fitting.
  if(!length(candidates)||length(candidates)>4L)stop('Candidate set outside the bounded analysis scope; inspect the screen before prediction.')
  components<-list()
  for(a in names(spec))for(s in c('DSP_block','RNA_Myeloid')){
    q<-d;q$.state<-stdz(q[[s]]);q$.feature<-stdz(mis[['CTAG1B.1']]);z<-robust_fit(q,'.state',c('.feature',spec[[a]]),'.feature',serum=TRUE)
    r<-z$row;r$component<-s;r$adjustment<-a;components[[length(components)+1L]]<-r;record(paste('component',a,s,sep='|'),z,TRUE)
  }
  components<-adjust_family(do.call(rbind,components),'adjustment')
  # Save the screen before the longer cross-validation calculation.
  for(t in c('screen','components','hits'))write_csv(get(t),file.path(out,'checks',paste0('completed_',t,'.csv')))
  plans<-figure4_read_plans(root,d$Barcode)
  for(role in names(plans))write_csv(plans[[role]],file.path(out,'resampling',paste0('CV_',role,'_folds.csv')))
  preddata<-d[,c('Barcode','Age','BMI','ParousFlag','TSLB_years','design_postpartum_serum'),drop=FALSE];preddata$.target<-d$myeloid
  for(f in unique(c(candidates,'CTAG1B.1')))preddata[[f]]<-mis[[f]]
  subsets<-unlist(lapply(seq_along(candidates),function(k)combn(candidates,k,simplify=FALSE)),recursive=FALSE)
  cmp<-figure4_cv(preddata,subsets,plans$comparison,'candidate');selected<-strsplit(cmp$winner,' + ',fixed=TRUE)[[1]]
  fixedsets<-list(selected);if(!identical(selected,'CTAG1B.1'))fixedsets[[2]]<-'CTAG1B.1'
  fixed<-figure4_cv(preddata,fixedsets,plans$fixed,'fixed');nested<-figure4_nested(preddata,subsets,plans)
  cv_audit<-dplyr::bind_rows(cmp$audit,fixed$audit,nested$audit);write_csv(cv_audit,file.path(out,'checks/CV_training_audit.csv'))
  write_csv(nested$inner_errors,file.path(out,'results/S9_nested_inner_errors.csv'));write_csv(nested$selection,file.path(out,'results/S9_nested_outer_selection.csv'))
  # Global locked serum-only component; covariate coefficients are never included
  # in the score plotted against parity or TSLB.
  X<-as.matrix(preddata[,selected,drop=FALSE]);mu<-colMeans(X);sc<-apply(X,2,sd);if(any(!is.finite(sc)|sc<=0))stop('Invalid lock feature SD')
  Z<-sweep(sweep(X,2,mu,'-'),2,sc,'/');colnames(Z)<-paste0('feature_',seq_along(selected))
  lockdata<-cbind(preddata,Z);lockdata$.target_z<-stdz(lockdata$.target)
  lm<-stats::lm(reformulate(c(figure4_config$serum_covariates,colnames(Z)),response='.target_z'),lockdata)
  if(lm$rank!=length(coef(lm)))stop('Rank-deficient full-cohort lock')
  b<-coef(lm)[colnames(Z)];fits$locked_full_covariate_model<-lm
  scores<-join_one(dat$clinical,dat$mispa$log);A<-as.matrix(scores[,selected,drop=FALSE]);scores$locked_score<-as.numeric(sweep(sweep(A,2,mu,'-'),2,sc,'/')%*%b)
  pars<-data.frame(feature=selected,log2_mean=mu,log2_sd=sc,beta=as.numeric(b),pseudocount=dat$mispa$pseudocounts[selected],training_n=nrow(d),selected_subset=cmp$winner,
    formula='sum(beta * (log2(seroreactivity + pseudocount)-log2_mean)/log2_sd); no clinical terms',row.names=NULL)
  clin<-list();clinfit<-list()
  for(role in c('parity','postpartum')){
    q<-if(role=='parity')scores[flag(scores$primary_serum)&is.finite(scores$ParousFlag),,drop=FALSE] else scores[flag(scores$design_postpartum_serum),,drop=FALSE]
    terms<-c(if(role=='parity')'ParousFlag' else 'TSLB_years','Age','BMI')
    z<-robust_fit(q,'locked_score',terms,terms[1]);r<-z$row;r$model<-role;r$units<-'serum-only score';clin[[length(clin)+1L]]<-r;clinfit[[role]]<-z;record(paste0('clinical_',role),z,TRUE)
    q$.measured_CTAG<-stdz(q[['CTAG1B.1']]);z<-robust_fit(q,'.measured_CTAG',terms,terms[1]);r<-z$row;r$model<-paste0('measured_CTAG_',role);r$units<-'within-subset SD';clin[[length(clin)+1L]]<-r;record(r$model,z,TRUE)
  }
  clin<-do.call(rbind,clin);stopifnot(figure4_one(clin,model='parity')$n==131L,figure4_one(clin,model='postpartum')$n==43L)
  norm<-clin;norm$lower<-norm$beta-qnorm(.975)*norm$se;norm$upper<-norm$beta+qnorm(.975)*norm$se;norm$p_raw<-2*pnorm(-abs(norm$beta/norm$se));norm$p_reference<-'normal sensitivity';norm$ci_reference<-'normal'
  P<-fixed$predictions[fixed$predictions$subset==cmp$winner,,drop=FALSE]
  oo<-aggregate(P[,c('target','pred','serum_only_score')],list(Barcode=P$Barcode),mean);oo<-join_one(oo,dat$clinical)
  cross<-list()
  for(role in c('parity','postpartum')){
    q<-if(role=='parity')oo else oo[flag(oo$design_postpartum_serum),,drop=FALSE]
    terms<-c(if(role=='parity')'ParousFlag' else 'TSLB_years','Age','BMI');z<-robust_fit(q,'serum_only_score',terms,terms[1]);r<-z$row;r$model<-role;cross[[role]]<-r;record(paste0('crossfit_clinical_',role),z,TRUE)
  };cross<-do.call(rbind,cross);stopifnot(figure4_one(cross,model='postpartum')$n==40L)
  # In-sample nominated-feature proxy; both sides residualized for the same covariates.
  pr<-d;pr$proxy<-rowMeans(do.call(cbind,lapply(candidates,function(f)stdz(mis[[f]]))))
  pr$proxy_residual<-residualize(pr,'proxy',figure4_config$serum_covariates);pr$tissue_residual<-residualize(pr,'myeloid',figure4_config$serum_covariates)
  proxy<-data.frame(n=nrow(pr),candidate_features=paste(candidates,collapse=' + '),adjusted_in_sample_r=cor(pr$proxy_residual,pr$tissue_residual))
  prfit<-robust_fit(pr,'tissue_residual','proxy_residual','proxy_residual');record('descriptive_proxy_line',prfit,TRUE)
  predline<-robust_fit(oo,'pred','target','target');record('descriptive_OOF_line',predline,TRUE)
  pp<-states[flag(states$design_postpartum_serum),,drop=FALSE];pp$.myeloid_z<-stdz(pp$myeloid)
  z<-robust_fit(pp,'.myeloid_z',c('TSLB_years','Age','BMI'),'TSLB_years');tissue<-z$row;tissue$state<-'myeloid';record('tissue_state_TSLB43',z,TRUE)
  ls<-figure4_lobule_screen(dat,record)
  lss<-do.call(rbind,lapply(figure4_config$endpoints,function(ep){q<-ls$screen[ls$screen$endpoint==ep,];data.frame(endpoint=ep,nominal=sum(q$p_raw<.01),BH=sum(q$q_BH209<.05),overall_BH=sum(q$q_BH2717<.05),donors=unique(q$donors),best_feature=q$feature[which.min(q$p_raw)])}))
  curves<-prediction_grid(clinfit$postpartum,'TSLB_years',seq(2,10,length.out=figure4_config$curve_points));names(curves)[1]<-'TSLB_years'
  pg<-prediction_grid(prfit,'proxy_residual',seq(min(pr$proxy_residual),max(pr$proxy_residual),length.out=161))
  cg<-prediction_grid(predline,'target',seq(min(oo$target),max(oo$target),length.out=161))
  top<-unique(unlist(lapply(figure4_config$states,function(s){q<-screen[screen$state==s&screen$adjustment=='primary',];q$feature[order(q$p_raw,match(q$feature,dat$mispa$features))][1:3]})))
  volc<-screen[screen$state=='myeloid'&screen$adjustment=='primary',];volc$minuslog10p<- -log10(volc$p_raw)
  membership<-rbind(data.frame(analysis='common_state',Barcode=d$Barcode),data.frame(analysis='serum_parity',Barcode=clinfit$parity$data$Barcode),data.frame(analysis='postpartum_serum',Barcode=clinfit$postpartum$data$Barcode),data.frame(analysis='crossfit_postpartum',Barcode=oo$Barcode[flag(oo$design_postpartum_serum)]))
  transformations<-data.frame(feature=dat$mispa$features,pseudocount=dat$mispa$pseudocounts[dat$mispa$features])
  logs<-dat$mispa$log
  write_csv(ls$membership,file.path(out,'checks/S10_membership.csv'));write_csv(ls$scales,file.path(out,'checks/S10_feature_scales.csv'))
  write_csv(dat$donors,file.path(out,'checks/rebuilt_MxIF_donors.csv'));write_csv(logs,file.path(out,'checks/MISPA_log2.csv'))
  write_csv(cv_audit,file.path(out,'checks/CV_training_audit.csv'))
  tables<-list(S8_state_serum_screen=screen,S8_state_hit_counts=hits,S9_primary_candidate_features=data.frame(feature=candidates),
    S9_CTAG_component_decomposition=components,S9_candidate_performance=cmp$performance,S9_fixed_performance=fixed$performance,
    S9_nested_metrics=nested$metrics,S9_nested_selection_frequency=nested$frequency,S9_candidate_mean_proxy=proxy,
    S9_locked_parameters=pars,S9_locked_clinical_models=clin,S9_clinical_normal_reference_sensitivity=norm,
    S9_crossfitted_clinical_sensitivity=cross,S9_tissue_state_TSLB43=tissue,
    S9_OOF_donor_means=oo[,c('Barcode','target','pred','serum_only_score','Age','BMI','ParousFlag','TSLB_years','design_postpartum_serum')],
    S9_locked_score_donors=scores[,c('Barcode','primary_serum','design_postpartum_serum','locked_score','Age','BMI','ParousFlag','TSLB_years')],
    S10_lobule_serum_screen=ls$screen,S10_hit_counts=lss,Serum_model_membership=membership,MISPA_transformations=transformations,
    Model_registry=do.call(rbind,register),Full_model_coefficients=do.call(rbind,coeff),
    S9_candidate_OOF_predictions=cmp$predictions,S9_fixed_OOF_predictions=fixed$predictions,S9_nested_OOF_predictions=nested$predictions)
  figure4_write_workbook(tables,out);saveRDS(fits,file.path(out,'models/Figure4_S8-S9_fitted_models.rds'))
  writeLines(c(paste('Nominated:',paste(candidates,collapse=', ')),paste('Selected:',cmp$winner),
    'Selection is computed from the screen and cross-validated model comparison.'),file.path(out,'results/SELECTION_STATUS.txt'))
  plots<-list(Fig4a_counts=hits,Fig4b_associations=screen[screen$state=='myeloid'&screen$feature%in%figure4_config$targeted_features,],
    Fig4c_donor_averages=tables$S9_OOF_donor_means,Fig4c_descriptive_line=cg,Fig4c_performance=fixed$performance[fixed$performance$subset==cmp$winner,],
    Fig4d_donors=clinfit$parity$data,Fig4e_donors=clinfit$postpartum$data,Fig4e_adjusted_line_and_HC3_t_interval=curves,Clinical_annotations=clin,
    FigS8a_all_state_results_for_selected_features=screen[screen$adjustment=='primary'&screen$feature%in%top,],FigS8a_row_order=data.frame(feature=top),
    FigS8b_screen=volc,FigS8c_proxy_donors=pr[,c('Barcode','proxy_residual','tissue_residual')],FigS8c_descriptive_line_interval=pg,
    FigS8c_annotation=proxy,FigS8d_error_summary=cmp$performance,FigS8e_component_tests=components,FigS9a_counts=lss,
    FigS9b_candidate_coefficients=ls$screen[ls$screen$feature%in%figure4_config$targeted_features,])
  for(n in names(plots))write_csv(plots[[n]],file.path(out,'plot_data',paste0(n,'.csv')))
  list(data=dat,dsp=dsp,states=states,results=tables,plots=plots)
}

