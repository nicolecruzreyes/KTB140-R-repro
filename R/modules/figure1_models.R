# Figure 1 and S1-S4: calculations from the common input data.
# Centroid-derived endpoints in data_common are inputs; their original 1.62M-cell
# reconstruction is NOT rerun here. Pair B is checked directly from its cells.

model_row <- function(obj) {
  r <- obj$row
  r$n_lobules_or_donors <- r$n
  r$residual_df <- r$df
  r
}

mixed_dual <- function(d, y, adjusted=TRUE) {
  terms <- c('ParousFlag', if (adjusted) c('Age','BMI'))
  dd <- d[finite_rows(d,c('Barcode',y,terms)),,drop=FALSE]
  form <- as.formula(paste(y,'~',paste(terms,collapse=' + '),'+ (1|Barcode)'))
  fit <- lmerTest::lmer(form, data=dd, REML=FALSE,
     control=lme4::lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=200000)))
  sm <- coef(summary(fit, ddf='Satterthwaite'))['ParousFlag',]
  b <- unname(sm['Estimate']); se <- unname(sm['Std. Error']); df <- unname(sm['df'])
  vc <- as.data.frame(lme4::VarCorr(fit))
  iv <- vc$vcov[vc$grp=='Barcode']; rv <- vc$vcov[vc$grp=='Residual']
  msg <- fit@optinfo$conv$lme4$messages
  base <- data.frame(outcome=y,term='ParousFlag',beta=100*b,se=100*se,
    n_lobules_or_donors=nrow(dd),donors=length(unique(dd$Barcode)),
    family=if(grepl('p63neg',y))'denominator_2' else 'composition_3',
    units='percentage points',ICC=iv/(iv+rv),logLik=as.numeric(logLik(fit)),
    donor_variance=iv,residual_variance=rv,
    model=if(adjusted)'ML age/BMI adjusted' else 'ML unadjusted',
    singular=lme4::isSingular(fit),
    convergence=if(is.null(msg))'OK' else paste(msg,collapse='; '),stringsAsFactors=FALSE)
  normal <- base; normal$lower <- 100*(b-qnorm(.975)*se); normal$upper <- 100*(b+qnorm(.975)*se)
  normal$p_raw <- 2*pnorm(-abs(b/se));normal$p_reference <- 'normal';normal$ci_reference <- 'normal'
  satt <- base; satt$lower <- 100*(b-qt(.975,df)*se);satt$upper <- 100*(b+qt(.975,df)*se)
  satt$p_raw <- unname(sm['Pr(>|t|)']);satt$df <- df
  satt$p_reference <- 'Satterthwaite';satt$ci_reference <- 'Satterthwaite'
  list(normal=normal,satterthwaite=satt,fit=fit)
}

bootstrap_from_plan <- function(d, root, model, out) {
  serial <- model!='simple'; raw <- model=='serial_raw_contact'
  mediator <- if(raw)'raw_erp63_contact' else 'erp63_mean'
  order <- read_csv(file.path(root,'resampling/bootstrap_donor_order.csv'))
  # All three decompositions use the same 120 complete donors.
  dd <- d[finite_rows(d,c('Barcode','ParousFlag','Age','BMI','er_pooled','p63_pooled',
                         'erp63_mean','raw_erp63_contact')),,drop=FALSE]
  if(!setequal(dd$Barcode,order$Barcode)||anyDuplicated(dd$Barcode))stop('Bootstrap donor membership differs.')
  dd <- dd[match(order$Barcode,dd$Barcode),,drop=FALSE]
  con <- gzfile(file.path(root,'resampling',paste0(model,'_indices.csv.gz')),'rt')
  on.exit(close(con),add=TRUE)
  plan <- as.matrix(utils::read.csv(con,header=FALSE,check.names=FALSE))
  storage.mode(plan) <- 'integer'
  B <- if(serial)config$bootstrap_serial else config$bootstrap_simple
  if(nrow(plan)!=B||ncol(plan)!=nrow(dd)||anyNA(plan)||any(plan<1L|plan>nrow(dd)))stop('Invalid bootstrap plan: ',model)
  X0 <- model.matrix(~ParousFlag+Age+BMI,dd); Y <- dd$er_pooled; P <- dd$p63_pooled; M <- dd[[mediator]]
  calc <- function(ix) {
    X <- X0[ix,,drop=FALSE];m <- M[ix];y <- Y[ix];p <- P[ix]
    cf <- function(x,v) {o<-lm.fit(x,v);if(o$rank!=ncol(x))stop('Rank-deficient supplied bootstrap draw; record and resolve, do not replace silently.');unname(o$coefficients)}
    tt<-cf(X,y)
    if(!serial){aa<-cf(X,m);bb<-cf(cbind(X,m),y);return(c(tt[2],bb[2],aa[2]*tail(bb,1)))}
    aa<-cf(X,p);ab<-cf(cbind(X,p),m);bb<-cf(cbind(X,p,m),y)
    c(tt[2],bb[2],aa[2]*bb[5],ab[2]*bb[6],aa[2]*ab[5]*bb[6])
  }
  point<-calc(seq_len(nrow(dd))); sims<-matrix(NA_real_,B,length(point))
  for(i in seq_len(B)) {sims[i,]<-calc(plan[i,]);if(i%%5000L==0L)message(model,': ',i,'/',B,' draws')}
  comp<-if(!serial)c('total','direct','indirect_enrichment') else if(raw)
    c('total','direct','indirect_p63_only','indirect_raw_contact_only','serial_p63_raw_contact') else
    c('total','direct','indirect_p63_only','indirect_enrichment_only','serial_p63_enrichment')
  colnames(sims)<-comp
  lo<-apply(sims,2,quantile,.025,type=7);hi<-apply(sims,2,quantile,.975,type=7)
  ps<-apply(sims,2,function(v)min(1,2*min((1+sum(v<=0))/(B+1),(1+sum(v>=0))/(B+1))))
  saveRDS(sims,file.path(out,'results',paste0('bootstrap_draws_',model,'.rds')))
  data.frame(model=model,component=comp,beta_pp=100*point,lower_pp=100*lo,upper_pp=100*hi,
    p_bootstrap_sign=ps,donors=nrow(dd),resamples=B,CI_method='percentile donor bootstrap',
    seed=20260930,plan_source='Saved PCG64 donor indices; 1-based',row.names=NULL)
}

cohort_plot_data <- function(dat,states) {
  c<-dat$clinical;q<-c
  q[['Four states']]<-states$all_four_available[match(q$Barcode,states$Barcode)]
  q[['Four-state serum']]<-states$serum_common[match(q$Barcode,states$Barcode)]
  q$Plasma<-!q$primary_serum
  cols<-c('source_has_MxIF','source_has_DSP','source_has_IO360','primary_serum','Plasma',
          'Four states','Four-state serum','design_postpartum_tissue','design_postpartum_serum')
  score<-rowSums(q[,cols,drop=FALSE]);q<-q[order(-score,q$Barcode),c('Barcode',cols)]
  vs<-c('Age','BMI','ParousFlag','Race','BiradDensity','NumPregnancies','NumLiveBirths','AgeAtFirstBirth','TSLB_years')
  labs<-c('Age','BMI','Parity','Race','BI-RADS density','Number of pregnancies','Number of live births','Age at first birth','Time since last birth')
  miss<-do.call(rbind,lapply(seq_along(vs),function(i){v<-vs[i];ok<-if(v%in%c('AgeAtFirstBirth','TSLB_years'))c$ParousFlag==1 else rep(TRUE,nrow(c));ok[is.na(ok)]<-FALSE
    data.frame(variable=v,label=labs[i],missing=sum(is.na(c[[v]][ok])),denominator=sum(ok),percent=100*sum(is.na(c[[v]][ok]))/sum(ok))}))
  ns<-c(nrow(c),sum(c$source_has_IO360),sum(c$primary_serum),sum(c$source_has_MxIF),sum(c$source_has_DSP),
    sum(c$source_has_MxIF&c$source_has_DSP&c$source_has_IO360),sum(states$all_four_available),sum(states$serum_common),
    sum(c$design_postpartum_tissue),sum(c$design_postpartum_serum))
  labs<-c('All donors / blood samples','IO360','Primary serum','MxIF','GeoMx DSP','All three tissue modalities',
    'All four corrected tissue states','Common four-state serum','Postpartum tissue composition','Postpartum primary serum')
  list(FigS1a=q,FigS1b=miss,FigS1c=data.frame(subset=labs,n=ns))
}

make_interaction_bundle <- function(d,y) {
  dd<-d[d$eligible_erp63_100&is.finite(d$erp63_enrichment_called),,drop=FALSE]
  dd$mean_enrichment<-ave(dd$erp63_enrichment_called,dd$Barcode,FUN=mean)
  dd$within_enrichment<-dd$erp63_enrichment_called-dd$mean_enrichment
  dd$interaction<-dd$within_enrichment*dd$ParousFlag
  obj<-robust_fit(dd,y,c('mean_enrichment','within_enrichment','ParousFlag','interaction','Age','BMI'),'interaction',TRUE)
  x<-seq(min(obj$data$within_enrichment),max(obj$data$within_enrichment),length.out=150)
  # Preserve the approved 150-point display grid; this is not a statistical resampling step.
  lines<-do.call(rbind,lapply(0:1,function(g){z<-obj$data[rep(1,length(x)),,drop=FALSE]
    for(v in c('mean_enrichment','Age','BMI'))z[[v]]<-mean(obj$data[[v]])
    z$within_enrichment<-x;z$ParousFlag<-g;z$interaction<-x*g
    data.frame(x=x,fit_percent=100*predict(obj$fit,newdata=z),group=c('Nulliparous','Parous')[g+1])}))
  list(row=model_row(obj),data=obj$data,lines=lines,fit=obj)
}

recompute_figure1 <- function(root,out) {
  dir.create(file.path(out,'results'),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(out,'plot_data'),recursive=TRUE,showWarnings=FALSE)
  dat<-read_common(root);d<-dat$lobules;a<-dat$donors
  # The shared state builder is used only to derive S1 availability here.
  # No DSP or serum hypothesis tests are run by Figure 1.
  dsp<-build_dsp(root,dat$clinical);states<-build_states(dat,dsp)
  tabs<-list();plots<-cohort_plot_data(dat,states);objs<-list()
  ys<-c('p63_fraction','er_fraction','ki67_fraction','er_p63neg_fraction','ki67_p63neg_fraction')
  mm<-lapply(ys,function(y)mixed_dual(d,y)); names(mm)<-ys
  tabs$parity_ML_mixed_composition<-adjust_family(dplyr::bind_rows(lapply(mm,`[[`,'normal')),'family')
  tabs$parity_ML_mixed_composition_Satterthwaite<-adjust_family(dplyr::bind_rows(lapply(mm,`[[`,'satterthwaite')),'family')
  # Preserve both existing correction columns; raw model inference is displayed in Fig.1b,c.
  tabs$parity_ML_mixed_composition_Satterthwaite<-adjust_family(tabs$parity_ML_mixed_composition_Satterthwaite,'family',method='bonferroni',name='p_bonferroni')
  tabs$parity_ML_mixed_unadjusted<-dplyr::bind_rows(lapply(ys,function(y)mixed_dual(d,y,FALSE)$satterthwaite))
  tabs$parity_clustered_composition<-dplyr::bind_rows(lapply(ys,function(y){r<-model_row(robust_fit(d,y,c('ParousFlag','Age','BMI'),'ParousFlag',TRUE,scale=100));r$family<-if(grepl('p63neg',y))'denominator_2' else 'composition_3';r}))
  tabs$parity_clustered_composition<-adjust_family(tabs$parity_clustered_composition,'family')
  tabs$parity_donor_composition<-dplyr::bind_rows(lapply(c('pooled','mean_lobule','p63negative'),function(version){
    markers<-if(version=='p63negative')c('er','ki67') else c('p63','er','ki67')
    dplyr::bind_rows(lapply(markers,function(marker){y<-if(version=='p63negative')paste0(marker,'_p63neg_pooled') else paste0(marker,'_',version)
      r<-model_row(robust_fit(a,y,c('ParousFlag','Age','BMI'),'ParousFlag',scale=100));r$version<-version;r}))}))
  tabs$parity_donor_composition<-adjust_family(tabs$parity_donor_composition,'version')
  sp<-nn<-list()
  for(pair in c('erp63','kp63','erki')) {
    gate<-switch(pair,erp63='eligible_erp63_100',kp63='eligible_kp63_10_100',erki='eligible_erki_10_100')
    dd<-d[d[[gate]],,drop=FALSE]
    sp[[pair]]<-model_row(robust_fit(dd,paste0(pair,'_enrichment_called'),c('ParousFlag','Age','BMI'),'ParousFlag',TRUE))
    nn[[pair]]<-model_row(robust_fit(dd,paste0(pair,'_nn_logratio'),c('ParousFlag','Age','BMI'),'ParousFlag',TRUE))
  }
  tabs$parity_neighborhood_family<-adjust_family(dplyr::bind_rows(sp))
  tabs$parity_nearest_neighbor_family<-adjust_family(dplyr::bind_rows(nn))
  four<-dplyr::bind_rows(tabs$parity_clustered_composition[tabs$parity_clustered_composition$family=='composition_3',],sp$erp63)
  tabs$parity_four_endpoint_family<-adjust_family(four)
  br<-list()
  for(sub in names(dat$burdens))for(ep in c('p63_high','ER_low','ERp63_high','composite')) {
    obj<-robust_fit(dat$burdens[[sub]],paste0(ep,'_burden'),c('ParousFlag','Age','BMI'),'ParousFlag',weight=paste0(ep,'_eligible_n'),scale=100)
    r<-model_row(obj);r$subset<-sub;r$cutpoints<-'corrected_quartile_rule';br[[length(br)+1]]<-r
  }
  tabs$parity_burdens<-adjust_family(dplyr::bind_rows(br),'subset')
  it<-lapply(c('er_fraction','p63_fraction'),function(y)make_interaction_bundle(d,y));names(it)<-c('er','p63')
  tabs$enrichment_composition_interactions<-adjust_family(dplyr::bind_rows(lapply(it,`[[`,'row')))
  tabs$enrichment_composition_interactions$exposure<-'ParousFlag'
  assoc<-robust_fit(a,'er_pooled',c('erp63_mean','ParousFlag','Age','BMI'),'erp63_mean')
  tabs$donor_enrichment_ER_association<-model_row(assoc)
  s4<-assoc$data
  s4$xr<-resid(lm(erp63_mean~ParousFlag+Age+BMI,s4));s4$yr<-100*resid(lm(er_pooled~ParousFlag+Age+BMI,s4))
  tabs$exploratory_statistical_decompositions<-dplyr::bind_rows(lapply(c('simple','serial_p63','serial_raw_contact'),function(m)bootstrap_from_plan(a,root,m,out)))
  # Additional sensitivity analyses are reported in the results tables, even though
  # they are not additional main/supplementary panels in this module.
  specs<-list('Age + BMI'=character(),'+ lobule p63 fraction'='p63_fraction',
    '+ donor pooled p63 fraction'='donor_p63_pooled','+ ER and p63 fractions'=c('er_fraction','p63_fraction'),
    '+ log actual cell count'='log_n_called','+ log actual area and actual density'=c('log_area_mm2','log_called_density_mm2'),
    '+ call availability'='call_availability')
  sr<-list();gt<-list()
  for(pair in c('erp63','restricted_erp63')) {
    for(lab in names(specs)){dd<-d[d[[paste0('eligible_',pair,'_100')]],,drop=FALSE]
      r<-model_row(robust_fit(dd,paste0(pair,'_enrichment_called'),c('ParousFlag','Age','BMI',specs[[lab]]),'ParousFlag',TRUE));r$population<-pair;r$model<-lab;sr[[length(sr)+1]]<-r}
    for(th in c(50,100))for(pop in c('called','all_record','hybrid_legacy')){dd<-d[d[[paste0('eligible_',pair,'_',th)]],,drop=FALSE]
      r<-model_row(robust_fit(dd,paste0(pair,'_enrichment_',pop),c('ParousFlag','Age','BMI'),'ParousFlag',TRUE));r$population<-pair;r$cell_threshold<-th;r$neighborhood_population<-pop;gt[[length(gt)+1]]<-r}
  }
  tabs$spatial_covariate_sensitivities<-dplyr::bind_rows(sr);tabs$spatial_threshold_population_audit<-dplyr::bind_rows(gt)
  plots$Fig1bc_donor_values<-a
  plots$Fig1bc_mixed_model_annotations<-tabs$parity_ML_mixed_composition_Satterthwaite[match(c('p63_fraction','er_fraction'),tabs$parity_ML_mixed_composition_Satterthwaite$outcome),]
  plots$Fig1d_neighborhood_models<-tabs$parity_neighborhood_family
  plots$Fig1e_burden_models<-tabs$parity_burdens[tabs$parity_burdens$subset%in%c('all_lobules','count_rich'),]
  for(m in c('p63','er')) {ag<-aggregate(d[[paste0(m,'_fraction')]],list(Barcode=d$Barcode,ParousFlag=d$ParousFlag),function(x)c(min=min(x),max=max(x),mean=mean(x)))
    h<-cbind(ag[1:2],as.data.frame(ag$x));plots[[if(m=='p63')'FigS2c' else 'FigS2d']]<-h}
  plots$FigS2e<-a[is.finite(a$erp63_mean),]
  plots$FigS3a<-tabs$parity_nearest_neighbor_family
  plots$`_S3b_comp`<-tabs$parity_four_endpoint_family[match(c('p63_fraction','er_fraction','ki67_fraction'),tabs$parity_four_endpoint_family$outcome),]
  plots$`_S3b_sp`<-tabs$parity_four_endpoint_family[tabs$parity_four_endpoint_family$outcome=='erp63_enrichment_called',]
  plots$FigS3c<-it$er$data;plots$FigS3d<-it$p63$data
  plots$FigS3c_fitted_lines<-it$er$lines;plots$FigS3d_fitted_lines<-it$p63$lines
  plots$FigS4a<-s4
  plots$enrichment_composition_interactions<-tabs$enrichment_composition_interactions
  plots$donor_enrichment_ER_association<-tabs$donor_enrichment_ER_association
  plots$exploratory_statistical_decompositions<-tabs$exploratory_statistical_decompositions
  for(n in names(plots))write_csv(plots[[n]],file.path(out,'plot_data',paste0(n,'.csv')))
  for(n in names(tabs))write_csv(tabs[[n]],file.path(out,'results',paste0(n,'.csv')))
  write_csv(a,file.path(out,'results','MxIF_donors_rebuilt_R.csv'))
  write_csv(d,file.path(out,'results','MxIF_lobules_rebuilt_R.csv'))
  write_csv(data.frame(endpoint=names(dat$thresholds),cutpoint=as.numeric(dat$thresholds)),file.path(out,'results','burden_cutpoints_R.csv'))
  write_csv(states,file.path(out,'results','Tissue_states_R.csv'))
  # Workbook is an inspection copy; CSVs and saved fit objects remain the programmatic outputs.
  wb<-openxlsx::createWorkbook()
  for(i in seq_along(tabs)){nm<-substr(paste0(sprintf('%02d',i),'_',names(tabs)[i]),1,31);openxlsx::addWorksheet(wb,nm);openxlsx::writeData(wb,nm,tabs[[i]]);openxlsx::freezePane(wb,nm,firstRow=TRUE);openxlsx::setColWidths(wb,nm,1:ncol(tabs[[i]]),'auto')}
  openxlsx::saveWorkbook(wb,file.path(out,'results','Figure1_S1-S4_R_results.xlsx'),overwrite=TRUE)
  saveRDS(list(mixed=lapply(mm,`[[`,'fit'),association=assoc,interactions=it),file.path(out,'results','fitted_models.rds'))
  list(plots=plots,tables=tabs,data=dat,states=states)
}
