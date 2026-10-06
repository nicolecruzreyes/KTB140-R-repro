# Figure 2 / S5: deterministic calculations from the common input data.
figure2_holm_steps <- function(d,group,key) {
  do.call(rbind,lapply(split(d,d[[group]]),function(q) {
    q<-q[order(q$p_raw,q[[key]]),];m<-nrow(q);mult<-(m:1)*q$p_raw
    data.frame(family=q[[group]],endpoint=q[[key]],rank=seq_len(m),p_raw=q$p_raw,
      multiplier=m:1,step_product=mult,p_holm=pmin(1,cummax(mult)),stringsAsFactors=FALSE)
  }))
}

recompute_figure2 <- function(root,out) {
  for(f in c('checks','results','plot_data','models'))dir.create(file.path(out,f),recursive=TRUE,showWarnings=FALSE)
  dat<-read_common(root);dsp<-build_dsp(root,dat$clinical);states<-build_states(dat,dsp)
  stopifnot(nrow(dsp$roi)==724L,nrow(dsp$donor)==124L,
    sum(states$all_four_available)==114L,sum(states$serum_common)==109L)
  res<-list();fits<-list()
  for(s in names(dsp$defs))for(model in c('ROI','donor','isotype')) {
    q<-if(model=='ROI')dsp$roi else dsp$donor
    terms<-c('ParousFlag','Age','BMI',if(model=='ROI')'epi_score' else 'mean_epi_score',dsp$instrument_terms,
      if(model=='isotype')c('mean_MsIgG2a','mean_RbIgG'))
    # Nuclei weighting is in the DONOR SIGNATURE SUMMARY, not in a second donor regression weight.
    obj<-robust_fit(q,paste0('DSP_',s),terms,'ParousFlag',cluster=model=='ROI')
    r<-obj$row;r$signature<-s;r$model<-model
    res[[length(res)+1L]]<-r;fits[[paste(model,s,sep='_')]]<-obj
  }
  parity<-adjust_family(do.call(rbind,res),'model')
  parous<-dsp$donor$mean_epi_score[dsp$donor$ParousFlag==1]
  nul<-dsp$donor$mean_epi_score[dsp$donor$ParousFlag==0]
  tt<-stats::t.test(parous,nul,var.equal=FALSE)
  content<-data.frame(beta=mean(parous)-mean(nul),p_raw=tt$p.value,n=length(parous)+length(nul),
    n_nulliparous=length(nul),n_parous=length(parous),lower=tt$conf.int[1],upper=tt$conf.int[2],
    df=unname(tt$parameter),p_reference='Welch t',stringsAsFactors=FALSE)
  stopifnot(content$n_nulliparous==61L,content$n_parous==63L)
  sr<-lapply(names(state_labels),function(s) {
    o<-robust_fit(states,s,c('ParousFlag','Age','BMI'),'ParousFlag');fits[[paste0('state_',s)]]<<-o
    r<-o$row;r$state<-s;r
  })
  state_parity<-adjust_family(do.call(rbind,sr))
  cc<-states[states$all_four_available,names(state_labels),drop=FALSE]
  corr<-stats::cor(cc,method='pearson');cor_d<-data.frame(state=rownames(corr),corr,check.names=FALSE,row.names=NULL)
  coverage<-data.frame(state=names(state_labels),donors=vapply(states[names(state_labels)],function(x)sum(is.finite(x)),integer(1)))
  stopifnot(identical(as.integer(coverage$donors),c(127L,120L,124L,124L)))
  membership<-states[,c('Barcode','all_four_available','serum_common'),drop=FALSE]
  hm<-rbind(figure2_holm_steps(parity,'model','signature'),
    figure2_holm_steps(transform(state_parity,family='states'),'family','state'))
  # Export adjustment models and intermediates to make state construction inspectable.
  tech<-list();addtech<-function(q,y,terms,layer) {
    m<-stats::lm(stats::reformulate(terms,response=y),q);cf<-coef(m)
    tech[[length(tech)+1L]]<<-data.frame(layer=layer,outcome=y,term=names(cf),beta=as.numeric(cf),
      n=stats::nobs(m),rank=m$rank,stringsAsFactors=FALSE)
    fits[[paste('technical',layer,y,sep='_')]]<<-m
  }
  for(s in names(dsp$defs))addtech(dsp$donor,paste0('DSP_',s),c('mean_epi_score','mean_log2_MsIgG2a',dsp$instrument_terms),'DSP')
  for(s in names(dat$io360$defs))addtech(dat$io360$donor,paste0('RNA_',s,'_raw'),c('binding_density','fov_pct_counted'),'IO360')
  defs<-rbind(do.call(rbind,lapply(names(dsp$defs),function(s)data.frame(layer='DSP',signature=s,member=dsp$defs[[s]]))),
    do.call(rbind,lapply(names(dat$io360$defs),function(s)data.frame(layer='IO360',signature=s,member=dat$io360$defs[[s]]))))
  tables<-list(S5_DSP_parity=parity,S5_epithelial_content=content,S6_state_parity=state_parity,
    S6_state_correlations=cor_d,DSP_ROIs=dsp$roi,DSP_donors=dsp$donor,
    IO360_signatures=dat$io360$donor,Tissue_states=states,State_coverage=coverage,
    State_membership=membership,Holm_calculations=hm,Technical_adjustment_coefficients=do.call(rbind,tech),
    Fixed_signature_members=defs)
  # Some tapply-derived columns retain a one-dimensional array attribute.
  # Export ordinary vectors so openxlsx stores numbers as numeric cells and
  # missingness checks compare values rather than array-versus-vector metadata.
  tables<-lapply(tables,function(d) {
    for(k in names(d))if(length(dim(d[[k]]))==1L)d[[k]]<-as.vector(d[[k]])
    d
  })
  for(n in names(tables))write_csv(tables[[n]],file.path(out,'results',paste0(n,'.csv')))
  write_csv(data.frame(gene=rownames(dat$io360$expression_log2),dat$io360$expression_log2,check.names=FALSE),
    file.path(out,'results/IO360_donor_log2_expression.csv'))
  saveRDS(fits,file.path(out,'models/Figure2_S5_fitted_models.rds'))
  plots<-list(Fig2b=parity[parity$model%in%c('ROI','donor'),],Fig2c=cor_d,Fig2d=state_parity,
    FigS5b=dsp$donor,FigS5b_annotation=content,FigS5c=parity[parity$model%in%c('donor','isotype'),],
    Fig2a_coverage=coverage)
  for(n in names(plots))write_csv(plots[[n]],file.path(out,'plot_data',paste0(n,'.csv')))
  wb<-openxlsx::createWorkbook(creator='KTB140 R pipeline')
  header<-openxlsx::createStyle(fontName='Arial',fontSize=10,textDecoration='bold',fgFill='#E8F0F0',wrapText=TRUE)
  sheet_names<-setNames(names(tables),names(tables));sheet_names['Technical_adjustment_coefficients']<-'Technical_adjustments'
  stopifnot(all(nchar(sheet_names)<=31L),!anyDuplicated(sheet_names))
  write_csv(data.frame(table=names(sheet_names),worksheet=unname(sheet_names)),file.path(out,'results/WORKBOOK_SHEET_MAP.csv'))
  for(n in names(tables)) {
    sn<-sheet_names[[n]]
    openxlsx::addWorksheet(wb,sn);openxlsx::writeData(wb,sn,tables[[n]],headerStyle=header,keepNA=FALSE)
    openxlsx::freezePane(wb,sn,firstRow=TRUE);openxlsx::setColWidths(wb,sn,cols=seq_along(tables[[n]]),widths=19)
  }
  file<-file.path(out,'results/Figure2_S5_R_results.xlsx');openxlsx::saveWorkbook(wb,file,overwrite=FALSE)
  for(n in names(tables)) {
    got<-openxlsx::read.xlsx(file,sheet=sheet_names[[n]],check.names=FALSE,skipEmptyCols=FALSE,skipEmptyRows=FALSE)
    expect<-as.data.frame(tables[[n]])
    if(!identical(names(got),names(expect))||nrow(got)!=nrow(expect))stop('Workbook schema round-trip mismatch: ',n)
    for(k in names(expect)) {
      a<-expect[[k]];b<-got[[k]]
      if(is.factor(a))a<-as.character(a)
      if(!identical(is.na(a),is.na(b)))stop('Workbook missingness mismatch: ',n,' / ',k)
      ok<-!is.na(a)
      if(is.numeric(a)) {if(any(abs(as.numeric(b[ok])-a[ok])>1e-12+1e-12*abs(a[ok])))stop('Workbook value mismatch: ',n,' / ',k)}
      else if(is.logical(a)) {if(!identical(flag(a),flag(b)))stop('Workbook flag mismatch: ',n,' / ',k)}
      else if(is.character(a)&&is.logical(b)&&all(tolower(a[ok])%in%c('true','false'))) {
        # openxlsx reads text such as source "True" as logical TRUE. Compare
        # only these explicit boolean spellings without changing stored data.
        if(!identical(tolower(a[ok]),tolower(as.character(b[ok]))))stop('Workbook boolean-text mismatch: ',n,' / ',k)
      }
      else if(!identical(as.character(a[ok]),as.character(b[ok])))stop('Workbook text mismatch: ',n,' / ',k)
    }
  }
  list(data=dat,dsp=dsp,states=states,results=tables,plots=plots)
}

