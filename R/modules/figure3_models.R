# Figure 3 / S6-S7: calculations from the common input data.
figure3_rcs <- function(x) {
  k<-figure3_config$knots;a<-k[1];b<-k[2];c<-k[3]
  (pmax(x-a,0)^3-(c-a)/(c-b)*pmax(x-b,0)^3+(b-a)/(c-b)*pmax(x-c,0)^3)/(c-a)^2
}
figure3_predict <- function(obj,x,values,modifier=NULL) {
  q<-obj$data;terms<-all.vars(stats::delete.response(stats::terms(obj$fit)))
  nd<-q[rep(1L,length(values)),,drop=FALSE]
  for(v in terms)nd[[v]]<-rep(if(is.numeric(q[[v]]))mean(q[[v]]) else q[[v]][1],length(values))
  nd[[x]]<-values;if(!is.null(modifier))nd<-modifier(nd)
  X<-stats::model.matrix(stats::delete.response(stats::terms(obj$fit)),nd)
  b<-as.numeric(X%*%coef(obj$fit));s<-sqrt(pmax(0,rowSums((X%*%obj$V)*X)))
  crit<-if(obj$row$ci_reference=='normal')qnorm(.975) else qt(.975,obj$row$df)
  data.frame(x=values,fit=b,lower=b-crit*s,upper=b+crit*s)
}
figure3_one <- function(d,key,value) {
  i<-which(d[[key]]==value);if(length(i)!=1L)stop('Nonunique/missing Figure 3 row: ',key,'=',value)
  d[i,,drop=FALSE]
}
# Workbook serialization only. Treat 1-d arrays and boolean text consistently,
# as in the tested Figure 2 exporter; never change scientific CSV values.
figure3_write_tables <- function(tables,out) {
  tables<-lapply(tables,function(d) {
    d<-as.data.frame(d);for(k in names(d))if(length(dim(d[[k]]))==1L)d[[k]]<-as.vector(d[[k]])
    d
  })
  stopifnot(all(nchar(names(tables))<=31L),!anyDuplicated(names(tables)))
  wb<-openxlsx::createWorkbook(creator='KTB140 R pipeline')
  header<-openxlsx::createStyle(fontName='Arial',fontSize=10,textDecoration='bold',fgFill='#E8F0F0',wrapText=TRUE)
  for(n in names(tables)) {
    write_csv(tables[[n]],file.path(out,'results',paste0(n,'.csv')))
    openxlsx::addWorksheet(wb,n);openxlsx::writeData(wb,n,tables[[n]],headerStyle=header,keepNA=FALSE)
    openxlsx::freezePane(wb,n,firstRow=TRUE);openxlsx::setColWidths(wb,n,seq_along(tables[[n]]),19)
  }
  file<-file.path(out,'results/Figure3_S6-S7_R_results.xlsx');openxlsx::saveWorkbook(wb,file,overwrite=FALSE)
  audits<-list()
  for(n in names(tables)) {
    a<-tables[[n]];b<-openxlsx::read.xlsx(file,sheet=n,check.names=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE)
    if(!identical(names(a),names(b))||nrow(a)!=nrow(b))stop('Workbook schema mismatch: ',n)
    for(k in names(a)) {
      x<-a[[k]];y<-b[[k]];if(is.factor(x))x<-as.character(x)
      if(!identical(is.na(x),is.na(y)))stop('Workbook missingness mismatch: ',n,' / ',k)
      ok<-!is.na(x)
      eq<-if(is.numeric(x))all(abs(as.numeric(y[ok])-x[ok])<=1e-12+1e-12*abs(x[ok])) else
        if(is.logical(x))identical(flag(x),flag(y)) else
        if(is.character(x)&&is.logical(y)&&all(tolower(x[ok])%in%c('true','false')))identical(tolower(x[ok]),tolower(as.character(y[ok]))) else
        identical(as.character(x[ok]),as.character(y[ok]))
      if(!eq)stop('Workbook value mismatch: ',n,' / ',k)
    }
    audits[[n]]<-data.frame(worksheet=n,rows=nrow(a),columns=ncol(a),pass=TRUE)
  }
  write_csv(do.call(rbind,audits),file.path(out,'checks/workbook_readback.csv'))
  invisible(tables)
}

recompute_figure3 <- function(root,out) {
  for(f in c('checks','results','plot_data','models'))dir.create(file.path(out,f),recursive=TRUE,showWarnings=FALSE)
  dat<-read_common(root);dsp<-build_dsp(root,dat$clinical);states<-build_states(dat,dsp)
  d<-join_one(states,dat$donors);pp<-d[flag(d$design_postpartum_tissue),,drop=FALSE]
  stopifnot(nrow(pp)==48L,all(pp$ParousFlag==1),all(pp$Age<45),all(is.finite(pp$TSLB_years)),
    min(pp$TSLB_years)==2,max(pp$TSLB_years)==10,sum(is.finite(pp$erp63_mean))==44L,sum(is.finite(pp$DSP_ecm))==48L)
  pp$late<-as.integer(pp$TSLB_years>3)
  fits<-list();member<-list();plot<-list();coeff<-list()
  record<-function(id,obj,unit='donor') {
    fits[[id]]<<-obj;dd<-obj$data
    member[[id]]<<-data.frame(model=id,Barcode=dd$Barcode,
      lobule_uid=if(unit=='lobule')dd$lobule_uid else NA_character_,unit=unit,
      weight=if('.w'%in%names(dd))dd$.w else 1,stringsAsFactors=FALSE)
    cf<-coef(obj$fit);coeff[[id]]<<-data.frame(model=id,term=names(cf),beta=unname(cf),
      se=sqrt(diag(obj$V)),n=nrow(dd),df=stats::df.residual(obj$fit),stringsAsFactors=FALSE)
  }
  # Eight nonlegacy summaries: three count-pooled, three equal-lobule means,
  # and two p63-negative-denominator sensitivities. Family corrections stay separate.
  comp<-list()
  sets<-list(pooled_counts=figure3_config$composition,equal_lobule_mean=figure3_config$equal_lobule,
    p63negative_pooled_counts=figure3_config$denominator)
  for(v in names(sets))for(y in sets[[v]]) {
    z<-robust_fit(pp,y,c('TSLB_years','Age','BMI'),'TSLB_years',scale=100);record(paste0('composition_',y),z)
    r<-z$row;r$version<-v;r$family<-if(v=='p63negative_pooled_counts')'denominator_2' else 'composition_3'
    comp[[length(comp)+1L]]<-r
  }
  comp<-adjust_family(do.call(rbind,comp),'version')
  q<-fits$composition_er_pooled
  plot$Fig3a_donors<-q$data;plot$Fig3a_curve<-figure3_predict(q,'TSLB_years',figure3_config$curve_years)
  plot$Fig3a_annotation<-figure3_one(comp,'outcome','er_pooled')
  # Burdens use actual eligible-lobule counts as regression weights. High density
  # stays a separately reported source-definition sensitivity, not a substitute.
  bur<-list();burdata<-list();preds<-list()
  for(s in c('all_lobules','count_rich','density_rich')) {
    dd<-dat$burdens[[s]];dd<-dd[flag(dd$design_postpartum_tissue),,drop=FALSE]
    dd$eligible_n<-dd$ER_low_eligible_n
    z<-robust_fit(dd,'ER_low_burden',c('TSLB_years','Age','BMI'),'TSLB_years',weight='eligible_n',scale=100)
    record(paste0('burden_',s),z);r<-z$row;r$subset<-s;r$cutpoints<-'corrected_quartile_rule'
    g<-figure3_predict(z,'TSLB_years',c(3,7))
    for(i in 1:2)for(f in c('fit','lower','upper')) {
      nm<-if(f=='fit')paste0('pred_',g$x[i],'_percent') else paste0('pred_',g$x[i],'_',f)
      r[[nm]]<-g[[f]][i]*100
    }
    bur[[s]]<-r;burdata[[s]]<-dd
    preds[[s]]<-data.frame(subset=s,TSLB_years=g$x,fit_percent=100*g$fit,lower_percent=100*g$lower,upper_percent=100*g$upper)
  }
  bur<-do.call(rbind,bur);bur$p_holm<-NA_real_;idx<-bur$subset%in%c('all_lobules','count_rich');bur$p_holm[idx]<-p.adjust(bur$p_raw[idx],'holm')
  for(i in 1:2) {
    nm<-c('Fig3b','Fig3c')[i];s<-c('all_lobules','count_rich')[i];z<-fits[[paste0('burden_',s)]]
    plot[[paste0(nm,'_donors')]]<-z$data;plot[[paste0(nm,'_curve')]]<-figure3_predict(z,'TSLB_years',figure3_config$curve_years)
    plot[[paste0(nm,'_annotation')]]<-figure3_one(bur,'subset',s)
  }
  # The entire four-endpoint family is fitted even though the summary displays
  # only neighborhood and ECM; displayed-only correction would be incorrect.
  tim<-list()
  for(y in figure3_config$neighborhood_dsp) {
    terms<-c('TSLB_years','Age','BMI',if(startsWith(y,'DSP_'))c('mean_epi_score',dsp$instrument_terms))
    z<-robust_fit(pp,y,terms,'TSLB_years');record(paste0('timing_',y),z);tim[[y]]<-z$row
  }
  tim<-adjust_family(do.call(rbind,tim));spl<-list();early<-list()
  for(y in figure3_config$spline) {
    extra<-if(y=='DSP_ecm')c('mean_epi_score',dsp$instrument_terms) else character()
    qq<-pp;qq$rcs_nonlinear<-figure3_rcs(qq$TSLB_years)
    terms<-c('TSLB_years','rcs_nonlinear','Age','BMI',extra)
    z<-robust_fit(qq,y,terms,'rcs_nonlinear');record(paste0('spline_',y),z)
    r<-z$row;tt<-c('TSLB_years','rcs_nonlinear');bv<-coef(z$fit)[tt];V<-z$V[tt,tt,drop=FALSE]
    # HC3 robust joint Wald F test with 2 numerator df.
    Fval<-as.numeric(crossprod(bv,solve(V,bv)))/2
    r$overall_spline_p<-stats::pf(Fval,2,z$row$df,lower.tail=FALSE);r$overall_F<-Fval
    spl[[y]]<-r
    g<-figure3_predict(z,'TSLB_years',figure3_config$curve_years,function(nd){nd$rcs_nonlinear<-figure3_rcs(nd$TSLB_years);nd})
    lin<-if(y=='er_pooled')fits$composition_er_pooled else fits[[paste0('timing_',y)]]
    stopifnot(setequal(lin$data$Barcode,z$data$Barcode));g$linear<-figure3_predict(lin,'TSLB_years',g$x)$fit
    plot[[paste0('Spline_',y,'_donors')]]<-z$data;plot[[paste0('Spline_',y,'_curves')]]<-g
    e<-robust_fit(pp,y,c('late','Age','BMI',extra),'late',scale=if(y=='er_pooled')100 else 1)
    record(paste0('earlylate_',y),e);early[[y]]<-e$row;plot[[paste0('Earlylate_',y,'_donors')]]<-e$data
  }
  spl<-adjust_family(do.call(rbind,spl));early<-adjust_family(do.call(rbind,early))
  # Donor-mean and within-donor variables are computed within the eligible
  # postpartum lobules, before fitting either endpoint.
  lob<-dat$lobules;il<-lob[flag(lob$eligible_erp63_100)&flag(lob$design_postpartum_tissue)&is.finite(lob$erp63_enrichment_called),,drop=FALSE]
  il$mean_enrichment<-ave(il$erp63_enrichment_called,il$Barcode,FUN=mean)
  il$within_enrichment<-il$erp63_enrichment_called-il$mean_enrichment
  il$interaction<-il$within_enrichment*il$TSLB_years
  stopifnot(nrow(il)==249L,length(unique(il$Barcode))==44L)
  ints<-list()
  for(y in figure3_config$interactions) {
    z<-robust_fit(il,y,c('mean_enrichment','within_enrichment','TSLB_years','interaction','Age','BMI'),'interaction',cluster=TRUE)
    record(paste0('interaction_',y),z,'lobule');r<-z$row;r$exposure<-'TSLB_years';ints[[y]]<-r
    plot[[paste0('Interaction_',y,'_lobules')]]<-z$data
    gs<-lapply(c(3,7),function(t) {
      g<-figure3_predict(z,'within_enrichment',seq(min(z$data$within_enrichment),max(z$data$within_enrichment),length.out=161L),
        function(nd){nd$TSLB_years<-t;nd$interaction<-nd$within_enrichment*t;nd})
      g$TSLB_years<-t;g
    });plot[[paste0('Interaction_',y,'_curves')]]<-do.call(rbind,gs)
  }
  ints<-adjust_family(do.call(rbind,ints))
  plot$Spline_annotations<-spl;plot$Earlylate_annotations<-early;plot$Interaction_annotations<-ints
  selected<-list(figure3_one(comp,'outcome','er_pooled'),figure3_one(bur,'subset','all_lobules'),figure3_one(bur,'subset','count_rich'),
    figure3_one(tim,'outcome','erp63_mean'),figure3_one(tim,'outcome','DSP_ecm'),figure3_one(ints,'outcome','p63_fraction'))
  labels<-c('ER-positive fraction','ER-low: all lobules','ER-low: high-cell-count lobules','ER–p63 neighborhood enrichment','DSP ECM/fibroblast signature','Within-donor p63 coupling × TSLB')
  axes<-c(rep('Percentage points/year',3),'Log-enrichment units/year','Mean-z signature units/year','Change in enrichment–fraction slope/year')
  plot$Fig3d_summary<-do.call(rbind,lapply(seq_along(selected),function(i) {
    r<-selected[[i]];data.frame(label=labels[i],beta=r$beta,lower=r$lower,upper=r$upper,p_raw=r$p_raw,p_holm=r$p_holm,donors=r$donors,axis_units=axes[i])
  }))
  hm<-list()
  addholm<-function(d,fam,key) {
    q<-d[order(d$p_raw,d[[key]]),,drop=FALSE];m<-nrow(q);prod<-q$p_raw*(m:1)
    hm[[fam]]<<-data.frame(family=fam,endpoint=q[[key]],rank=seq_len(m),p_raw=q$p_raw,multiplier=m:1,step_product=prod,p_holm=pmin(1,cummax(prod)))
  }
  for(v in names(sets))addholm(comp[comp$version==v,],v,'outcome')
  addholm(bur[idx,],'ER_low_burdens','subset');addholm(tim,'neighborhood_DSP','outcome');addholm(spl,'nonlinearity','outcome');addholm(early,'earlylate','outcome');addholm(ints,'interactions','outcome')
  tables<-list(S7_composition=comp,S7_ER_low_burdens=bur,S7_neighborhood_DSP_family=tim,S7_spline_family=spl,S7_earlylate_family=early,
    S7_interactions=ints,S7_burden_predictions=do.call(rbind,preds),S7_model_membership=do.call(rbind,member),
    S7_Holm_calculations=do.call(rbind,hm),Postpartum_donors=pp,Postpartum_burdens=do.call(rbind,burdata),
    Interaction_lobules=il,Full_model_coefficients=do.call(rbind,coeff),Summary_Fig3d=plot$Fig3d_summary)
  figure3_write_tables(tables,out)
  saveRDS(fits,file.path(out,'models/Figure3_S6-S7_fitted_models.rds'))
  for(n in names(plot))write_csv(plot[[n]],file.path(out,'plot_data',paste0(n,'.csv')))
  # Preserve freshly built upstream scores to check their reuse, not to substitute
  # stored Figure 2 donor states for the current computation.
  write_csv(dsp$donor,file.path(out,'checks/rebuilt_DSP_donors.csv'))
  write_csv(states,file.path(out,'checks/rebuilt_Tissue_states.csv'))
  list(data=dat,dsp=dsp,states=states,results=tables,plots=plot)
}

