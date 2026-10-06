# Models reported only in the supplementary tables (S3B-S3E, S4, S7 slopes).
# Uses the same endpoints, covariance and correction families; no figure fitting,
# bootstrap, cross-validation, cell reconstruction, or spatial randomization.

tableonly_prepare <- function(dat) {
  d <- dat$lobules
  assert_cols(d,c('Race','Hispanic','p63_fraction','er_fraction','ki67_fraction','log_n_called'))
  race <- tolower(trimws(as.character(d$Race)))
  his <- tolower(trimws(as.character(d$Hispanic)))
  # Exact original sensitivity coding. This is not a new ancestry classification.
  d$NonwhiteFlag <- ifelse(is.na(d$Race)|race=='',NA_real_,as.numeric(race!='white'))
  d$HispanicFlag <- ifelse(his=='yes',1,ifelse(his=='no',0,NA_real_))
  stopifnot(all(is.finite(d$log_n_called)),all(abs(d$log_n_called-log(d$n_called))<1e-12))
  d
}

tableonly_mixed <- function(d,y,terms,reml=FALSE) {
  dd <- d[finite_rows(d,c('Barcode',y,terms)),,drop=FALSE]
  f <- as.formula(paste(y,'~',paste(terms,collapse=' + '),'+ (1|Barcode)'))
  X <- model.matrix(reformulate(terms),dd)
  if(qr(X)$rank!=ncol(X))stop('Rank-deficient table-only mixed model: ',deparse(f))
  fit <- lmerTest::lmer(f,data=dd,REML=reml,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000)))
  vc <- as.data.frame(lme4::VarCorr(fit));iv <- vc$vcov[vc$grp=='Barcode'];rv <- vc$vcov[vc$grp=='Residual']
  if(length(iv)!=1L||length(rv)!=1L||!is.finite(iv)||!is.finite(rv)||iv<0||rv<=0)stop('Invalid variance components.')
  cm <- fit@optinfo$conv$lme4$messages;code <- fit@optinfo$conv$opt
  ok <- is.null(cm)&&(is.null(code)||all(code==0))
  sm <- coef(summary(fit,ddf='Satterthwaite'));r <- sm['ParousFlag',]
  b <- unname(r['Estimate']);se <- unname(r['Std. Error']);df <- unname(r['df'])
  base <- data.frame(outcome=y,beta=100*b,se=100*se,n=nrow(dd),donors=length(unique(dd$Barcode)),
    donor_variance=iv,residual_variance=rv,ICC=iv/(iv+rv),fit=if(reml)'REML' else 'ML',
    singular=lme4::isSingular(fit),convergence=if(ok)'OK' else paste(c(cm,code),collapse='; '),
    units='percentage points',stringsAsFactors=FALSE)
  normal <- base;normal$lower <- 100*(b-qnorm(.975)*se);normal$upper <- 100*(b+qnorm(.975)*se)
  normal$p_raw <- 2*pnorm(-abs(b/se));normal$p_reference <- 'normal';normal$ci_reference <- 'normal'
  sat <- base;sat$lower <- 100*(b-qt(.975,df)*se);sat$upper <- 100*(b+qt(.975,df)*se)
  sat$p_raw <- unname(r['Pr(>|t|)']);sat$df <- df;sat$p_reference <- 'Satterthwaite';sat$ci_reference <- 'Satterthwaite'
  list(fit=fit,data=dd,V=as.matrix(vcov(fit)),normal=normal,satterthwaite=sat,
    row=base,formula=paste(deparse(f),collapse=' '))
}

# Numerics-only workbook writer; no graph/font-device dependency is required.
# Uses the established openxlsx read-back conventions, including 1-d arrays.
tableonly_write_workbook <- function(tables,out) {
  wb <- openxlsx::createWorkbook(creator='KTB140 R table-only completion')
  openxlsx::modifyBaseFont(wb,fontName='Arial',fontSize=10)
  head <- openxlsx::createStyle(fontName='Arial',fontSize=10,textDecoration='bold',fgFill='#E8F0F0',wrapText=TRUE)
  saved <- list()
  for(n in names(tables)) {
    if(nchar(n)>31L)stop('Worksheet name too long: ',n)
    d <- as.data.frame(tables[[n]],stringsAsFactors=FALSE)
    for(k in names(d)){if(length(dim(d[[k]]))==1L)d[[k]]<-as.vector(d[[k]]);if(is.factor(d[[k]]))d[[k]]<-as.character(d[[k]])}
    saved[[n]] <- d
    write_csv(d,file.path(out,'results',paste0(n,'.csv')))
    openxlsx::addWorksheet(wb,n);openxlsx::writeData(wb,n,d,headerStyle=head,keepNA=FALSE)
    openxlsx::freezePane(wb,n,firstRow=TRUE)
    widths <- pmax(18,pmin(28,nchar(names(d))+2))
    for(k in intersect(c('outcome','model','model_id','formula','detail','term','item'),names(d)))
      widths[match(k,names(d))] <- switch(k,outcome=22,model=44,model_id=48,formula=64,detail=94,term=30,item=24)
    openxlsx::setColWidths(wb,n,seq_along(d),widths)
    openxlsx::setRowHeights(wb,n,1,32)
    if(nrow(d)){
      openxlsx::addStyle(wb,n,openxlsx::createStyle(valign='center',wrapText=TRUE),rows=2:(nrow(d)+1),cols=seq_along(d),gridExpand=TRUE,stack=TRUE)
      numeric_cols <- which(vapply(d,is.numeric,logical(1)))
      if(length(numeric_cols))openxlsx::addStyle(wb,n,openxlsx::createStyle(numFmt='0.000000'),rows=2:(nrow(d)+1),cols=numeric_cols,gridExpand=TRUE,stack=TRUE)
      for(j in numeric_cols) {
        is_count <- names(d)[j]%in%c('n','donors','lobules','n_lobules_or_donors','fixed_rank','clusters','at_low','at_high')
        if(is_count)openxlsx::addStyle(wb,n,openxlsx::createStyle(numFmt='0'),rows=2:(nrow(d)+1),cols=j,gridExpand=TRUE,stack=TRUE)
        # Preserve nonzero small values, especially P values, in the display.
        small <- which(is.finite(d[[j]]) & d[[j]]!=0 & abs(d[[j]])<1e-6)
        if(grepl('^p_',names(d)[j]))small <- seq_len(nrow(d))
        if(length(small))openxlsx::addStyle(wb,n,openxlsx::createStyle(numFmt='0.000000E+00'),rows=small+1L,cols=j,gridExpand=TRUE,stack=TRUE)
      }
      # Explicit heights avoid clipped wrapped descriptions across spreadsheet readers.
      lines <- vapply(seq_len(nrow(d)),function(i)max(vapply(seq_along(d),function(j){
        s <- if(is.na(d[[j]][i]))'' else as.character(d[[j]][i])
        ceiling(max(1,nchar(s))/(widths[j]*0.8))
      },numeric(1))),numeric(1))
      openxlsx::setRowHeights(wb,n,2:(nrow(d)+1),pmax(20,15*lines+5))
    }
    # openxlsx creates unused drawing/VML relationships even for plain tables.
    # Strict OOXML readers reject these absent parts. There is no artwork here.
    i <- length(saved)
    stopifnot(length(wb$drawings[[i]])==0L,length(wb$vml[[i]])==0L)
    wb$worksheets_rels[[i]] <- wb$worksheets_rels[[i]][!grepl('/(drawing|vmlDrawing)"',wb$worksheets_rels[[i]])]
    wb$worksheets[[i]]$drawing <- character(0)
  }
  f <- file.path(out,'results/KTB140_table_only_R_results.xlsx')
  openxlsx::saveWorkbook(wb,f,overwrite=FALSE)
  audit <- list()
  for(n in names(saved)) {
    a <- saved[[n]];b <- openxlsx::read.xlsx(f,sheet=n,check.names=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE)
    if(!identical(names(a),names(b))||nrow(a)!=nrow(b))stop('Workbook schema differs: ',n)
    for(k in names(a)) {
      x <- a[[k]];y <- b[[k]]
      if(!identical(as.vector(is.na(x)),as.vector(is.na(y))))stop('Workbook missingness differs: ',n,'/',k)
      use <- !is.na(x)
      good <- if(is.numeric(x))all(abs(as.numeric(y[use])-x[use])<=1e-12+1e-12*abs(x[use])) else
        if(is.logical(x))identical(flag(x),flag(y)) else
        if(is.character(x)&&is.logical(y)&&all(tolower(x[use])%in%c('true','false')))identical(tolower(x[use]),tolower(as.character(y[use]))) else
        identical(as.character(x[use]),as.character(y[use]))
      if(!good)stop('Workbook differs from R table: ',n,'/',k)
    }
    audit[[n]] <- data.frame(worksheet=n,rows=nrow(a),columns=ncol(a),pass=TRUE)
  }
  write_csv(do.call(rbind,audit),file.path(out,'checks/workbook_readback.csv'))
  invisible(saved)
}

recompute_tableonly <- function(root,out) {
  for(n in c('results','models','checks','source_tables'))dir.create(file.path(out,n),recursive=TRUE,showWarnings=FALSE)
  dat <- read_common(root);L <- tableonly_prepare(dat);clinical <- dat$clinical
  fits <- list();reg <- list();member <- list();coeff <- list();covars <- list();sensitivity <- list();normals <- list();iccs <- list()
  record <- function(id,obj,type,unit) {
    fits[[id]] <<- obj
    dd <- obj$data;cf <- if(type=='mixed')lme4::fixef(obj$fit) else coef(obj$fit);V <- as.matrix(obj$V)
    reg[[id]] <<- data.frame(model_id=id,type=type,unit=unit,formula=paste(deparse(formula(obj$fit)),collapse=' '),
      n=nrow(dd),donors=length(unique(dd$Barcode)),fixed_rank=length(cf),
      weight=if('.w'%in%names(dd))'eligible-lobule count' else 'none',
      convergence=if(type=='mixed')obj$row$convergence else 'full-rank',
      singular=if(type=='mixed')obj$row$singular else FALSE,stringsAsFactors=FALSE)
    member[[id]] <<- data.frame(model_id=id,Barcode=dd$Barcode,unit_id=if(unit=='lobule')dd$lobule_uid else dd$Barcode,
      weight=if('.w'%in%names(dd))dd$.w else 1,stringsAsFactors=FALSE)
    coeff[[id]] <<- data.frame(model_id=id,term=names(cf),beta=unname(cf),se=sqrt(diag(V)),stringsAsFactors=FALSE)
    vv <- expand.grid(term_row=names(cf),term_col=names(cf),stringsAsFactors=FALSE)
    vv$model_id <- id;vv$covariance <- as.vector(V);covars[[id]] <<- vv[,c('model_id','term_row','term_col','covariance')]
    saveRDS(obj,file.path(out,'models',paste0(id,'.rds')))
  }
  labels <- c(race_ethnicity='Age/BMI + race/ethnicity',cell_count='Age/BMI + log actual called-cell count',unadjusted='Unadjusted')
  # Table S3B: three established specifications for each composition endpoint.
  # The unadjusted fit is one of the three reported specifications.
  for(y in c('p63_fraction','er_fraction','ki67_fraction'))for(spec in names(labels)) {
    terms <- switch(spec,race_ethnicity=c('ParousFlag','Age','BMI','NonwhiteFlag','HispanicFlag'),
      cell_count=c('ParousFlag','Age','BMI','log_n_called'),unadjusted='ParousFlag')
    z <- tableonly_mixed(L,y,terms,FALSE);id <- paste('S3B',spec,y,sep='_');record(id,z,'mixed','lobule')
    a <- z$normal;a$model <- labels[[spec]];a$model_id <- id;normals[[id]] <- a
    b <- z$satterthwaite;b$model <- labels[[spec]];b$model_id <- id;sensitivity[[id]] <- b
  }
  # Table S3C is adjusted REML, not the primary ML fit and not an intercept-only ICC.
  for(y in c('p63_fraction','er_fraction','ki67_fraction')) {
    z <- tableonly_mixed(L,y,c('ParousFlag','Age','BMI'),TRUE);id <- paste0('S3C_REML_',y);record(id,z,'mixed','lobule')
    r <- z$row[,c('outcome','fit','ICC','donor_variance','residual_variance','n','donors','singular','convergence')]
    r$model_id <- id;r$fixed_covariates <- 'ParousFlag + Age + BMI';r$variance_units <- 'fraction squared';iccs[[id]] <- r
  }
  # Table S3D/E: equal weight for donors, means of eligible lobule log endpoints.
  spatial <- list();spdata <- list()
  for(pair in c('erp63','kp63','erki'))for(kind in c('neighborhood','nearest_neighbor')) {
    gate <- switch(pair,erp63='eligible_erp63_100',kp63='eligible_kp63_10_100',erki='eligible_erki_10_100')
    col <- paste0(pair,if(kind=='neighborhood')'_enrichment_called' else '_nn_logratio')
    dd <- L[flag(L[[gate]])&is.finite(L[[col]]),,drop=FALSE]
    aa <- aggregate(dd[[col]],list(Barcode=dd$Barcode),mean);names(aa)[2] <- 'mean_spatial'
    aa <- join_one(aa,clinical)
    z <- robust_fit(aa,'mean_spatial',c('ParousFlag','Age','BMI'),'ParousFlag')
    id <- paste('S3',kind,pair,sep='_');record(id,z,'HC3','donor')
    r <- z$row;r$pair <- pair;r$kind <- kind;r$lobules <- nrow(dd);r$source_endpoint <- col;r$model_id <- id
    spatial[[id]] <- r
    spdata[[id]] <- data.frame(model_id=id,Barcode=z$data$Barcode,mean_spatial=z$data$mean_spatial,
      eligible_lobules=as.integer(table(dd$Barcode)[z$data$Barcode]))
  }
  spatial <- adjust_family(do.call(rbind,spatial),'kind')
  # Tables S4/S7 display summaries: quick refits only, no resampling or new tests.
  burdens <- list()
  for(sub in names(dat$burdens))for(ep in c('p63_high','ER_low','ERp63_high','composite')) {
    y <- paste0(ep,'_burden');z <- robust_fit(dat$burdens[[sub]],y,c('ParousFlag','Age','BMI'),'ParousFlag',weight=paste0(ep,'_eligible_n'),scale=100)
    id <- paste('S4_burden',sub,ep,sep='_');record(id,z,'HC3','donor');r <- z$row
    g <- prediction_grid(z,'ParousFlag',c(0,1));r$subset <- sub;r$endpoint <- ep;r$model_id <- id
    r$model_pred_null_percent <- g$fit[1]*100;r$model_pred_parous_percent <- g$fit[2]*100
    r$mean_Age <- mean(z$data$Age);r$mean_BMI <- mean(z$data$BMI);burdens[[id]] <- r
  }
  burdens <- adjust_family(do.call(rbind,burdens),'subset')
  slopes <- list()
  for(exposure in c('ParousFlag','TSLB_years'))for(y in c('er_fraction','p63_fraction')) {
    use <- flag(L$eligible_erp63_100)&is.finite(L$erp63_enrichment_called)
    if(exposure=='TSLB_years')use <- use&flag(L$design_postpartum_tissue)
    dd <- L[use,,drop=FALSE];dd$mean_enrichment <- ave(dd$erp63_enrichment_called,dd$Barcode,FUN=mean)
    dd$within_enrichment <- dd$erp63_enrichment_called-dd$mean_enrichment;dd$interaction <- dd$within_enrichment*dd[[exposure]]
    z <- robust_fit(dd,y,c('mean_enrichment','within_enrichment',exposure,'interaction','Age','BMI'),'interaction',cluster=TRUE)
    id <- paste('coupling',exposure,y,sep='_');record(id,z,'clustered HC1','lobule')
    at <- if(exposure=='ParousFlag')c(0,1) else c(3,7);cf <- coef(z$fit);r <- z$row
    r$exposure <- exposure;r$slope_low <- unname(cf['within_enrichment']+at[1]*cf['interaction'])
    r$slope_high <- unname(cf['within_enrichment']+at[2]*cf['interaction']);r$at_low <- at[1];r$at_high <- at[2]
    r$model_id <- id;slopes[[id]] <- r
  }
  slopes <- adjust_family(do.call(rbind,slopes),'exposure')
  rd <- do.call(rbind,reg);stopifnot(nrow(rd)==34L,length(fits)==34L)
  notes <- data.frame(item=c('Edition','Scope','S3B inference','S3C specification','Spatial averaging','Burden display','Slopes','Not repeated','Run platform','Publication status'),
    detail=c('R calculations from the common input data',
      '9 ML sensitivities, 3 adjusted REML ICCs, 6 donor spatial models; 12 burden prediction and 4 slope refits',
      'Satterthwaite P/CI are supplied with separate normal-Wald references; raw sensitivity P values, no new correction family',
      'Parity + age + BMI and donor random intercept; REML variance components, not primary fixed-effect inference',
      'Mean of eligible lobule log endpoints; unweighted donor regression, HC3 t inference; Holm across three pairs within kind',
      'At arithmetic mean age/BMI within the model subset; WLS uses endpoint-eligible lobule counts',
      'Within-donor fraction/log-enrichment slope; parity at 0/1 and TSLB at 3/7 years',
      'No production bootstrap, CV, serum screen, figure generation, whole-cell reconstruction or conditional randomization',
      paste(R.version.string,Sys.info()[['sysname']]),'Source workbook for the supplementary tables'),stringsAsFactors=FALSE)
  tabs <- list(README=notes,S3B_ML_Satterthwaite=do.call(rbind,sensitivity),S3B_ML_normal_reference=do.call(rbind,normals),
    S3C_REML_ICC=do.call(rbind,iccs),S3DE_donor_spatial=spatial,S4_burden_predictions=burdens,
    S4_S7_local_slopes=slopes,Model_registry=rd,Full_model_coefficients=do.call(rbind,coeff))
  tableonly_write_workbook(tabs,out)
  write_csv(do.call(rbind,covars),file.path(out,'results/Full_covariance_elements.csv'))
  write_csv(do.call(rbind,member),file.path(out,'source_tables/Model_membership.csv'))
  write_csv(do.call(rbind,spdata),file.path(out,'source_tables/Spatial_donor_values.csv'))
  write_csv(L[,c('Barcode','lobule_uid','Race','Hispanic','NonwhiteFlag','HispanicFlag','n_called','log_n_called')],file.path(out,'source_tables/S3B_coding.csv'))
  write_csv(dat$donors,file.path(out,'source_tables/MxIF_donors_rebuilt_R.csv'))
  saveRDS(fits,file.path(out,'models/tableonly_all_models.rds'))
  invisible(list(tables=tabs,models=fits,data=dat))
}
