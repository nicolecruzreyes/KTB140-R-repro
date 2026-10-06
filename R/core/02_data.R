# A single, explicit input layer replaces the five duplicated source-workbook builders.
read_common <- function(root) {
  indir<-file.path(root,config$input_dir)
  clin<-read_csv(file.path(indir,'donors.csv'))
  assert_cols(clin,c('Barcode','Age','BMI','ParousFlag','TSLB_years','primary_serum','design_postpartum_tissue'))
  stopifnot(!anyDuplicated(clin$Barcode),nrow(clin)==140L)
  for(v in grep('^source_has_|^primary_serum$|^design_postpartum',names(clin),value=TRUE))clin[[v]]<-flag(clin[[v]])
  clin$parity_label<-factor(ifelse(is.na(clin$ParousFlag),NA,ifelse(clin$ParousFlag==1,'Parous','Nulliparous')),
                         levels=c('Nulliparous','Parous'))
  lob<-read_csv(file.path(indir,'mxif_lobules.csv'))
  needed<-c('Barcode','lobule_uid','n_epithelial_records','n_called','n_uncalled','area_mm2',
    'n_er','n_p63','n_ki67','n_erk','n_erp63','n_kp63','n_triple','n_p63neg','n_er_p63neg','n_ki67_p63neg',
    'erp63_enrichment_called','kp63_enrichment_called','erki_enrichment_called','restricted_erp63_enrichment_called')
  assert_cols(lob,needed,'MxIF');stopifnot(!anyDuplicated(lob$lobule_uid),nrow(lob)==979L,
    all(lob$n_called+lob$n_uncalled==lob$n_epithelial_records),sum(lob$n_epithelial_records)==1617202,
    sum(lob$n_called)==1577721,all(lob$area_mm2>0),all(lob$n_called==round(lob$n_called)),
    all(lob$n_er_p63neg==lob$n_er-lob$n_erp63),all(lob$n_ki67_p63neg==lob$n_ki67-lob$n_kp63),
    all(lob$n_p63neg==lob$n_called-lob$n_p63))
  for(m in c('er','p63','ki67','erk')) {
    lob[[paste0(m,'_fraction')]]<-lob[[paste0('n_',m)]]/lob$n_called
    lob[[paste0(m,'_density_mm2')]]<-lob[[paste0('n_',m)]]/lob$area_mm2
  }
  lob$called_density_mm2<-lob$n_called/lob$area_mm2
  lob$all_record_density_mm2<-lob$n_epithelial_records/lob$area_mm2
  lob$call_availability<-lob$n_called/lob$n_epithelial_records
  lob$er_p63neg_fraction<-lob$n_er_p63neg/lob$n_p63neg
  lob$ki67_p63neg_fraction<-lob$n_ki67_p63neg/lob$n_p63neg
  for(th in c(50L,100L)) {
    lob[[paste0('eligible_erp63_',th)]]<-lob$n_er>=th&lob$n_p63>=th
    lob[[paste0('eligible_restricted_erp63_',th)]]<-lob$n_er_p63neg>=th&lob$n_p63>=th
  }
  lob$eligible_kp63_10_100<-lob$n_ki67>=10&lob$n_p63>=100
  lob$eligible_erki_10_100<-lob$n_ki67>=10&lob$n_er>=100
  for(v in c('n_called','area_mm2','called_density_mm2'))lob[[paste0('log_',v)]]<-log(lob[[v]])
  lob$count_rich<-ave(log1p(lob$n_called),lob$Barcode,FUN=function(x)as.integer(x>=quantile(x,2/3,type=7)))==1
  lob$density_rich<-ave(log1p(lob$called_density_mm2),lob$Barcode,FUN=function(x)as.integer(x>=quantile(x,2/3,type=7)))==1
  lob<-join_one(lob,clin)
  counts<-c('n_epithelial_records','n_called','n_uncalled','n_er','n_p63','n_ki67','n_erk','n_erp63','n_kp63',
            'n_triple','n_p63neg','n_er_p63neg','n_ki67_p63neg','area_mm2')
  donor<-aggregate(lob[,counts,drop=FALSE],list(Barcode=lob$Barcode),sum)
  donor<-join_one(donor,clin)
  donor$n_lobules<-as.numeric(table(lob$Barcode)[donor$Barcode])
  for(m in c('er','p63','ki67')) {
    donor[[paste0(m,'_pooled')]]<-donor[[paste0('n_',m)]]/donor$n_called
    donor[[paste0(m,'_mean_lobule')]]<-unname(tapply(lob[[paste0(m,'_fraction')]],lob$Barcode,mean)[donor$Barcode])
  }
  donor$er_p63neg_pooled<-donor$n_er_p63neg/donor$n_p63neg
  donor$ki67_p63neg_pooled<-donor$n_ki67_p63neg/donor$n_p63neg
  donor$call_availability<-donor$n_called/donor$n_epithelial_records
  ok<-lob$eligible_erp63_100&is.finite(lob$erp63_enrichment_called)
  donor$erp63_mean<-unname(tapply(lob$erp63_enrichment_called[ok],lob$Barcode[ok],mean)[donor$Barcode])
  if('raw_erp63_contact'%in%names(lob))donor$raw_erp63_contact<-unname(tapply(lob$raw_erp63_contact[ok],lob$Barcode[ok],mean)[donor$Barcode])
  lob$donor_p63_pooled<-donor$p63_pooled[match(lob$Barcode,donor$Barcode)]
  thresholds<-c(p63_high=unname(quantile(lob$p63_fraction,.75)),ER_low=unname(quantile(lob$er_fraction,.25)),
                ERp63_high=unname(quantile(lob$erp63_enrichment_called[ok],.75)))
  lob$p63_high<-lob$p63_fraction>=thresholds['p63_high'];lob$ER_low<-lob$er_fraction<=thresholds['ER_low']
  lob$ERp63_high<-ifelse(ok,lob$erp63_enrichment_called>=thresholds['ERp63_high'],NA)
  lob$composite<-ifelse(ok,(lob$p63_high+lob$ER_low+lob$ERp63_high)>=2,NA)
  burdens<-lapply(c('all_lobules','count_rich','density_rich'),function(subset) {
    rows<-lapply(donor$Barcode,function(b) {
      dd<-lob[lob$Barcode==b,,drop=FALSE];if(subset!='all_lobules')dd<-dd[dd[[subset]],,drop=FALSE]
      z<-data.frame(Barcode=b,subset=subset)
      for(ep in c('p63_high','ER_low','ERp63_high','composite')) {
        v<-dd[[ep]];z[[paste0(ep,'_eligible_n')]]<-sum(!is.na(v))
        z[[paste0(ep,'_burden')]]<-if(any(!is.na(v)))mean(v,na.rm=TRUE) else NA_real_
      };z
    });join_one(do.call(rbind,rows),clin)
  });names(burdens)<-c('all_lobules','count_rich','density_rich')
  list(clinical=clin,lobules=lob,donors=donor,burdens=burdens,thresholds=thresholds,
       io360=build_io360(indir),mispa=build_mispa(indir))
}
build_io360 <- function(indir) {
  expr<-read_csv(file.path(indir,'io360_normalized_expression.csv'))
  samples<-read_csv(file.path(indir,'io360_samples.csv'));qc<-read_csv(file.path(indir,'io360_qc.csv'))
  stopifnot(!anyDuplicated(expr$gene),!anyDuplicated(samples$sample_id),!anyDuplicated(qc$Barcode))
  X<-as.matrix(expr[,samples$sample_id,drop=FALSE]);storage.mode(X)<-'double'
  if(any(!is.finite(X))||any(X<0))stop('IO360 normalized matrix contains unavailable/nonfinite/negative values; reconcile before proceeding.')
  X<-log2(X+1);ids<-unique(samples$Barcode)
  mat<-sapply(ids,function(id)rowMeans(X[,samples$Barcode==id,drop=FALSE]));rownames(mat)<-expr$gene
  defs<-list(ER=c('ESR1','EPCAM','ERBB2','BCL2','SFRP1'),
    Myeloid=c('CD68','CSF1R','ITGAX','LYZ','HLA-DRA','HLA-DPA1','CTSS','C1QB','MRC1'),
    Stromal=c('PECAM1','FLT1','PLOD2','CD36','TGFBR2','ZEB2','PPARG'))
  d<-join_one(data.frame(Barcode=ids),qc)
  for(nm in names(defs)) {
    if(!all(defs[[nm]]%in%rownames(mat)))stop('Missing predefined IO360 genes: ',nm)
    M<-mat[defs[[nm]],,drop=FALSE];Z<-t(apply(M,1,stdz))
    d[[paste0('RNA_',nm,'_raw')]]<-colMeans(Z)
    d[[paste0('RNA_',nm)]]<-stdz(residualize(d,paste0('RNA_',nm,'_raw'),c('binding_density','fov_pct_counted')))
  }
  list(donor=d,expression_log2=mat,defs=defs,samples=samples)
}
build_mispa <- function(indir) {
  raw<-read_csv(file.path(indir,'mispa_raw.csv'));map<-read_csv(file.path(indir,'mispa_feature_map.csv'))
  feats<-map$feature_id;assert_cols(raw,c('Barcode','Sample_Type',feats));stopifnot(length(feats)==209L,!anyDuplicated(raw$Barcode))
  logd<-raw[,c('Barcode','Sample_Type'),drop=FALSE];eps<-numeric(length(feats));names(eps)<-feats
  for(f in feats) {
    x<-as.numeric(raw[[f]]);if(any(!is.finite(x))||any(x<0))stop('Invalid MISPA values for ',f)
    pos<-x[x>0];if(!length(pos))stop('No positive MISPA values for ',f)
    eps[f]<-min(pos)/2;logd[[f]]<-log2(x+eps[f])
  }
  list(raw=raw,log=logd,features=feats,map=map,pseudocounts=eps)
}
dsp_definitions <- function()list(
  immune=c('CD45','CD3','CD8','CD4','CD20','CD56'),
  checkpoint=c('PD.1','PD.L1','PD.L2','CTLA4','LAG3','Tim.3','VISTA','B7.H3'),
  myeloid=c('CD68','CD163','CD14','CD11c','HLA.DR','ARG1','IDO1'),
  ecm=c('FAP.alpha','SMA','Fibronectin','CD34'))
build_dsp <- function(root,clinical,review=FALSE) {
  if(!review) {
    if(config$dsp_status!='source_reconciled')stop('Final run requires the source-reconciled GeoMx DSP export.')
    p<-read_csv(file.path(root,config$dsp_provenance_file));assert_cols(p,c('item','confirmed','detail'))
    required<-c('normalized_values','identifier_reconciliation','instrument_metadata','retention_policy')
    if(!all(required%in%p$item)||!all(flag(p$confirmed[match(required,p$item)])))stop('Required source reconciliation has not been completed; complete DSP source reconciliation before analysis.')
    d<-read_csv(file.path(root,config$dsp_file));assert_cols(d,c('roi_key','Barcode','included',config$dsp_instrument_column))
    d<-d[flag(d$included),,drop=FALSE];status<-'GeoMx DSP export; source-reconciled; all supplied study ROIs retained'
  } else {d<-read_csv(file.path(root,config$dsp_review_file));status<-'PROVISIONAL: unreconciled DSP matrix; instrument omitted'}
  assert_cols(d,c('roi_key','Barcode','AOI.nuclei.count','PanCk','EpCAM','Ms.IgG2a','Rb.IgG'))
  stopifnot(!anyDuplicated(d$roi_key));defs<-dsp_definitions()
  targets<-unique(c(unlist(defs),'PanCk','EpCAM','Ms.IgG2a','Rb.IgG'))
  assert_cols(d,targets);for(v in targets)if(any(!is.finite(d[[v]])|d[[v]]<=0))stop('DSP log2 requires positive finite values; reconcile target ',v,' rather than add an undocumented pseudocount.')
  if(any(!is.finite(d$AOI.nuclei.count)|d$AOI.nuclei.count<=0))stop('Nonpositive or missing DSP nuclei weights.')
  # Clinical covariates come from one donor source, never competing ROI columns.
  d<-d[,setdiff(names(d),c('Age','BMI','ParousFlag','parity_label','TSLB_years')),drop=FALSE];d<-join_one(d,clinical)
  if(anyNA(d$Age)|anyNA(d$BMI))stop('Unmatched DSP donors.')
  d$epi_score<-(stdz(log2(d$PanCk))+stdz(log2(d$EpCAM)))/2
  for(s in names(defs))d[[paste0('DSP_',s)]]<-rowMeans(sapply(defs[[s]],function(v)stdz(log2(d[[v]]))))
  instcols<-character()
  if(!review) {
    instrument<-as.character(d[[config$dsp_instrument_column]])
    if(anyNA(instrument)|any(!nzchar(instrument)))stop('Missing confirmed instrument labels.')
    lev<-sort(unique(instrument))
    if(length(lev)>1)for(k in seq_along(lev)[-1]) {
      nm<-paste0('instrument_',k);d[[nm]]<-as.numeric(instrument==lev[k]);instcols<-c(instcols,nm)
    }
  }
  dr<-lapply(unique(d$Barcode),function(b) {
    q<-d[d$Barcode==b,,drop=FALSE];w<-q$AOI.nuclei.count
    r<-data.frame(Barcode=b,n_ROIs=nrow(q),mean_epi_score=mean(q$epi_score),
        mean_log2_MsIgG2a=mean(log2(q$Ms.IgG2a)),mean_MsIgG2a=mean(q$Ms.IgG2a),mean_RbIgG=mean(q$Rb.IgG))
    for(s in names(defs))r[[paste0('DSP_',s)]]<-safe_weighted_mean(q[[paste0('DSP_',s)]],w)
    for(v in instcols)r[[v]]<-safe_weighted_mean(q[[v]],w)
    r
  });donor<-join_one(do.call(rbind,dr),clinical)
  for(s in names(defs))donor[[paste0('DSP_',s,'_adj')]]<-stdz(residualize(donor,paste0('DSP_',s),c('mean_epi_score','mean_log2_MsIgG2a',instcols)))
  donor$DSP_block<-(donor$DSP_immune_adj+donor$DSP_myeloid_adj)/2
  list(roi=d,donor=donor,defs=defs,instrument_terms=instcols,status=status)
}
build_states <- function(dat,dsp) {
  d<-dat$clinical
  d<-join_one(d,dat$donors[,c('Barcode','er_pooled','p63_pooled','erp63_mean')])
  d<-join_one(d,dat$io360$donor)
  d$epithelial_balance<-stdz((stdz(d$er_pooled)-stdz(d$p63_pooled)+d$RNA_ER)/3)
  d$neighborhood<-stdz(d$erp63_mean)
  if(is.null(dsp)) {
    d$myeloid<-NA_real_;d$stromal<-NA_real_
    d$all_four_available<-NA;d$serum_common<-NA
    d$anticipated_all_four<-is.finite(d$epithelial_balance)&is.finite(d$neighborhood)&d$source_has_DSP
    d$anticipated_serum_common<-d$anticipated_all_four&d$primary_serum&is.finite(d$ParousFlag)
  } else {
    d<-join_one(d,dsp$donor)
    # Strict missingness: both assay blocks must be available; never rowMeans(...,na.rm=TRUE).
    d$myeloid<-stdz((d$DSP_block+d$RNA_Myeloid)/2)
    d$stromal<-stdz((d$DSP_ecm_adj+d$RNA_Stromal)/2)
    d$all_four_available<-apply(d[,names(state_labels)],1,function(x)all(is.finite(x)))
    d$serum_common<-d$all_four_available&d$primary_serum&finite_rows(d,c('Age','BMI','ParousFlag'))
  };d
}
