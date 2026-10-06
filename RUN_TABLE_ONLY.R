#!/usr/bin/env Rscript
# Fits the 34 models reported only in the supplementary tables (S3B-S3E, S4, S7 slopes).
# No bootstrap, cross-validation or figure drawing is run here.
#
# Usage (from the project root):
#   Rscript RUN_TABLE_ONLY.R [--run-id my_tables] [--root PATH]
arg <- commandArgs(TRUE);opts<-list(root=NULL,`run-id`=NULL)
i<-1L
while(i<=length(arg)){
  k<-sub('^--','',arg[i]);if(!k%in%names(opts)||i==length(arg))stop('Usage: Rscript RUN_TABLE_ONLY.R [--run-id <new_id>] [--root <project>]')
  opts[[k]]<-arg[i+1L];i<-i+2L
}
a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(!is.null(opts$root))normalizePath(opts$root,mustWork=TRUE) else if(length(a))dirname(normalizePath(sub('^--file=','',a[1]),mustWork=TRUE)) else getwd()
runid<-if(is.null(opts$`run-id`))paste0('tables_',format(Sys.time(),'%Y%m%d_%H%M%S')) else opts$`run-id`
if(!grepl('^[A-Za-z0-9_-]+$',runid))stop('Use a simple run-id (letters, numbers, "_" or "-").')
out<-file.path(root,'outputs',runid);if(file.exists(out))stop('Output folder already exists; choose a new --run-id: ',out)
dir.create(out,recursive=TRUE)
for(d in c('checks','results','models','source_tables'))dir.create(file.path(out,d),showWarnings=FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
source(file.path(root,'config/analysis_config.R'),encoding='UTF-8')
for(f in c('00_project.R','01_numeric.R','02_data.R'))source(file.path(root,'R/core',f),encoding='UTF-8')
source(file.path(root,'R/modules/tableonly_models.R'),encoding='UTF-8')
need<-c('lme4','lmerTest','sandwich','dplyr','jsonlite','openxlsx')
miss<-need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss))stop('Missing packages: ',paste(miss,collapse=', '),'. Run 00_INSTALL_PACKAGES.R --yes.')
log<-file(file.path(out,'execution.log'),'wt');sink(log,split=TRUE)
warnings<-list();start<-Sys.time()
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
err<-NULL
result<-tryCatch(withCallingHandlers({
  # Cohort definitions used by the supplementary tables.
  d <- read_common(root);L <- tableonly_prepare(d)
  ok <- finite_rows(L,c('p63_fraction','ParousFlag','Age','BMI','NonwhiteFlag','HispanicFlag'))
  stopifnot(nrow(d$clinical)==140L,nrow(L)==979L,nrow(d$donors)==127L,
    sum(ok)==974L,length(unique(L$Barcode[ok]))==126L,
    sum(L$eligible_erp63_100)==727L,sum(L$eligible_kp63_10_100)==447L,sum(L$eligible_erki_10_100)==481L)
  cat('Fitting the 34 supplementary-table models.\n')
  r<-recompute_tableonly(root,out)
  cat('Supplementary-table models completed.\n')
  TRUE
},warning=function(w){
  warnings[[length(warnings)+1L]]<<-data.frame(message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=' '))
  cat('RECORDED WARNING: ',conditionMessage(w),'\n',sep='');invokeRestart('muffleWarning')
}),error=function(e){err<<-conditionMessage(e);cat('STOPPED: ',err,'\n',sep='');FALSE})
warn<-if(length(warnings))do.call(rbind,warnings) else data.frame(message=character(),call=character())
write_csv(warn,file.path(out,'warnings.csv'))
ktb_write_status(file.path(out,'RUN_STATUS.json'),if(result)'completed' else 'failed',
  list(module='tables',error=err,warnings=nrow(warn),elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),
    platform=Sys.info()[['sysname']],R_version=R.version.string))
sink();close(log)
if(!result)stop(err,call.=FALSE)
cat('Results: ',out,'\n',sep='')
