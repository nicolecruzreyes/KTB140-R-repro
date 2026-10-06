# Project setup helpers: locate the project root, read runner options, check packages.
ktb_find_root <- function(path=NULL) {
  if(!is.null(path))return(normalizePath(path,mustWork=TRUE))
  filearg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
  start<-if(length(filearg))dirname(normalizePath(sub('^--file=','',filearg[1]),mustWork=TRUE)) else getwd()
  for(i in 0:6){if(file.exists(file.path(start,'config/analysis_config.R')))return(normalizePath(start,mustWork=TRUE));start<-dirname(start)}
  stop('Run from the project root (the folder containing config/analysis_config.R), or supply --root.')
}
# Options for RUN_MODULE.R. --module is required so a run never falls back to a default figure.
ktb_module_options <- function(args) {
  allowed<-c('module','root','out','run-id','pdf-device','powerpoint')
  a<-list(module=NULL,root=NULL,out=NULL,`run-id`=NULL,`pdf-device`='cairo',powerpoint='yes')
  i<-1L
  while(i<=length(args)) {
    k<-sub('^--','',args[i]);if(!k%in%allowed||i==length(args))stop('Unknown option or missing value: ',args[i])
    a[[k]]<-args[i+1L];i<-i+2L
  }
  usage<-'Usage: Rscript RUN_MODULE.R --module figure1|figure2|figure3|figure4 [--run-id ID] [--pdf-device cairo|quartz] [--powerpoint yes|no]'
  if(is.null(a$module)||!a$module%in%c('figure1','figure2','figure3','figure4'))stop(usage)
  if(!a$`pdf-device`%in%c('cairo','quartz'))stop('--pdf-device: cairo (all platforms) or quartz (macOS only)')
  if(!a$powerpoint%in%c('yes','no'))stop('--powerpoint: yes or no')
  a
}
ktb_required <- function() c('ggplot2','systemfonts','ragg','svglite','jsonlite','dplyr','sandwich','lme4','lmerTest','openxlsx')
ktb_require <- function() {
  p<-ktb_required();miss<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss))stop('Missing R packages: ',paste(miss,collapse=', '),'. Run 00_INSTALL_PACKAGES.R --yes (only missing packages are installed).')
  if(packageVersion('ggplot2')<'3.4.0')stop('ggplot2 >=3.4.0 is required for linewidth support.')
  invisible(p)
}
ktb_write_status <- function(path,status,details=list()) {
  jsonlite::write_json(c(list(status=status,timestamp_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),details),path,pretty=TRUE,auto_unbox=TRUE,na='null')
}
