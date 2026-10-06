#!/usr/bin/env Rscript
# Runs one figure module: fits the models, writes result tables and plot data,
# and draws the figures, individual panels and (optionally) PowerPoint assemblies.
#
# Usage (from the project root):
#   Rscript RUN_MODULE.R --module figure2 [--run-id my_run] [--pdf-device cairo|quartz] [--powerpoint yes|no]
# Modules: figure1 (Figure 1, S1-S4), figure2 (Figure 2, S5), figure3 (Figure 3, S6-S7), figure4 (Figure 4, S8-S9).
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
entry<-if(length(arg))dirname(normalizePath(sub('^--file=','',arg[1]),mustWork=TRUE)) else getwd()
source(file.path(entry,'R/core/00_project.R'),encoding='UTF-8')
opt<-ktb_module_options(commandArgs(TRUE));root<-ktb_find_root(opt$root)
source(file.path(root,'config/analysis_config.R'),encoding='UTF-8')
source(file.path(root,'config/figure3_config.R'),encoding='UTF-8')
source(file.path(root,'config/figure4_config.R'),encoding='UTF-8')
for(f in c('R/core/01_numeric.R','R/core/02_data.R','R/core/03_graphics.R','R/core/05_portable_exports.R',
  'R/modules/figure1_models.R','R/modules/figure1_graphics.R',
  'R/modules/figure2_models.R','R/modules/figure2_graphics.R',
  'R/modules/figure3_models.R','R/modules/figure3_graphics.R',
  'R/modules/figure4_models.R','R/modules/figure4_graphics.R'))source(file.path(root,f),encoding='UTF-8')
options(ktb.pdf_device=opt$`pdf-device`)
runid<-if(is.null(opt$`run-id`))paste0(opt$module,'_',format(Sys.time(),'%Y%m%d_%H%M%S')) else opt$`run-id`
if(!grepl('^[A-Za-z0-9_.-]+$',runid))stop('Use a simple run-id (letters, numbers, ".", "_" or "-").')
outbase<-if(is.null(opt$out))file.path(root,'outputs',runid) else opt$out
if(dir.exists(outbase))stop('Output folder already exists; choose a new --run-id: ',outbase)
dir.create(outbase,recursive=TRUE);outbase<-normalizePath(outbase,mustWork=TRUE)
logcon<-file(file.path(outbase,'execution.log'),'wt');sink(logcon,split=TRUE);sink(logcon,type='message')
failed<-FALSE;warnings<-list()
tryCatch(withCallingHandlers({
  cat('KTB140 analysis\nModule:',opt$module,'\nPDF device:',opt$`pdf-device`,'\nRoot:',root,'\n')
  ktb_require();ktb_portable_require();font<-ktb_check_font()
  packages<-unique(c(ktb_required(),ktb_portable_packages()))
  write_csv(data.frame(package=packages,version=vapply(packages,function(p)as.character(packageVersion(p)),character(1))),file.path(outbase,'package_versions.csv'))
  writeLines(capture.output(sessionInfo()),file.path(outbase,'sessionInfo.txt'))
  registry<-list(figure1=list(recompute=recompute_figure1,render=render_figure1),
    figure2=list(recompute=recompute_figure2,render=render_figure2),
    figure3=list(recompute=recompute_figure3,render=render_figure3),
    figure4=list(recompute=recompute_figure4,render=render_figure4))
  module<-registry[[opt$module]]
  dir.create(file.path(outbase,'checks'),showWarnings=FALSE);dir.create(file.path(outbase,'plot_data'),showWarnings=FALSE)
  result<-module$recompute(root,outbase)
  module$render(root,outbase,result$plots,'recompute')
  if(opt$powerpoint=='yes')ktb_export_powerpoint(root,outbase)
  cat('Run complete. Results, plot data, figures and panels are in:',outbase,'\n')
},warning=function(w) {
  warnings[[length(warnings)+1L]]<<-data.frame(message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=' '))
}),error=function(e) {failed<<-TRUE;message('ERROR: ',conditionMessage(e));writeLines(conditionMessage(e),file.path(outbase,'ERROR.txt'))})
write_csv(if(length(warnings))do.call(rbind,warnings) else data.frame(message=character(),call=character()),file.path(outbase,'warnings.csv'))
ktb_write_status(file.path(outbase,'RUN_STATUS.json'),if(failed)'failed' else 'completed',
  list(module=opt$module,R_version=R.version.string,OS=.Platform$OS.type,pdf_device=opt$`pdf-device`,powerpoint=opt$powerpoint))
sink(type='message');sink();close(logcon)
cat('Output:',outbase,'\n')
# quit() only when run with Rscript, so running this file inside RStudio does not end the session.
if(!interactive())quit(save='no',status=if(failed)1L else 0L)
