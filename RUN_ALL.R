#!/usr/bin/env Rscript
# Runs the complete analysis: the four figure modules, then the supplementary-table models.
# Each step runs in its own R process; console output for each step is saved in a log.
#
# Usage (from the project root):
#   Rscript RUN_ALL.R [--run-id my_run]
# Outputs: outputs/<run-id>_fig1 ... _fig4, outputs/<run-id>_tables, and outputs/<run-id>_logs.
a<-commandArgs(TRUE);prefix<-paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S'))
if(length(a)){if(length(a)!=2L||a[1]!='--run-id')stop('Usage: Rscript RUN_ALL.R [--run-id my_run]');prefix<-a[2]}
if(!grepl('^[A-Za-z0-9_.-]+$',prefix))stop('Use a simple run ID (letters, numbers, ".", "_" or "-").')
f<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),mustWork=TRUE)) else getwd()
if(!file.exists(file.path(root,'config/analysis_config.R')))stop('Run from the project root.')
setwd(root);Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
rs<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
if(!file.exists(rs))stop('Rscript not found at ',rs)
mods<-c('figure1','figure2','figure3','figure4','tableonly')
runids<-paste0(prefix,'_',c('fig1','fig2','fig3','fig4','tables'))
logs<-file.path('outputs',paste0(prefix,'_logs'))
if(dir.exists(logs)||any(dir.exists(file.path('outputs',runids))))stop('An output folder for this run ID already exists. Choose a new --run-id.')
dir.create(logs,recursive=TRUE)
utils::write.csv(data.frame(module=mods,run=runids,path=file.path('outputs',runids)),file.path(logs,'run_map.csv'),row.names=FALSE)
stepno<-0L
step<-function(script,args=character()){
 stepno<<-stepno+1L;log<-file.path(logs,sprintf('%02d_%s.log',stepno,sub('\\.R$','',basename(script))))
 cat('\nRunning:',script,paste(args,collapse=' '),'\nLog:',log,'\n')
 rc<-system2(rs,args=vapply(c(script,args),shQuote,character(1)),stdout=log,stderr=log)
 cat(paste(tail(readLines(log,warn=FALSE),8),collapse='\n'),'\n')
 if(!identical(as.integer(rc),0L))stop('Step failed; see ',log)
}
for(i in 1:4)step('RUN_MODULE.R',c('--module',mods[i],'--pdf-device','cairo','--powerpoint','yes','--run-id',runids[i]))
step('RUN_TABLE_ONLY.R',c('--run-id',runids[5]))
cat('\nComplete. Output folders are listed in',file.path(logs,'run_map.csv'),'\n')
