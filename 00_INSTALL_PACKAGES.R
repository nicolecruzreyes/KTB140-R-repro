#!/usr/bin/env Rscript
# Checks for, and optionally installs, missing R packages. Installed packages are not updated.
#   Rscript 00_INSTALL_PACKAGES.R          # report missing packages only
#   Rscript 00_INSTALL_PACKAGES.R --yes    # install the missing packages from CRAN
a<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(a))dirname(normalizePath(sub('^--file=','',a[1]))) else getwd()
source(file.path(root,'R/core/00_project.R'),encoding='UTF-8')
p<-unique(c(ktb_required(),'Cairo','pdftools','officer','zip','xml2'))
miss<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
cat('Missing packages:',if(length(miss))paste(miss,collapse=', ') else 'none','\n')
if(length(miss)) {
  if(!'--yes'%in%commandArgs(TRUE))stop('Run again with --yes to install the missing packages.')
  repo<-getOption('repos');if(!length(repo)||any(repo=='@CRAN@'))repo<-c(CRAN='https://cloud.r-project.org')
  install.packages(miss,repos=repo)
}
if(any(!vapply(p,requireNamespace,logical(1),quietly=TRUE)))stop('Some dependencies are still missing.')
cat('All required packages are available.\n')
