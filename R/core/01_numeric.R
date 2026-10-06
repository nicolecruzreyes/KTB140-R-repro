# Shared numerical helpers retained from the revision-2 project.
read_csv <- function(path) {
  if(!file.exists(path))stop('Missing input: ',path,call.=FALSE)
  utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,
                 na.strings=c('','NA','NaN','nan'),encoding='UTF-8')
}
write_csv <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(as.data.frame(x),path,row.names=FALSE,na='',fileEncoding='UTF-8')
}
assert_cols <- function(d,cols,label='data') {
  miss <- setdiff(cols,names(d));if(length(miss))stop(label,' missing: ',paste(miss,collapse=', '))
}
flag <- function(x) !is.na(x)&tolower(as.character(x))%in%c('true','1','1.0','yes')
stdz <- function(x) {
  x<-as.numeric(x); good<-is.finite(x); y<-rep(NA_real_,length(x));
  if(sum(good)<2)return(y)
  s<-stats::sd(x[good]);if(!is.finite(s)||s==0){y[good]<-0;return(y)}
  y[good]<-(x[good]-mean(x[good]))/s;y
}
safe_mean <- function(x) if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
safe_weighted_mean <- function(x,w) {
  use<-is.finite(x)&is.finite(w)&w>0
  if(!any(use))return(NA_real_);sum(x[use]*w[use])/sum(w[use])
}
join_one <- function(left,right,key='Barcode') {
  if(anyDuplicated(right[[key]]))stop('Non-unique right key: ',key)
  i<-match(left[[key]],right[[key]]);new<-setdiff(names(right),names(left))
  cbind(left,right[i,new,drop=FALSE])
}
finite_rows <- function(d,vars) {
  assert_cols(d,vars);ok<-stats::complete.cases(d[,vars,drop=FALSE])
  for(v in vars)if(is.numeric(d[[v]]))ok<-ok&is.finite(d[[v]])
  ok
}
# Return residuals aligned to the original donor rows, including NA for omitted rows.
residualize <- function(d,y,terms) {
  vars<-unique(c(y,terms));ok<-finite_rows(d,vars);out<-rep(NA_real_,nrow(d))
  if(sum(ok)<=length(terms)+2)stop('Insufficient observations for ',y)
  m<-stats::lm(stats::reformulate(terms,response=y),data=d[ok,,drop=FALSE])
  if(m$rank<length(stats::coef(m)))stop('Rank-deficient residualization: ',y)
  out[ok]<-stats::residuals(m);out
}
robust_fit <- function(d,y,terms,term,cluster=FALSE,weight=NULL,scale=1,
                       serum=FALSE,ci_reference=NULL) {
  vars<-unique(c('Barcode',y,terms,weight));ok<-finite_rows(d,vars)
  if(!is.null(weight))ok<-ok&d[[weight]]>0
  dd<-d[ok,,drop=FALSE];if(nrow(dd)<=length(terms)+2)stop('Insufficient rows: ',y)
  f<-stats::reformulate(terms,response=y)
  if(is.null(weight))m<-stats::lm(f,data=dd) else {dd$.w<-dd[[weight]];m<-stats::lm(f,data=dd,weights=.w)}
  if(m$rank!=length(stats::coef(m)))stop('Rank-deficient model: ',deparse(f))
  V<-if(cluster)sandwich::vcovCL(m,cluster=dd$Barcode,type='HC1',cadjust=TRUE) else sandwich::vcovHC(m,type='HC3')
  b<-unname(stats::coef(m)[term]);se<-sqrt(V[term,term]);df<-stats::df.residual(m)
  if(serum){p<-2*stats::pnorm(-abs(b/se));q<-stats::qnorm(.975);ci<-'normal';pr<-'normal'}
  else {p<-2*stats::pt(-abs(b/se),df);ci<-if(is.null(ci_reference))if(cluster)config$cluster_ci else 'residual_t' else ci_reference
        q<-if(ci=='normal')stats::qnorm(.975) else stats::qt(.975,df);pr<-'residual_t'}
  r<-data.frame(outcome=y,term=term,beta=scale*b,se=scale*se,lower=scale*(b-q*se),upper=scale*(b+q*se),
    p_raw=p,n=nrow(dd),donors=length(unique(dd$Barcode)),df=df,
    covariance=if(cluster)'donor-clustered HC1' else 'HC3',p_reference=pr,ci_reference=ci,
    units=if(scale==100)'percentage points' else 'native',stringsAsFactors=FALSE)
  list(row=r,fit=m,V=V,data=dd)
}
mixed_fit <- function(d,y,adjusted=TRUE) {
  vars<-c('Barcode',y,'ParousFlag',if(adjusted)c('Age','BMI'));dd<-d[finite_rows(d,vars),,drop=FALSE]
  f<-stats::as.formula(paste(y,'~ ParousFlag',if(adjusted)'+ Age + BMI' else '', '+ (1|Barcode)'))
  m<-lmerTest::lmer(f,data=dd,REML=FALSE)
  co<-coef(summary(m))['ParousFlag',];b<-co['Estimate'];se<-co['Std. Error'];df<-co['df'];q<-qt(.975,df)
  conv<-m@optinfo$conv$lme4$messages
  data.frame(outcome=y,model=if(adjusted)'ML age/BMI adjusted' else 'ML unadjusted',beta=100*b,se=100*se,
    lower=100*(b-q*se),upper=100*(b+q*se),p_raw=co['Pr(>|t|)'],df=df,n=nrow(dd),
    donors=length(unique(dd$Barcode)),singular=lme4::isSingular(m),
    convergence=if(is.null(conv))'OK' else paste(conv,collapse='; '),units='percentage points',row.names=NULL)
}
adjust_family <- function(d,groups=character(),method='holm',p='p_raw',name='p_holm',n_family=NULL) {
  if(!nrow(d))return(d)
  k<-if(length(groups))do.call(interaction,c(d[,groups,drop=FALSE],list(drop=TRUE,lex.order=TRUE))) else factor(rep('all',nrow(d)))
  d[[name]]<-NA_real_
  for(ix in split(seq_len(nrow(d)),k))d[[name]][ix]<-p.adjust(d[[p]][ix],method=method,n=if(is.null(n_family))length(ix) else n_family)
  d
}
prediction_grid <- function(obj,x,values) {
  dd<-obj$data; vars<-all.vars(stats::delete.response(stats::terms(obj$fit)))
  new<-dd[rep(1,length(values)),,drop=FALSE]
  for(v in vars) {
    if(is.numeric(dd[[v]]))new[[v]]<-rep(mean(dd[[v]]),length(values))
    else new[[v]]<-rep(dd[[v]][1],length(values))
  }
  new[[x]]<-values;X<-model.matrix(delete.response(terms(obj$fit)),new)
  b<-as.numeric(X%*%coef(obj$fit));se<-sqrt(pmax(0,rowSums((X%*%obj$V)*X)))
  crit<-if(obj$row$ci_reference=='normal')qnorm(.975) else qt(.975,obj$row$df)
  data.frame(x=values,fit=b,lower=b-crit*se,upper=b+crit*se)
}
