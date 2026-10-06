# Native R drawing at the approved physical sizes. Panel d is redrawn at 122 mm,
# never scaled down from 183 mm. The tested portable device/PPTX adapter is reused.
render_figure3 <- function(root,out,bundle,mode) {
  C<-ktb_style;epochs<-figure3_config$spline
  get<-function(n) {x<-bundle[[n]];if(is.null(x))stop('Missing Figure 3 plotting source: ',n);x}
  qp<-function(x)ifelse(x<.001,formatC(x,format='e',digits=1),sprintf('%.3f',x))
  one<-figure3_one
  yspec<-function(ep)switch(ep,er_pooled=list(limits=c(0,60),ticks=c(0,20,40,60),mult=100,label='ER+ epithelial cells (%)'),
    erp63_mean=list(limits=c(-.05,.40),ticks=seq(0,.4,.1),mult=1,label='ER–p63 log-enrichment'),
    DSP_ecm=list(limits=c(-1.5,1.6),ticks=-1:1,mult=1,label='DSP ECM/fibroblast\n(mean-z signature units)'))
  checkbounds<-function(v,limits,label)if(any(!is.finite(v))||any(v<limits[1]-1e-10|v>limits[2]+1e-10))stop(label,' has values outside approved plotting range; do not clip.')
  points<-function(d,x,y,diam_mm=.95,col=C$teal,alpha=.62) {
    d$.x<-d[[x]];d$.y<-d[[y]];d$.diam<-rep(diam_mm,length.out=nrow(d))/.75
    ggplot2::ggplot(d,ggplot2::aes(x=.x,y=.y))+ggplot2::geom_point(ggplot2::aes(size=.diam),colour=col,alpha=alpha,stroke=0)+ggplot2::scale_size_identity()
  }
  smalllegend<-function(x,y,labels,cols,ltys=rep(1,length(labels)),shapes=NULL) {
    for(i in seq_along(labels)) {
      xx<-x+(i-1)*26
      if(is.null(shapes))ktb_line(xx,y,xx+6,y,cols[i],if(i==1)1.2 else .85,ltys[i]) else ktb_symbol(xx+3,y,shapes[i],cols[i],if(i==1)cols[i] else 'white',1.59)
      ktb_text(xx+8,y,labels[i])
    }
  }
  trend<-function(name,x,y) {
    d<-get(paste0(name,'_donors'));r<-get(paste0(name,'_annotation'));g<-get(paste0(name,'_curve'))
    weighted<-name!='Fig3a';outcome<-if(weighted)'ER_low_burden' else 'er_pooled'
    d$.value<-100*d[[outcome]];g$fit<-100*g$fit;g$lower<-100*g$lower;g$upper<-100*g$upper
    yl<-if(weighted)c(-20,107) else c(0,60);ax<-x+14;ay<-y+19;aw<-42;ah<-43
    # Marker areas are 3*n points^2 (burdens) or 8 points^2.
    # Literal diameter in mm is sqrt(area)*25.4/72, then the shared point factor.
    size<-if(weighted)sqrt(3*d$eligible_n)*25.4/72 else sqrt(8)*25.4/72
    checkbounds(d$.value,yl,name);checkbounds(c(g$lower,g$upper),yl,paste(name,'interval'))
    p<-ggplot2::ggplot()+ggplot2::geom_ribbon(data=g,ggplot2::aes(x=x,ymin=lower,ymax=upper),fill=C$teal,alpha=.13)+
      ggplot2::geom_line(data=g,ggplot2::aes(x=x,y=fit),colour=C$teal,linewidth=1.2*25.4/72)
    d$.size<-size/.75
    if(weighted)p<-p+ggplot2::geom_hline(yintercept=0,colour=C$grid,linewidth=.5*25.4/72)
    p<-p+ggplot2::geom_point(data=d,ggplot2::aes(x=TSLB_years,y=.value,size=.size),colour=C$teal,alpha=.62,stroke=0)+ggplot2::scale_size_identity()
    ktb_marks(p,ax,ay,aw,ah,c(1.7,10.3),yl)
    ktb_axes(ax,ay,aw,ah,c(1.7,10.3),yl,xticks=seq(2,10,2),yticks=if(weighted)c(0,25,50,75,100) else c(0,20,40,60),
      xlabel='TSLB (years)',ylabel=if(weighted)'ER-low lobules (%)' else 'ER+ epithelial cells (%)',xlab_offset=6.5,ylab_offset=if(weighted)9.3 else 7.6)
    ktb_text(x+31,y+2,paste0('β = ',ktb_num(r$beta,2,TRUE),' pp/year'),just='center')
    ktb_text(x+31,y+6.6,paste0('95% CI, ',ktb_num(r$lower,2),' to ',ktb_num(r$upper,2)),just='center')
    ktb_text(x+31,y+11.2,paste0('P = ',qp(r$p_raw),'; Holm P = ',qp(r$p_holm)),just='center')
  }
  summary<-function(x,y) {
    d<-get('Fig3d_summary');stopifnot(nrow(d)==6L)
    ktb_title('d',x+2,y+5.5,'Summary of postpartum associations')
    ktb_text(x+100,y+12.5,'Donors',bold=TRUE,just='center');ktb_text(x+115,y+12.5,'Holm P',bold=TRUE,just='center')
    settings<-list(list(top=18,h=22,xl=c(-10,4),ticks=c(-10,-5,0,4),lab='Percentage points/year',ix=1:3),
      list(top=53,h=6,xl=c(-.02,.02),ticks=c(-.02,0,.02),lab='Log-enrichment units/year',ix=4),
      list(top=73,h=6,xl=c(-.12,.12),ticks=c(-.1,0,.1),lab='Mean-z signature units/year',ix=5),
      list(top=93,h=6,xl=c(-.1,.06),ticks=c(-.1,-.05,0,.05),lab='Change in enrichment–fraction\nslope/year',ix=6))
    labs<-c('ER-positive fraction','ER-low: all lobules','ER-low:\nhigh-cell-count lobules','ER–p63 neighborhood\nenrichment','DSP ECM/fibroblast\nsignature','Within-donor p63\ncoupling × TSLB')
    for(s in settings) {
      ax<-x+49;ay<-y+s$top;aw<-43;ah<-s$h;yl<-c(-.5,length(s$ix)-.5)
      ktb_axes(ax,ay,aw,ah,s$xl,yl,xticks=s$ticks,yleft=FALSE,xlabel=s$lab,xlab_offset=if(length(s$ix)==1L&&s$ix==6L)8.4 else 6.5)
      ktb_zero(ax,ay,aw,ah,s$xl)
      for(j in seq_along(s$ix)) {
        i<-s$ix[j];yy<-length(s$ix)-j;Y<-ay+ah*(yl[2]-yy)/diff(yl);sig<-d$p_holm[i]<.05;col<-if(sig)C$teal else C$gray
        ktb_ci(d$beta[i],d$lower[i],d$upper[i],yy,ax,ay,aw,ah,s$xl,yl,col,if(sig)16 else 21,if(sig)col else 'white',4.4*25.4/72)
        ktb_text(x+2,Y,labs[i]);ktb_text(x+100,Y,as.character(d$donors[i]),just='center');ktb_text(x+115,Y,qp(d$p_holm[i]),just='center')
      }
    }
  }
  main<-function() {
    for(i in 1:3) {
      x<-c(0,62,124)[i];name<-c('Fig3a','Fig3b','Fig3c')[i];n<-nrow(get(paste0(name,'_donors')))
      ktb_title(letters[i],x+2,6,if(i==1)'ER-positive fraction' else 'ER-low burden')
      ktb_text(x+31,11.5,paste0(c('Count-pooled','All lobules','High cell count')[i],'; n = ',n),just='center',col=C$gray)
      trend(name,x,15)
    };summary(30.5,96)
  }
  spline<-function(ep,x,y) {
    d<-get(paste0('Spline_',ep,'_donors'));g<-get(paste0('Spline_',ep,'_curves'));r<-one(get('Spline_annotations'),'outcome',ep);sp<-yspec(ep)
    d$.value<-d[[ep]]*sp$mult;g$fit<-g$fit*sp$mult;g$linear<-g$linear*sp$mult
    checkbounds(d$.value,sp$limits,paste('Spline',ep));checkbounds(c(g$fit,g$linear),sp$limits,'Spline/linear curve')
    p<-points(d,'TSLB_years','.value',3*25.4/72,C$gray,.55)+
      ggplot2::geom_line(data=g,ggplot2::aes(x=x,y=fit),inherit.aes=FALSE,colour=C$teal,linewidth=1.25*25.4/72)+
      ggplot2::geom_line(data=g,ggplot2::aes(x=x,y=linear),inherit.aes=FALSE,colour=C$ink,linewidth=.85*25.4/72,linetype='dashed')
    ktb_marks(p,x+20,y+22,65,47,c(1.7,10.3),sp$limits)
    ktb_axes(x+20,y+22,65,47,c(1.7,10.3),sp$limits,xticks=seq(2,10,2),yticks=sp$ticks,
      xlabel='TSLB (years)',ylabel=sp$label,xlab_offset=6.5,ylab_offset=switch(ep,er_pooled=7.8,erp63_mean=8.7,DSP_ecm=9.6))
    ktb_text(x+53,y+2,paste0('n = ',nrow(d),'; nonlinear P = ',qp(r$p_raw)),just='center')
    ktb_text(x+53,y+6.5,paste0('Holm P = ',qp(r$p_holm),'; overall P = ',qp(r$overall_spline_p)),just='center')
    smalllegend(x+28,y+14,c('Spline','Linear'),c(C$teal,C$ink),c(1,2))
  }
  nonlinear<-function(x,y) {
    d<-get('Spline_annotations');d<-d[match(epochs,d$outcome),];ax<-x+28;ay<-y+22;aw<-54;ah<-47;xl<-c(-.02,1.03);yl<-c(-.5,2.5)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(0,.25,.5,.75,1),yleft=FALSE,xlabel='Nonlinearity P value',xlab_offset=6.5)
    ktb_line(ax+aw*(.05-xl[1])/diff(xl),ay,ax+aw*(.05-xl[1])/diff(xl),ay+ah,C$gray,.65,2)
    labs<-c('ER-positive\nfraction','ER–p63\nenrichment','DSP ECM/\nfibroblast')
    for(i in 1:3) {
      yy<-3-i;Y<-ay+ah*(yl[2]-yy)/diff(yl);ktb_text(ax-4,Y,labs[i],just='right')
      ktb_symbol(ax+aw*(d$p_raw[i]-xl[1])/diff(xl),ay+ah*(yl[2]-yy-.11)/diff(yl),16,C$gray,C$gray,4.5*25.4/72)
      ktb_symbol(ax+aw*(d$p_holm[i]-xl[1])/diff(xl),ay+ah*(yl[2]-yy+.11)/diff(yl),22,C$teal,'white',4.5*25.4/72)
    }
    ktb_text(x+53,y+3,'Three-endpoint nonlinearity family',just='center',col=C$gray)
    smalllegend(x+29,y+14,c('Raw P','Holm P'),c(C$gray,C$teal),shapes=c(16,22))
  }
  supp6<-function() {
    for(i in 1:3) {x<-c(0,94,0)[i];y<-c(0,0,98)[i];ktb_title(letters[i],x+2,y+6,c('ER-positive fraction','ER–p63 enrichment','DSP ECM/fibroblast')[i]);spline(epochs[i],x,y+10)}
    ktb_title('d',96,104,'Nonlinearity tests');nonlinear(94,108)
  }
  # Same IQR/whisker definitions as the shared parity boxplots, with early/late
  # group colors and saved donor-keyed horizontal jitter rather than random draws.
  earlybox<-function(ep,panel,x,y) {
    d<-get(paste0('Earlylate_',ep,'_donors'));r<-one(get('Earlylate_annotations'),'outcome',ep);sp<-yspec(ep)
    jit<-read_csv(file.path(root,'config/figure3_visual_jitter.csv'));jit<-jit[jit$panel==panel,]
    d$.off<-jit$offset[match(d$Barcode,jit$Barcode)];stopifnot(!anyNA(d$.off),!anyDuplicated(jit$Barcode))
    d$.v<-d[[ep]]*sp$mult;d$.g<-d$late;d$.c<-ifelse(d$.g==0,C$gray,C$teal)
    st<-do.call(rbind,lapply(0:1,function(k)data.frame(g=k,t(ktb_box_stats(d$.v[d$.g==k])),col=if(k==0)C$gray else C$teal)))
    p<-ggplot2::ggplot()+ggplot2::geom_segment(data=st,ggplot2::aes(x=g,xend=g,y=lo,yend=hi),colour=C$ink,linewidth=.65*25.4/72)+
      ggplot2::geom_segment(data=st,ggplot2::aes(x=g-.125,xend=g+.125,y=lo,yend=lo),colour=C$ink,linewidth=.65*25.4/72)+
      ggplot2::geom_segment(data=st,ggplot2::aes(x=g-.125,xend=g+.125,y=hi,yend=hi),colour=C$ink,linewidth=.65*25.4/72)+
      ggplot2::geom_rect(data=st,ggplot2::aes(xmin=g-.25,xmax=g+.25,ymin=q1,ymax=q3,colour=col,fill=col),alpha=.08,linewidth=.8*25.4/72)+
      ggplot2::geom_segment(data=st,ggplot2::aes(x=g-.25,xend=g+.25,y=med,yend=med,colour=col),linewidth=1.25*25.4/72)+
      ggplot2::geom_point(data=d,ggplot2::aes(x=.g+.off,y=.v,colour=.c),size=3*25.4/72/.75,alpha=.68,stroke=0)+
      ggplot2::scale_colour_identity()+ggplot2::scale_fill_identity()
    checkbounds(d$.v,sp$limits,paste(panel,'donors'));ktb_marks(p,x+15,y+21,40,44,c(-.5,1.5),sp$limits)
    ktb_axes(x+15,y+21,40,44,c(-.5,1.5),sp$limits,xticks=0:1,xticklabels=c('',''),yticks=sp$ticks,
      ylabel=sp$label,ylab_offset=switch(ep,er_pooled=7.6,erp63_mean=8.5,DSP_ecm=9.4),xlabel='Time since last birth',xlab_offset=10)
    ns<-table(factor(d$late,levels=0:1));stopifnot(identical(as.integer(ns),if(ep=='erp63_mean')c(15L,29L) else c(17L,31L)))
    for(i in 1:2)ktb_text(x+15+40*(i-.5)/2,y+69.5,paste0(c('≤3 years','>3 years')[i],'\nn = ',ns[i]),just='center')
    digs<-if(ep=='er_pooled')2 else 3;unit<-if(ep=='er_pooled')' pp' else ''
    ktb_text(x+31,y+3,paste0('Late − early = ',ktb_num(r$beta,digs,TRUE),unit),just='center')
    ktb_text(x+31,y+7.6,paste0('95% CI, ',ktb_num(r$lower,digs),' to ',ktb_num(r$upper,digs)),just='center')
    ktb_text(x+31,y+12.2,paste0('P = ',qp(r$p_raw),'; Holm P = ',qp(r$p_holm)),just='center')
  }
  interaction<-function(ep,x,y) {
    d<-get(paste0('Interaction_',ep,'_lobules'));g<-get(paste0('Interaction_',ep,'_curves'));r<-one(get('Interaction_annotations'),'outcome',ep)
    d$.value<-100*d[[ep]];g$fit<-100*g$fit;xl<-c(-.21,.285);yl<-c(0,65)
    checkbounds(d$within_enrichment,xl,'Interaction x');checkbounds(d$.value,yl,'Interaction y');checkbounds(g$fit,yl,'Interaction line')
    p<-points(d,'within_enrichment','.value',sqrt(7.5)*25.4/72,C$gray,.25)+
      ggplot2::geom_line(data=g[g$TSLB_years==3,],ggplot2::aes(x=x,y=fit),inherit.aes=FALSE,colour=C$gray,linewidth=1.2*25.4/72,linetype='dashed')+
      ggplot2::geom_line(data=g[g$TSLB_years==7,],ggplot2::aes(x=x,y=fit),inherit.aes=FALSE,colour=C$teal,linewidth=1.2*25.4/72)
    ktb_marks(p,x+17,y+26,68,46,xl,yl)
    ktb_axes(x+17,y+26,68,46,xl,yl,xticks=seq(-.2,.2,.1),yticks=c(0,20,40,60),
      xlabel='Within-donor centered\nER–p63 log-enrichment',ylabel=paste0(if(ep=='er_fraction')'ER' else 'p63','+ epithelial cells (%)'),xlab_offset=8.4,ylab_offset=7.8)
    ktb_text(x+51,y+2,paste0(nrow(d),' lobules from ',length(unique(d$Barcode)),' donors'),just='center',col=C$gray)
    ktb_text(x+51,y+6.8,paste0('Interaction β = ',ktb_num(r$beta,3,TRUE),' per year'),just='center')
    ktb_text(x+51,y+11.6,paste0('P = ',qp(r$p_raw),'; Holm P = ',qp(r$p_holm)),just='center')
    smalllegend(x+26,y+19,c('3 years','7 years'),c(C$gray,C$teal),c(2,1))
  }
  supp7<-function() {
    for(i in 1:3) {x<-c(0,62,124)[i];ktb_title(letters[i],x+2,6,c('ER-positive fraction','ER–p63 enrichment','DSP ECM/fibroblast')[i]);earlybox(epochs[i],paste0('FigS7',letters[i]),x,10)}
    ktb_title('d',2,102,'p63–enrichment relationship');ktb_title('e',96,102,'ER–enrichment relationship')
    interaction('p63_fraction',0,106);interaction('er_fraction',94,106)
  }
  .graphics_state$audit<-list()
  dims<-rbind(export_ktb_page(main,'Figure3',183,212,out,crops=list(Fig3a=c(0,0,60,95),Fig3b=c(62,0,60,95),Fig3c=c(124,0,59,95),Fig3d=c(30.5,96,122,114))),
    export_ktb_page(supp6,'FigureS6',183,193,out,crops=list(FigS6a=c(0,0,90,97),FigS6b=c(94,0,89,97),FigS6c=c(0,98,90,95),FigS6d=c(94,98,89,95))),
    export_ktb_page(supp7,'FigureS7',183,200,out,crops=list(FigS7a=c(0,0,60,95),FigS7b=c(62,0,60,95),FigS7c=c(124,0,59,95),FigS7d=c(0,97,90,103),FigS7e=c(94,97,89,103))))
  write_csv(dims,file.path(out,'checks/figure_dimensions.csv'))
  audit<-do.call(rbind,.graphics_state$audit);write_csv(audit,file.path(out,'checks/text_bounds_review.csv'))
  if(any(!audit$approx_in_bounds))warning('Some text bounds require visual inspection; keep approved font sizes and repair layout locally.')
  invisible(dims)
}
