# Figure 4 / S8-S9 at publication dimensions. Statistical values come from the
# plot data produced by recompute_figure4().
render_figure4 <- function(root,out,bundle,mode) {
  C<-ktb_style;mm<-25.4/72;one<-figure4_one
  get<-function(n){d<-bundle[[n]];if(is.null(d))stop('Missing Figure4 plot data: ',n);d}
  selected<-get('Fig4c_performance')$subset
  if(length(selected)!=1L||selected!='CTAG1B.1')stop('The selected predictor is not CTAG1B.1, but the figure labels are written for CTAG1B.1. Results are saved; review them before drawing this figure.')
  qp<-function(x)ifelse(x<.001,formatC(x,format='e',digits=1),sprintf('%.3f',x))
  bounds<-function(v,lim,label)if(any(!is.finite(v))||any(v<lim[1]-1e-10|v>lim[2]+1e-10))stop(label,' lies outside the approved axis; inspect, do not clip silently.')
  pts<-function(d,x,y,area=10,col=C$teal,alpha=.62){d$.x<-d[[x]];d$.y<-d[[y]];ggplot2::ggplot(d,ggplot2::aes(x=.x,y=.y))+ggplot2::geom_point(colour=col,alpha=alpha,size=sqrt(area)*mm/.75,stroke=0)}
  map<-function(x,y,w,h,xl,yl)list(X=function(z)x+w*(z-xl[1])/diff(xl),Y=function(z)y+h*(yl[2]-z)/diff(yl))
  # ReportLab's approved titles use a baseline; compensate locally in grid.
  title4<-function(letter,x,y,title){ktb_text(x,y-1.77,letter,14,TRUE);ktb_text(x+7,y-1.263,title,10,TRUE)}
  legend<-function(x,y,spread=33.85){ktb_line(x-2,y,x+2,y,C$teal,.85);ktb_symbol(x,y,16,C$teal,C$teal,4.4*mm);ktb_text(x+4.2,y,'Age + BMI + parity');ktb_line(x+spread-2,y,x+spread+2,y,C$gray,.85);ktb_symbol(x+spread,y,22,C$gray,'white',4.2*mm);ktb_text(x+spread+4.2,y,'Age + BMI')}
  # Draw base-10 subscript with ordinary Arial digits. Cairo's Unicode subscript
  # glyphs were missing even though its font inventory reported Arial.
  log_label<-function(x,y){
    parts<-c('−log','10','(raw P)');sizes<-c(9,6.2,9)
    widths<-vapply(seq_along(parts),function(i)grid::convertWidth(grid::grobWidth(grid::textGrob(parts[i],gp=grid::gpar(fontfamily=C$font,fontsize=sizes[i]))),'mm',TRUE),numeric(1))
    start<-y+sum(widths)/2;before<-0
    for(i in seq_along(parts)){ktb_text(x+if(i==2).65 else 0,start-before-widths[i]/2,parts[i],sizes[i],just='center',rot=90);before<-before+widths[i]}
  }
  # Colour encodes signed coefficients, not significance. Text/shape supplies the
  # additional significance cue. Colour values saturate at the declared limits.
  colours<-grDevices::colorRamp(c(C$blue,'#FFFFFF',C$rust))
  colour<-function(v,lim){z<-pmax(0,pmin(1,(v+lim)/(2*lim)));a<-colours(z)/255;grDevices::rgb(a[,1],a[,2],a[,3])}
  colorbar<-function(x,y,w,h,lim,ticks,label=NULL){n<-160L;for(i in 1:n)ktb_rect(x,y+(i-1)*h/n,w,h/n+.02,colour(lim*(1-2*(i-.5)/n),lim),NA,0)
    ktb_rect(x,y,w,h,NA,C$ink,.5)
    for(v in ticks){yy<-y+h*(lim-v)/(2*lim);ktb_line(x+w,yy,x+w+.6,yy,lwd=.5);ktb_text(x+w+1.7,yy,ktb_num(v,if(lim<.1)2 else 1))}
    if(!is.null(label))ktb_text(x+w+12.7,y+h/2,label,9,just='center',rot=90)
  }
  hitplot<-function(x,y){d<-get('Fig4a_counts');d<-d[d$adjustment=='primary',];d<-d[match(figure4_config$states,d$state),];stopifnot(nrow(d)==4L,!anyNA(d$hits))
    ax<-x+37;ay<-y+10;aw<-22;ah<-48;xl<-c(0,2.6);yl<-c(-.55,3.55);m<-map(ax,ay,aw,ah,xl,yl)
    labs<-c('ER–p63\nbalance','ER–p63\nenrichment','Myeloid/antigen-\npresentation','Stromal/\nfibrovascular')
    for(i in 1:4){yy<-4-i;cc<-if(d$state[i]=='myeloid')C$teal else C$gray
      if(d$hits[i]>0)ktb_rect(m$X(0),m$Y(yy+.21),m$X(d$hits[i])-m$X(0),.42*ah/diff(yl),cc,NA,0)
      ktb_text(ax-2,m$Y(yy),labs[i],just='right');ktb_text(m$X(d$hits[i])+1.015,m$Y(yy),as.character(d$hits[i]))
    };ktb_axes(ax,ay,aw,ah,xl,yl,xticks=0:2,yleft=FALSE,xlabel='Significant features',xlab_offset=7)
  }
  association<-function(x,y){d<-get('Fig4b_associations');features<-c('CTAG1B.1','CTSD','FOXP3');ax<-x+18;ay<-y+14;aw<-51;ah<-39;xl<-c(-.025,.59);yl<-c(-.55,2.55);m<-map(ax,ay,aw,ah,xl,yl)
    ktb_zero(ax,ay,aw,ah,xl);ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(0,.2,.4),xticklabels=c('0.0','0.2','0.4'),yleft=FALSE,xlabel='Standardized association coefficient',xlab_offset=7)
    for(i in 1:3){yy<-3-i;ktb_text(ax-2,m$Y(yy),sub('.1','',features[i],fixed=TRUE),just='right')
      for(a in c('primary','sensitivity')){r<-one(d,feature=features[i],adjustment=a);p<-a=='primary';cc<-if(p)C$teal else C$gray
        ktb_ci(r$beta,r$lower,r$upper,yy+if(p).12 else -.12,ax,ay,aw,ah,xl,yl,cc,if(p)16 else 22,if(p)cc else 'white',4.4*mm)
        ktb_text(x+if(p)82 else 104,m$Y(yy),qp(r$q_BH209),just='center')
      }}
    ktb_text(x+82,y+8,'Primary q',bold=TRUE,just='center');ktb_text(x+104,y+8,'Sensitivity q',bold=TRUE,just='center');legend(x+35.25,y+70.7)
  }
  prediction<-function(x,y){d<-get('Fig4c_donor_averages');g<-get('Fig4c_descriptive_line');r<-get('Fig4c_performance');xl<-c(-2.3,3.5);yl<-c(-.7,1.55)
    bounds(d$target,xl,'Cross-fitted tissue targets');bounds(d$pred,yl,'Cross-fitted predictions')
    p<-pts(d,'target','pred',10,alpha=.6)+ggplot2::geom_line(data=g,ggplot2::aes(x=x,y=fit),inherit.aes=FALSE,colour=C$teal,linewidth=1.15*mm)
    ktb_marks(p,x+14,y+20,42,49,xl,yl);ktb_axes(x+14,y+20,42,49,xl,yl,xticks=c(-2,0,2),yticks=seq(-.5,1.5,.5),yticklabels=sprintf('%.1f',seq(-.5,1.5,.5)),xlabel='Adjusted tissue target\n(cross-fitted residual)',ylabel='Predicted tissue residual',ylab_offset=10.8,xlab_offset=8.7)
    ktb_text(x+31,y+2,paste0('n = ',nrow(d),'; r = ',sprintf('%.3f',r$cv_r)),just='center');ktb_text(x+31,y+7,'100 × 10-fold cross-fitting',just='center')
  }
  parity<-function(x,y){d<-get('Fig4d_donors');r<-one(get('Clinical_annotations'),model='parity');yl<-c(-.7,1.55);xl<-c(-.6,1.6);bounds(d$locked_score,yl,'Parity scores')
    jitter<-read_csv(file.path(root,'config/figure4_visual_jitter.csv'));d<-d[order(match(d$Barcode,jitter$Barcode)),]
    p<-ktb_boxplot(d,'locked_score',jitter,point_mm=sqrt(7)*mm)
    # Approved neutral outlines, coloured medians/points; retain donor jitter.
    for(i in 1:4){p$layers[[i]]$mapping$colour<-NULL;p$layers[[i]]$aes_params$colour<-C$ink;p$layers[[i]]$aes_params$linewidth<-(if(i==4).75 else .65)*mm}
    p$layers[[4]]$aes_params$alpha<-.09
    p$layers[[5]]$aes_params$alpha<-.62
    p$layers[[6]]$aes_params$linewidth<-1.3*mm
    ktb_marks(p,x+14,y+20,42,49,xl,yl);ns<-table(factor(d$ParousFlag,levels=0:1))
    ktb_axes(x+14,y+20,42,49,xl,yl,xticks=numeric(),yticks=seq(-.5,1.5,.5),yticklabels=sprintf('%.1f',seq(-.5,1.5,.5)),ylabel='CTAG1B serum-only score',ylab_offset=10.8)
    for(i in 0:1)ktb_text(x+14+42*(i-xl[1])/diff(xl),y+72.8,paste0(c('Nulliparous','Parous')[i+1],'\nn = ',as.integer(ns)[i+1]),just='center')
    ktb_text(x+31,y+2,paste0('Adjusted β = ',ktb_num(r$beta,3,TRUE)),just='center');ktb_text(x+31,y+7,paste0('P = ',qp(r$p_raw)),just='center')
  }
  postpartum<-function(x,y){d<-get('Fig4e_donors');g<-get('Fig4e_adjusted_line_and_HC3_t_interval');r<-one(get('Clinical_annotations'),model='postpartum');xl<-c(1.7,10.3);yl<-c(-.7,1.55)
    bounds(d$locked_score,yl,'Postpartum score');bounds(c(g$lower,g$upper),yl,'Postpartum confidence band')
    p<-ggplot2::ggplot()+ggplot2::geom_ribbon(data=g,ggplot2::aes(x=TSLB_years,ymin=lower,ymax=upper),fill=C$teal,alpha=.13)+ggplot2::geom_line(data=g,ggplot2::aes(x=TSLB_years,y=fit),colour=C$teal,linewidth=1.2*mm)+ggplot2::geom_point(data=d,ggplot2::aes(x=TSLB_years,y=locked_score),colour=C$teal,alpha=.65,size=sqrt(11)*mm/.75,stroke=0)
    ktb_marks(p,x+14,y+20,42,49,xl,yl);ktb_axes(x+14,y+20,42,49,xl,yl,xticks=seq(2,10,2),yticks=seq(-.5,1.5,.5),yticklabels=sprintf('%.1f',seq(-.5,1.5,.5)),xlabel='TSLB (years)',ylabel='CTAG1B serum-only score',ylab_offset=10.8,xlab_offset=7)
    ktb_text(x+31,y+2,paste0('β = ',ktb_num(r$beta,3,TRUE),'/year; n = ',nrow(d)),just='center');ktb_text(x+31,y+7,paste0('95% CI, ',ktb_num(r$lower,3),' to ',ktb_num(r$upper,3)),just='center');ktb_text(x+31,y+12,paste0('P = ',qp(r$p_raw)),just='center')
  }
  main<-function(){title4('a',2,6,'Serum hits by state');ktb_text(35,11.5,'109 donors; BH q < 0.05',just='center',col=C$gray);hitplot(0,14)
    title4('b',71,6,'Serum–tissue associations');ktb_text(124,11.5,'Myeloid/antigen-presentation state',just='center',col=C$gray);association(69,14)
    title4('c',2,98,'Cross-fitted prediction');prediction(0,102);title4('d',64,98,'Score by parity');parity(62,102);title4('e',126,98,'Postpartum association');postpartum(124,102)
  }
  leading<-function(x,y){d<-get('FigS8a_all_state_results_for_selected_features');features<-get('FigS8a_row_order')$feature;states<-figure4_config$states;n<-length(features);ax<-x+20;ay<-y+3;cw<-54/4;ch<-50/n
    for(i in seq_along(features)){ktb_text(ax-1.6,ay+(i-.5)*ch,sub('.1','',features[i],fixed=TRUE),just='right')
      for(j in seq_along(states)){r<-one(d,feature=features[i],state=states[j]);ktb_rect(ax+(j-1)*cw,ay+(i-1)*ch,cw,ch,colour(r$beta,.4),NA,0);ktb_text(ax+(j-.5)*cw,ay+(i-.5)*ch,paste0(sprintf('%.2f',r$beta),if(r$q_BH209<.05)'*' else ''),just='center',col=if(abs(r$beta)>.31)'white' else C$ink)}}
    labs<-c('ER–p63\nbalance','ER–p63\nenrich.','Myeloid','Stromal/\nvascular');for(j in 1:4)ktb_text(ax+(j-.5)*cw,ay+55,labs[j],just='center')
    colorbar(x+76,y+3,2.5,50,.4,c(-.4,0,.4));ktb_text(x+47,y+71,'Standardized coefficient; * q < 0.05',just='center')
  }
  volcano<-function(x,y){d<-get('FigS8b_screen');d$.hit<-d$q_BH209<.05;ax<-x+16;ay<-y+4;aw<-68;ah<-54;xl<-c(-.12,.46);yl<-c(-.1,3.95);m<-map(ax,ay,aw,ah,xl,yl)
    bounds(d$beta,xl,'Volcano coefficients');bounds(d$minuslog10p,yl,'Volcano P values')
    p<-pts(d[!d$.hit,],'beta','minuslog10p',8,C$gray,.55);ktb_marks(p,ax,ay,aw,ah,xl,yl)
    for(i in which(d$.hit))ktb_symbol(m$X(d$beta[i]),m$Y(d$minuslog10p[i]),22,C$teal,C$teal,sqrt(20)*mm,lwd=0)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=seq(-.1,.4,.1),yticks=0:3,xlabel='Standardized association coefficient',ylabel=NULL,xlab_offset=7,ylab_offset=10);log_label(ax-5.5,ay+ah/2)
    offsets<-list(CTAG1B.1=c(7,-12),CTSD=c(-35,10),FOXP3=c(-39,-13))
    for(f in names(offsets)){r<-one(d,feature=f);o<-offsets[[f]]*mm;xx<-m$X(r$beta)+o[1];yy<-m$Y(r$minuslog10p)-o[2];lab<-sub('.1','',f,fixed=TRUE);tw<-grid::convertWidth(grid::grobWidth(grid::textGrob(lab,gp=grid::gpar(fontfamily=C$font,fontsize=8.5))),'mm',TRUE);ex<-if(o[1]<0)xx+tw+.8 else xx-.8;ktb_line(m$X(r$beta),m$Y(r$minuslog10p),ex,yy-.4,C$gray,.45);ktb_text(xx,yy-.8,lab)}
    ktb_symbol(x+20.85,y+73.7,22,C$teal,C$teal,4.3*mm);ktb_text(x+24,y+73.7,'Within-state q < 0.05');ktb_symbol(x+55.65,y+73.7,16,C$gray,C$gray,3.2*mm);ktb_text(x+58.75,y+73.7,'Other features')
  }
  proxy<-function(x,y){d<-get('FigS8c_proxy_donors');g<-get('FigS8c_descriptive_line_interval');r<-get('FigS8c_annotation');xl<-c(-1.9,3.2);yl<-c(-2,3.2)
    bounds(d$proxy_residual,xl,'Proxy residual');bounds(d$tissue_residual,yl,'Proxy tissue target')
    p<-ggplot2::ggplot()+ggplot2::geom_ribbon(data=g,ggplot2::aes(x=x,ymin=lower,ymax=upper),fill=C$teal,alpha=.13)+ggplot2::geom_line(data=g,ggplot2::aes(x=x,y=fit),colour=C$teal,linewidth=1.15*mm)+ggplot2::geom_point(data=d,ggplot2::aes(x=proxy_residual,y=tissue_residual),colour=C$teal,alpha=.6,size=sqrt(10)*mm/.75,stroke=0)
    ktb_marks(p,x+16,y+12,66,42,xl,yl);ktb_axes(x+16,y+12,66,42,xl,yl,xticks=-1:3,yticks=c(-2,0,2),xlabel='CTAG1B/CTSD proxy residual',ylabel='Tissue-state residual',xlab_offset=7,ylab_offset=10)
    ktb_text(x+47,y+3,paste0('n = ',nrow(d),'; adjusted r = ',sprintf('%.3f',r$adjusted_in_sample_r)),just='center')
  }
  candidate<-function(x,y){d<-get('FigS8d_error_summary');d<-d[order(d$mean_cv_mse,d$n_features),];ax<-x+33;ay<-y+12;aw<-41;ah<-40;xl<-c(.74,.88);yl<-c(-.55,2.55);m<-map(ax,ay,aw,ah,xl,yl)
    stopifnot(nrow(d)==3L);ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(.75,.80,.85),xticklabels=c('0.75','0.80','0.85'),yleft=FALSE,xlabel='Held-out mean squared error',xlab_offset=7)
    for(i in 1:3){yy<-3-i;cc<-if(i==1)C$teal else C$gray;ktb_ci(d$mean_cv_mse[i],d$repeat_mse_q025[i],d$repeat_mse_q975[i],yy,ax,ay,aw,ah,xl,yl,cc,if(i==1)16 else 21,if(i==1)cc else 'white',4.4*mm)
      lab<-gsub('.1','',d$subset[i],fixed=TRUE);lab<-sub(' + ','\n+ ',lab,fixed=TRUE);ktb_text(ax-2,m$Y(yy),lab,just='right');ktb_text(x+85,m$Y(yy),sprintf('%.3f',d$mean_cv_mse[i]),just='center')}
    ktb_text(x+85,y+5,'Mean',bold=TRUE,just='center');ktb_text(x+48,y+65,'50 repeats × 10 folds',just='center')
  }
  component<-function(x,y){d<-get('FigS8e_component_tests');ref<-get('Fig4b_associations');ax<-x+63;ay<-y+13;aw<-62;ah<-32;xl<-c(-.035,.64);yl<-c(-.6,2.6);m<-map(ax,ay,aw,ah,xl,yl)
    ktb_zero(ax,ay,aw,ah,xl);ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(0,.2,.4,.6),xticklabels=c('0.0','0.2','0.4','0.6'),yleft=FALSE,xlabel='Standardized CTAG1B association',xlab_offset=7)
    labs<-c('Integrated state\n(reference)','DSP immune/myeloid\nblock','IO360 myeloid/antigen-\npresentation')
    for(i in 1:3){yy<-3-i;ktb_text(ax-3,m$Y(yy),labs[i],just='right');for(a in c('primary','sensitivity')){
      r<-if(i==1)one(ref,feature='CTAG1B.1',adjustment=a) else one(d,component=c('DSP_block','RNA_Myeloid')[i-1],adjustment=a);p<-a=='primary';cc<-if(p)C$teal else C$gray
      ktb_ci(r$beta,r$lower,r$upper,yy+if(p).12 else -.12,ax,ay,aw,ah,xl,yl,cc,if(p)16 else 22,if(p)cc else 'white',4.4*mm)
      ktb_text(x+if(p)142 else 170,m$Y(yy),if(i==1)'—' else qp(r$p_holm),just='center')}}
    ktb_text(x+142,y+6,'Primary',bold=TRUE,just='center');ktb_text(x+170,y+6,'Sensitivity',bold=TRUE,just='center');ktb_text(x+156,y+10,'Holm P',just='center');legend(x+70.8,y+58.7)
  }
  s8<-function(){title4('a',2,6,'Leading serum associations');leading(0,12);title4('b',94,6,'Myeloid-state screen');volcano(92,12);title4('c',2,98,'Two-feature proxy');proxy(0,104);title4('d',94,98,'Candidate-model comparison');candidate(92,104);title4('e',2,182,'CTAG1B association across tissue components');component(0,186)}
  eps<-figure4_config$endpoints;elabs<-c('Lobule p63 fraction','Lobule ER fraction','Lobule ER–p63 enrichment','All lobules: p63-high burden','All lobules: ER-low burden','All lobules: ER–p63-enriched burden','All lobules: composite burden','High count: p63-high burden','High count: ER-low burden','High count: ER–p63-enriched burden','High count: composite burden','Within-donor ER coupling','Within-donor p63 coupling')
  exploratory_hits<-function(x,y){d<-get('FigS9a_counts');d<-d[match(eps,d$endpoint),];ax<-x+73;ay<-y+13;aw<-71;ah<-76;xl<-c(-.25,4.5);yl<-c(-.6,12.6);m<-map(ax,ay,aw,ah,xl,yl)
    bounds(c(d$nominal,d$BH),xl,'Exploratory hit counts');ktb_axes(ax,ay,aw,ah,xl,yl,xticks=0:4,yleft=FALSE,xlabel='Number of serum features',xlab_offset=7)
    for(v in c(9.5,5.5,1.5))ktb_line(ax,m$Y(v),ax+aw,m$Y(v),C$grid,.55)
    for(i in seq_along(eps)){yy<-13-i;ktb_text(ax-2,m$Y(yy),elabs[i],just='right');ktb_symbol(m$X(d$nominal[i]),m$Y(yy+.11),16,C$gray,C$gray,4*mm);ktb_symbol(m$X(d$BH[i]),m$Y(yy-.11),22,C$teal,'white',4.3*mm);ktb_text(x+171,m$Y(yy),as.character(d$donors[i]),just='center')}
    ktb_text(x+171,y+7,'Donors',bold=TRUE,just='center');ktb_symbol(x+58.4,y+2.7,16,C$gray,C$gray,4*mm);ktb_text(x+61.65,y+2.7,'Nominal P < 0.01');ktb_symbol(x+92,y+2.7,22,C$teal,'white',4.3*mm);ktb_text(x+95.6,y+2.7,'Within-endpoint BH q < 0.05');ktb_text(x+95,y+103,'209 features per endpoint; 2,717 tests overall',just='center',col=C$gray)
  }
  exploratory_matrix<-function(x,y){d<-get('FigS9b_candidate_coefficients');features<-c('FOXP3','CTSD','CTAG1B.1');ax<-x+73;ay<-y+11;cw<-65/3;ch<-76/13
    for(i in seq_along(eps)){ktb_text(ax-2,ay+(i-.5)*ch,elabs[i],just='right');for(j in 1:3){r<-one(d,endpoint=eps[i],feature=features[j]);bounds(r$beta,c(-.08,.08),'Exploratory coefficient colour range');ktb_rect(ax+(j-1)*cw,ay+(i-1)*ch,cw,ch,colour(r$beta,.08),NA,0);ktb_text(ax+(j-.5)*cw,ay+(i-.5)*ch,sprintf('%+.3f',r$beta),just='center',col=if(abs(r$beta)>.065)'white' else C$ink)}}
    for(i in c(3,7,11))ktb_line(ax,ay+i*ch,ax+65,ay+i*ch,C$line,.65)
    for(j in 1:3)ktb_text(ax+(j-.5)*cw,ay+80,sub('.1','',features[j],fixed=TRUE),just='center')
    colorbar(x+146,y+11,3,76,.08,c(-.08,-.04,0,.04,.08),'Coefficient (native outcome units)')
    ktb_text(x+96,y+101,'Serum features standardized; outcomes retain endpoint-specific units.',just='center');ktb_text(x+96,y+106,'No displayed association met within-endpoint BH q < 0.05.',just='center',col=C$gray)
  }
  s9<-function(){title4('a',2,6,'Serum screen across lobule-aware endpoints');exploratory_hits(0,12);title4('b',2,126,'Candidate serum features and epithelial endpoints');exploratory_matrix(0,130)}
  export_ktb_page(main,'Figure4',183,190,out,list(Fig4a=c(0,0,67,93),Fig4b=c(69,0,114,93),Fig4c=c(0,94,60,96),Fig4d=c(62,94,60,96),Fig4e=c(124,94,59,96)))
  export_ktb_page(s8,'FigureS8',183,248,out,list(FigS8a=c(0,0,90,93),FigS8b=c(92,0,91,93),FigS8c=c(0,94,90,84),FigS8d=c(92,94,91,84),FigS8e=c(0,178,183,70)))
  export_ktb_page(s9,'FigureS9',183,240,out,list(FigS9a=c(0,0,183,121),FigS9b=c(0,122,183,118)))
}
