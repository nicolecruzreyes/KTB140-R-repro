# Figure 1 / S1-S4 layouts.
# All five pages are composed at 183 mm width. No approved PDF is used as artwork.

render_figure1 <- function(root,out,bundle,mode) {
  C<-ktb_style;dget<-function(n){if(is.null(bundle[[n]]))stop('Missing plotting table: ',n);bundle[[n]]}
  cells<-read_csv(file.path(root,'data_cells/PairB_cells.csv'));cells$complete_calls<-flag(cells$complete_calls)
  lookup<-read_csv(file.path(root,'data_cells/PairB_lobule_lookup.csv'))
  jitter<-read_csv(file.path(root,'config/visual_jitter.csv'))
  one<-function(q,key,val){i<-which(q[[key]]==val);if(length(i)!=1L)stop('Nonunique/missing row ',key,'=',val);q[i,,drop=FALSE]}
  U<-lookup$lobule_uid[match(c('Nulliparous','Parous'),lookup$group)]
  if(anyNA(U))stop('Missing selected Pair B fields')
  S<-cells[cells$lobule_uid==U[2]&cells$complete_calls,,drop=FALSE]
  index<-which(S$cell_id==22718);if(length(index)!=1L)stop('Missing fixed neighborhood example cell 22718')
  distance<-sqrt((S$x_um-S$x_um[index])^2+(S$y_um-S$y_um[index])^2)
  nb<-distance<=15&seq_len(nrow(S))!=index
  stopifnot(sum(nb)==9L,sum(S$p63_call[nb]==1)==2L)

  cellmap<-function(uid,x,y,w,h,ki67=FALSE) {
    q<-cells[cells$lobule_uid==uid,,drop=FALSE];center<-c(mean(range(q$x_um)),mean(range(q$y_um)))
    H<-670;W<-H*w/h;xx<-q$x_um-center[1];yy<-q$y_um-center[2]-14
    if(any(abs(xx)>=W/2)|any(abs(yy)>=H/2))stop('Cell falls outside fixed map extent')
    X<-x+w*(xx+W/2)/W;Y<-y+h*(yy+H/2)/H
    e<-q$ER_call==1;p<-q$p63_call==1;k<-q$Ki67_call==1;ok<-q$complete_calls
    mark<-function(sel,shape=16,col=C$other,fill=col,size=.4,lwd=.25){sel<-which(!is.na(sel)&sel);if(length(sel))ktb_symbol(X[sel],Y[sel],shape,col,fill,size,lwd)}
    sz<-(if(w<60)1.03 else 1.45)*25.4/72
    mark(ok,16,ktb_alpha(C$other,.82),size=sz)
    mark(!ok,4,C$gray,size=sz,lwd=.25)
    mark(e&!p,16,C$rust,size=sz);mark(p&!e,16,C$teal,size=sz)
    mark(e&p,21,C$rust,C$teal,size=sz,lwd=if(w<60).25 else .33)
    if(ki67)mark(k,1,C$ink,NA,sz+1.35*25.4/72,.42)
    bx<-x+90*w/W;by<-y+h*(H-21)/H
    ktb_line(bx,by,bx+100*w/W,by,C$ink,1)
    ktb_text(bx+50*w/W,by-2,'100 µm',8,just='center')
  }
  mapkey<-function(y,ki67=FALSE){
    xx<-if(ki67)c(10,38,67,96,124,155) else c(25,59,93,133,164)
    labs<-if(ki67)c('ER+','p63+','Both','Other','Uncalled','Ki67+') else c('ER+ / p63−','p63+ / ER−','Double-positive','Other called','Uncalled')
    cols<-c(C$rust,C$teal,C$rust,C$other,C$gray,C$ink)
    shapes<-c(16,16,21,16,4,1);fills<-c(C$rust,C$teal,C$teal,C$other,NA,NA)
    for(i in seq_along(xx)){ktb_symbol(xx[i],y,shapes[i],cols[i],fills[i],1.4);ktb_text(xx[i]+2,y,labs[i],if(ki67)8.5 else 8)}
  }
  forest<-function(q,labels,x,y,w,h,xlim,xlabel,pcol='p_holm',left=30,xticks=NULL,scale=1) {
    right<-if(is.null(pcol))3 else 17;ax<-x+left;ay<-y+7;aw<-w-left-right;ah<-h-22
    yl<-c(-.5,nrow(q)-.5);yy<-rev(seq_len(nrow(q))-1)
    ktb_axes(ax,ay,aw,ah,xlim,yl,xticks=if(is.null(xticks))pretty(xlim,n=3) else xticks,
      yleft=FALSE,xlabel=xlabel,xlab_offset=11)
    ktb_zero(ax,ay,aw,ah,xlim)
    for(i in seq_len(nrow(q))){r<-q[i,];Y<-ay+ah*(yl[2]-yy[i])/diff(yl)
      ktb_text(ax-2,Y,labels[i],just='right')
      ktb_ci(r$beta*scale,r$lower*scale,r$upper*scale,yy[i],ax,ay,aw,ah,xlim,yl,size=1.6)
      if(!is.null(pcol))ktb_text(ax+aw+4,Y,sprintf('%.3f',r[[pcol]]),just='left')}
    if(!is.null(pcol))ktb_text(ax+aw+4,y+2,'Holm P',bold=TRUE)
  }

  main <- function() {
    ktb_title('a',2,5.2,'Representative lobules')
    for(i in seq_along(c(23,76,129)))ktb_text(c(23,76,129)[i]+25.5,11,c('H&E','MxIF composite','Cell classification')[i],bold=TRUE,just='center')
    for(i in 1:2){y<-c(15.5,52.5)[i];ktb_text(11.5,y+17,c('Nulliparous','Parous')[i],bold=TRUE,just='center',col=c(C$blue,C$rust)[i])
      for(x in c(23,76)){ktb_rect(x,y,51,34);ktb_text(x+25.5,y+17,'Image pending',just='center',col=C$gray)}
      cellmap(U[i],129,y,51,34);ktb_rect(129,y,51,34,fill=NA)}
    mapkey(87.9)
    donors<-dget('Fig1bc_donor_values');mixed<-dget('Fig1bc_mixed_model_annotations')
    for(i in 1:2){marker<-c('p63','er')[i];ax<-c(15,111)[i];letter<-c('b','c')[i]
      ktb_title(letter,c(2,96)[i],95,paste0(if(i==1)'p63' else 'ER','-positive epithelium'))
      r<-one(mixed,'outcome',paste0(marker,'_fraction'))
      label<-if(identical(as.character(r$ci_reference),'Satterthwaite'))'P (Satt.)' else 'Wald P'
      ktb_text(ax,102,paste0('Adjusted difference: ',ktb_num(r$beta,1,TRUE),' pp'))
      ktb_text(ax,106.5,paste0('95% CI ',ktb_num(r$lower),' to ',ktb_num(r$upper),'; ',label,' = ',ktb_p(r$p_raw)))
      dd<-donors;dd$.value<-dd[[paste0(marker,'_pooled')]]*100
      if(any(dd$.value<0|dd$.value>70))stop('Composition data outside approved axis')
      p<-ktb_boxplot(dd,'.value',jitter[jitter$panel=='Fig1bc',])
      ktb_marks(p,ax,112,67,39,c(-.5,1.5),c(0,70),ygrid=c(0,20,40,60))
      ktb_axes(ax,112,67,39,c(-.5,1.5),c(0,70),yticks=c(0,20,40,60),
        ylabel=paste0(if(i==1)'p63' else 'ER','+ epithelial cells (%)'),ylab_offset=11)
      for(g in 0:1)ktb_text(ax+67*(g+.5)/2,157,paste0(c('Nulliparous','Parous')[g+1],'\nn = ',sum(dd$ParousFlag==g)),just='center')
    }
    ktb_title('d',2,170,'Neighborhood enrichment')
    q<-dget('Fig1d_neighborhood_models');keys<-c('erp63_enrichment_called','kp63_enrichment_called','erki_enrichment_called')
    ax<-24;ay<-183;aw<-39;ah<-36;yl<-c(-.55,2.55);xl<-c(-.15,.15)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(-.1,0,.1),yleft=FALSE);ktb_zero(ax,ay,aw,ah,xl)
    ktb_text(72,179,'Holm P',bold=TRUE,just='center')
    for(i in 1:3){r<-one(q,'outcome',keys[i]);yy<-3-i;Y<-ay+ah*(yl[2]-yy)/diff(yl)
      ktb_ci(r$beta,r$lower,r$upper,yy,ax,ay,aw,ah,xl,yl,size=2.1)
      ktb_text(21.2,Y-1.65,c('ER–p63','Ki67–p63','ER–Ki67')[i],just='right')
      ktb_text(21.2,Y+1.7,paste0(r$n_lobules_or_donors,' lobules'),just='right',col=C$gray)
      ktb_text(72,Y,sprintf('%.3f',r$p_holm),just='center')}
    ktb_text(43.5,227,'Parity difference\n(log-enrichment units)',9,just='center')
    ktb_title('e',82,170,'Burden of altered lobules')
    ktb_symbol(94,177,21,C$gray,'white',2);ktb_text(98,177,'All eligible')
    ktb_symbol(123,177,22,C$teal,C$teal,2);ktb_text(127,177,'High-cell-count')
    q<-dget('Fig1e_burden_models');ax<-113;ay<-183;aw<-47;ah<-36;xl<-c(-12,36);yl<-c(-.6,3.6)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(-10,0,10,20,30),yleft=FALSE);ktb_zero(ax,ay,aw,ah,xl)
    ktb_text(174,179,'Holm P',bold=TRUE,just='center')
    keys<-c('p63_high_burden','ER_low_burden','ERp63_high_burden','composite_burden')
    for(i in 1:4){Y<-ay+ah*(yl[2]-(4-i))/diff(yl);ktb_text(109,Y,c('p63-high','ER-low','ER–p63-enriched','Composite-altered')[i],just='right')
      for(j in 1:2){r<-q[q$outcome==keys[i]&q$subset==c('all_lobules','count_rich')[j],];stopifnot(nrow(r)==1)
        yy<-4-i+c(.19,-.19)[j];col<-c(C$gray,C$teal)[j]
        ktb_ci(r$beta,r$lower,r$upper,yy,ax,ay,aw,ah,xl,yl,col,c(21,22)[j],c('white',C$teal)[j],1.85)
        ktb_text(174,ay+ah*(yl[2]-yy)/diff(yl),sprintf('%.3f',r$p_holm),just='center')}}
    ktb_text(137,227,'Parity difference\n(percentage points)',9,just='center')
  }

  s1 <- function() {
    ktb_title('a',2,5,'Assay and analysis coverage')
    ktb_rect(108,5,2.2,2.2,C$teal,NA);ktb_text(112,6.1,'Available')
    ktb_rect(139,5,2.2,2.2,C$pale);ktb_text(143,6.1,'Unavailable')
    q<-dget('FigS1a');cols<-setdiff(names(q),'Barcode');stopifnot(length(cols)==9,nrow(q)==140)
    labs<-c('MxIF','GeoMx DSP','IO360','Serum','Plasma','Four tissue states','Four-state serum','Postpartum tissue','Postpartum serum')
    x<-48;y<-14;w<-130;h<-48
    for(i in 1:9){Y<-y+(i-.5)*h/9;ktb_text(x-2,Y,labs[i],just='right')
      for(j in 1:140)ktb_rect(x+(j-1)*w/140,y+(i-1)*h/9,w/140,h/9,if(flag(q[[cols[i]]][j]))C$teal else C$pale,NA)}
    for(j in c(1,35,70,105,140))ktb_text(x+(j-.5)*w/140,y+h+3,as.character(j),just='center')
    ktb_text(x+w/2,y+h+9,'Donors ordered by data availability',9,just='center')
    ktb_title('b',2,83.5,'Clinical and reproductive data availability');ktb_text(181,85.5,'Missing / eligible',just='right',col=C$gray)
    q<-dget('FigS1b');x<-56;y<-92;w<-105;h<-50;N<-nrow(q)
    for(i in seq_len(N)){Y<-y+(i-.5)*h/N;ktb_text(x-2,Y,q$label[i],just='right');ktb_rect(x,Y-h/N*.29,w*q$percent[i]/100,h/N*.58,C$gray,NA);ktb_text(x+w+5,Y,paste0(q$missing[i],'/',q$denominator[i]))}
    ktb_axes(x,y,w,h,c(0,100),c(0,N),xticks=c(0,25,50,75,100),yleft=FALSE,xlabel='Missing among applicable donors (%)')
    ktb_text(55,157,'Age at first birth and TSLB: parous donors only.',col=C$gray)
    ktb_title('c',2,168.5,'Assay coverage and analysis subsets');q<-dget('FigS1c')
    labs<-c('All donors / blood samples','IO360','Primary serum','MxIF','GeoMx DSP','All three tissue modalities','Four corrected tissue states','Common four-state serum','Postpartum tissue composition','Postpartum primary serum')
    x<-69;y<-177;w<-103;h<-53;N<-nrow(q)
    for(i in seq_len(N)){Y<-y+(i-.5)*h/N;ktb_text(x-2,Y,labs[i],just='right');ktb_rect(x,Y-h/N*.29,w*q$n[i]/153,h/N*.58,C$teal,NA);ktb_text(x+w*q$n[i]/153+2,Y,as.character(q$n[i]))}
    ktb_axes(x,y,w,h,c(0,153),c(0,N),xticks=c(0,50,100,150),yleft=FALSE,xlabel='Donors')
  }

  heterogeneity <- function(q,marker,x,y,w=89,h=58) {
    ax<-x+15;ay<-y+10;aw<-w-17;ah<-h-20
    dd<-do.call(rbind,lapply(0:1,function(g){z<-q[q$ParousFlag==g,];z<-z[order(z$mean,z$Barcode),];z$.x<-seq_len(nrow(z))-1+if(g==0)0 else 65;z$.c<-c(C$blue,C$rust)[g+1];z}))
    p<-ggplot2::ggplot(dd)+ggplot2::geom_segment(ggplot2::aes(x=.x,xend=.x,y=min*100,yend=max*100,colour=.c),linewidth=.18,alpha=.45)+
      ggplot2::geom_point(ggplot2::aes(x=.x,y=mean*100,colour=.c),size=.72/.75,stroke=0)+ggplot2::scale_colour_identity()
    ktb_marks(p,ax,ay,aw,ah,c(-3,134),c(0,80),ygrid=c(0,20,40,60,80));ktb_axes(ax,ay,aw,ah,c(-3,134),c(0,80),yticks=c(0,20,40,60,80),
      ylabel=paste0(marker,'+ epithelial cells (%)'),xlabel='Donors ordered by mean',xlab_offset=6,ylab_offset=11)
    for(g in 0:1){z<-dd[dd$ParousFlag==g,];ktb_text(ax+aw*(mean(z$.x)+3)/137,y+3,paste0(c('Nulliparous','Parous')[g+1],'\nn = ',nrow(z)),just='center')}
  }

  s2 <- function() {
    ktb_title('a',2,4.8,'Marker calls in the selected parous lobule')
    ktb_text(48,11.8,'MxIF composite',bold=TRUE,just='center');ktb_text(137,11.8,'Cell classification',bold=TRUE,just='center')
    ktb_rect(8,17,80,52);ktb_text(48,38,'MxIF image pending',just='center',col=C$gray)
    cellmap(U[2],97,17,80,52,TRUE);ktb_rect(97,17,80,52,NA);mapkey(74,TRUE)
    ktb_title('b',2,84.8,'From cell neighborhoods to lobule enrichment')
    z<-S[distance<=43,,drop=FALSE];xx<-z$x_um-S$x_um[index];yy<-z$y_um-S$y_um[index]
    X<-9+50*(xx+42)/84;Y<-91+50*(yy+42)/84;ok<-abs(xx)<=42&abs(yy)<=42
    e<-z$ER_call==1;p<-z$p63_call==1;col<-ifelse(e&!p,C$rust,ifelse(p,C$teal,C$other));border<-ifelse(e&p,C$rust,col)
    ktb_symbol(X[ok],Y[ok],21,border[ok],col[ok],1.2,.5)
    ktb_symbol(34,116,21,C$ink,C$rust,1.9,.8)
    # A true 15-um circle on an 84-um square field.
    grid::grid.circle(x=grid::unit(34,'mm'),y=grid::unit(.graphics_state$height-116,'mm'),
      r=grid::unit(15*50/84,'mm'),gp=grid::gpar(col=C$ink,fill=NA,lwd=.8*4/3,lty=2))
    ktb_line(34,116,34+15*50/84,116);ktb_text(34+7.5*50/84,113,'15 µm',just='center')
    ktb_text(69,96,'Observed: mean p63+ neighbors around ER+ cells')
    ktb_text(69,104,'Expected: mean called epithelial neighbors');ktb_text(86,109,'× lobule p63+ fraction')
    ktb_text(69,118,'Directional index: log[(observed + ε)/(expected + ε)]')
    ktb_text(69,126,'Average the ER→p63 and p63→ER directions.')
    # Use Arial's ordinary digits at superscript position; the Unicode
    # superscript glyphs can otherwise trigger a font fallback.
    prefix<-'Radius = 15 µm; ε = 10'
    ktb_text(69,134,prefix,col=C$gray)
    tw<-function(label,size)grid::convertWidth(grid::grobWidth(grid::textGrob(label,
      gp=grid::gpar(fontfamily=config$font,fontsize=size))),'mm',TRUE)
    sx<-69+tw(prefix,8.5)
    ktb_text(sx,132.9,'−6',size=6.2,col=C$gray)
    ktb_text(sx+tw('−6',6.2),134,'; index cell excluded.',col=C$gray)
    ktb_text(11,142,'Illustrated: 2 p63+ of 9 neighboring cells.',col=C$gray)
    ktb_title('c',2,151.8,'p63 variation within donors');ktb_title('d',95,151.8,'ER variation within donors')
    heterogeneity(dget('FigS2c'),'p63',0,157);heterogeneity(dget('FigS2d'),'ER',94,157)
    ktb_title('e',2,224.8,'Donor-level ER–p63 neighborhood enrichment')
    q<-dget('FigS2e');q<-q[is.finite(q$erp63_mean)&is.finite(q$ParousFlag),];xl<-c(min(-.2,min(q$erp63_mean)-.025),max(.45,max(q$erp63_mean)+.025))
    p<-ktb_boxplot(q,'erp63_mean',jitter[jitter$panel=='FigS2e',],TRUE,.78)
    ktb_marks(p,41,229,135,15,xl,c(-.48,1.48));ktb_axes(41,229,135,15,xl,c(-.48,1.48),xticks=c(-.2,0,.2,.4),yleft=FALSE,xlabel='Mean ER–p63 log-enrichment',xlab_offset=7)
    for(g in 0:1)ktb_text(39,229+15*(1.48-(1-g))/1.96,paste0(c('Nulliparous','Parous')[g+1],' (n=',sum(q$ParousFlag==g),')'),just='right')
  }

  interaction <- function(name,marker,x,y) {
    q<-dget(name);ln<-dget(paste0(name,'_fitted_lines'));r<-dget('enrichment_composition_interactions')
    r<-r[r$exposure=='ParousFlag'&r$outcome==marker,];stopifnot(nrow(r)==1)
    xl<-range(q$within_enrichment)+c(-.05,.05);ax<-x+16;ay<-y+12;aw<-70;ah<-48
    ln$.c<-ifelse(ln$group=='Nulliparous',C$blue,C$rust);ln$.l<-ifelse(ln$group=='Nulliparous','solid','dashed')
    q$.yy<-100*q[[marker]]
    p<-ggplot2::ggplot()+ggplot2::geom_point(data=q,ggplot2::aes(x=within_enrichment,y=.yy),colour=C$gray,alpha=.3,size=.72/.75,stroke=0)+
      ggplot2::geom_line(data=ln,ggplot2::aes(x=x,y=fit_percent,colour=.c,linetype=.l),linewidth=.39)+
      ggplot2::scale_colour_identity()+ggplot2::scale_linetype_identity()
    ktb_marks(p,ax,ay,aw,ah,xl,c(0,85));ktb_axes(ax,ay,aw,ah,xl,c(0,85),xticks=c(-.2,0,.2),yticks=c(0,20,40,60,80),
      xlabel='Within-donor centered\nER–p63 log-enrichment',ylabel=paste0(if(marker=='er_fraction')'ER' else 'p63','+ epithelial cells (%)'),xlab_offset=12,ylab_offset=12)
    ktb_text(ax,y+5,paste0('Interaction P = ',ktb_p(r$p_raw),'; Holm P = ',ktb_p(r$p_holm)))
    for(g in 0:1){Y<-ay+4+g*4.8;ktb_line(ax+2,Y,ax+7,Y,c(C$blue,C$rust)[g+1],1.1,g+1);ktb_text(ax+9,Y,c('Nulliparous','Parous')[g+1])}
  }
  s3 <- function() {
    ktb_title('a',2,4.8,'Nearest-neighbor sensitivity');ktb_title('b',91,4.8,'Four-endpoint sensitivity')
    forest(dget('FigS3a'),c('ER–p63','Ki67–p63','ER–Ki67'),0,11,88,69,c(-.1,.1),'Parity difference\n(log-distance ratio)',left=22,xticks=c(-.1,0,.1))
    forest(dget('_S3b_comp'),c('p63+ fraction','ER+ fraction','Ki67+ fraction'),91,11,92,43,c(-10,10),'Percentage points',left=29,xticks=c(-10,0,10))
    forest(dget('_S3b_sp'),'ER–p63',91,57,92,27,c(-.02,.1),'Log-enrichment units',left=29,xticks=c(0,.05,.1))
    if(!.graphics_state$cropped)ktb_text(2,84,'Primary lobule-level models; Holm adjustment within each stated family.',col=C$gray)
    ktb_title('c',2,93.8,'ER and neighborhood enrichment');ktb_title('d',95,93.8,'p63 and neighborhood enrichment')
    interaction('FigS3c','er_fraction',0,98);interaction('FigS3d','p63_fraction',94,98)
  }
  s4 <- function() {
    ktb_title('a',2,4.8,'Donor-level association');ktb_title('b',96,4.8,'Simple decomposition')
    q<-dget('FigS4a');r<-dget('donor_enrichment_ER_association');q$.c<-ifelse(q$ParousFlag==0,C$blue,C$rust)
    cf<-coef(lm(yr~xr,q));ln<-data.frame(xr=seq(min(q$xr),max(q$xr),length.out=100));ln$fit<-cf[1]+cf[2]*ln$xr
    xl<-range(q$xr)+diff(range(q$xr))*c(-.05,.05);yl<-range(q$yr)+diff(range(q$yr))*c(-.05,.05)
    p<-ggplot2::ggplot(q,ggplot2::aes(x=xr,y=yr))+ggplot2::geom_point(ggplot2::aes(colour=.c),size=1.1/.75,stroke=0,alpha=.75)+
      ggplot2::geom_line(data=ln,ggplot2::aes(x=xr,y=fit),inherit.aes=FALSE,colour=C$ink,linewidth=.35)+ggplot2::scale_colour_identity()
    ktb_marks(p,15,21,73,52,xl,yl);ktb_axes(15,21,73,52,xl,yl,xticks=c(-.2,0,.2),yticks=pretty(yl,n=4),
      xlabel='ER–p63 enrichment residual\n(log-enrichment units)',ylabel='ER fraction residual\n(percentage points)',xlab_offset=11,ylab_offset=10)
    ktb_text(15,14.5,paste0('n = ',nrow(q),'; model P = ',ktb_p(r$p_raw)))
    for(g in 0:1){ktb_symbol(58,25+g*4.8,16,c(C$blue,C$rust)[g+1],size=1.5);ktb_text(61,25+g*4.8,c('Nulliparous','Parous')[g+1])}
    raw<-dget('exploratory_statistical_decompositions')
    labels<-c(total='Total parity\nassociation',direct='Direct\ncomponent',indirect_enrichment='Indirect via\nenrichment',
      indirect_p63_only='Indirect via\np63 only',indirect_enrichment_only='Indirect via\nenrichment only',serial_p63_enrichment='Serial p63 →\nenrichment',
      indirect_raw_contact_only='Indirect via\nraw contact only',serial_p63_raw_contact='Serial p63 →\nraw contact')
    for(i in 1:3){m<-c('simple','serial_p63','serial_raw_contact')[i];z<-raw[raw$model==m,]
      z$beta<-z$beta_pp;z$lower<-z$lower_pp;z$upper<-z$upper_pp
      x<-c(94,0,94)[i];y<-c(10,103,103)[i];w<-c(89,90,89)[i];h<-c(78,77,77)[i]
      if(i>1)ktb_title(c('c','d')[i-1],x+2,97.8,c('Serial decomposition','Raw-contact sensitivity')[i-1])
      forest(z,unname(labels[z$component]),x,y,w,h,c(-10,7),'Component of ER difference\n(percentage points)',pcol=NULL,left=34,xticks=c(-10,-5,0,5))}
    if(!.graphics_state$cropped)ktb_text(3,183.8,'Descriptive decompositions; intervals are 95% percentile donor-bootstrap intervals.',col=C$gray)
  }
  .graphics_state$audit<-list()
  man<-list(
    export_ktb_page(main,'Figure1',183,234,out,list(Fig1a=c(0,0,183,90),Fig1b=c(0,90,92,74),Fig1c=c(92,90,91,74),Fig1d=c(0,164,80,70),Fig1e=c(80,164,103,70))),
    export_ktb_page(s1,'FigureS1',183,244,out,list(FigS1a=c(0,0,183,77),FigS1b=c(0,79,183,83),FigS1c=c(0,164,183,80))),
    export_ktb_page(s2,'FigureS2',183,253,out,list(FigS2a=c(0,0,183,78),FigS2b=c(0,79,183,67),FigS2c=c(0,147,92,72),FigS2d=c(94,147,89,72),FigS2e=c(0,220,183,33))),
    export_ktb_page(s3,'FigureS3',183,177,out,list(FigS3a=c(0,0,90,87),FigS3b=c(91,0,92,87),FigS3c=c(0,90,92,87),FigS3d=c(94,90,89,87))),
    export_ktb_page(s4,'FigureS4',183,187,out,list(FigS4a=c(0,0,93,90),FigS4b=c(94,0,89,90),FigS4c=c(0,91,93,92),FigS4d=c(94,91,89,92))))
  manifest<-do.call(rbind,man);manifest$mode<-mode;write_csv(manifest,file.path(out,'checks/figure_manifest.csv'))
  audit<-unique(do.call(rbind,.graphics_state$audit));write_csv(audit,file.path(out,'checks/text_bounds_review.csv'))
  if(any(!audit$approx_in_bounds))warning('Some approximate text bounds require inspection. See checks/text_bounds_review.csv; this is not a visual PASS.')
  invisible(manifest)
}
