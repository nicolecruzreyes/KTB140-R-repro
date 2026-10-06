# Native R implementation of the approved 183-mm Figure 2 and S5 layouts.
# The two publication pages and seven standalone panels share one drawing function.
render_figure2 <- function(root,out,bundle,mode) {
  C<-ktb_style;S<-c('immune','checkpoint','myeloid','ecm');K<-names(state_labels)
  labels<-c('Immune infiltration','Checkpoint /\nexhaustion','Myeloid','ECM / fibroblast')
  state_names<-c('ER–p63 epithelial\nbalance','ER–p63\nneighborhood\nenrichment',
    'Myeloid / antigen-\npresentation','Stromal /\nfibrovascular')
  get<-function(n){x<-bundle[[n]];if(is.null(x))stop('Missing figure source: ',n);x}
  jitter<-read_csv(file.path(root,'config/figure2_visual_jitter.csv'))
  qp<-function(x)ifelse(x<.001,formatC(x,format='e',digits=1),sprintf('%.3f',x))
  one<-function(d,k,v){i<-which(d[[k]]==v);if(length(i)!=1L)stop('Missing/nonunique plot row: ',k,'=',v);d[i,,drop=FALSE]}
  cardlines<-function(x,y,labs,size=8.5,bold=FALSE,leading=4.5,col=C$ink,just='left') {
    for(i in seq_along(labs))ktb_text(x,y+(i-1)*leading,labs[i],size,bold,just,col)
  }
  dsp_forest<-function(name,x,y,first,second,labs) {
    q<-get(name);ax<-x+43;ay<-y+13;aw<-96;ah<-34;xl<-c(-.15,.70);yl<-c(-.55,3.55)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(-.1,0,.25,.5),
      xticklabels=c('−0.1','0','0.25','0.5'),yleft=FALSE,
      xlabel='Parous − nulliparous difference\n(mean-z signature units)',xlab_offset=11)
    ktb_zero(ax,ay,aw,ah,xl)
    secondshape<-if(second=='isotype')23 else 22
    for(g in seq_along(c(first,second))) {
      model<-c(first,second)[g];r<-q[q$model==model,,drop=FALSE];r<-r[match(S,r$signature),,drop=FALSE]
      stopifnot(nrow(r)==4L,!anyNA(r$signature));off<-if(g==1).14 else -.14
      col<-if(g==1)C$teal else C$gray;shape<-if(g==1)16 else secondshape
      for(i in 1:4) {
        yy<-4-i;Y<-ay+ah*(yl[2]-yy)/diff(yl)
        ktb_ci(r$beta[i],r$lower[i],r$upper[i],yy+off,ax,ay,aw,ah,xl,yl,
          col=col,shape=shape,fill=if(g==1)col else 'white',size=1.55)
        ktb_text(x+c(150,173)[g],Y,qp(r$p_holm[i]),just='center')
        if(g==1)ktb_text(ax-5,Y,labels[i],just='right')
      }
    }
    # Model legend occupies the left of the same header row as the P-value columns.
    for(i in 1:2) {
      xx<-x+c(47,80)[i];col<-if(i==1)C$teal else C$gray
      ktb_line(xx-2,y+5,xx+2,y+5,col,.9)
      ktb_symbol(xx,y+5,if(i==1)16 else secondshape,col,if(i==1)col else 'white',1.55)
      ktb_text(xx+4,y+5,labs[i])
    }
    ktb_text(x+161.5,y+4,'Holm P',bold=TRUE,just='center')
    for(i in 1:2)ktb_text(x+c(150,173)[i],y+9,if(first=='ROI')c('ROI','Donor')[i] else c('Base','+ Isotypes')[i],just='center')
  }
  heatmap<-function(x,y) {
    d<-get('Fig2c');M<-as.matrix(d[match(K,d$state),K,drop=FALSE]);storage.mode(M)<-'double'
    ax<-x+28;ay<-y+7;sz<-11
    ramp<-grDevices::colorRamp(c(C$blue,'#FFFFFF',C$rust))
    cmap<-function(v){rgb<-ramp((v+1)/2)/255;grDevices::rgb(rgb[,1],rgb[,2],rgb[,3])}
    for(i in 1:4)for(j in 1:4) {
      v<-M[i,j];ktb_rect(ax+(j-1)*sz,ay+(i-1)*sz,sz,sz,cmap(v),'white',.6)
      ktb_text(ax+(j-.5)*sz,ay+(i-.5)*sz,sprintf('%.2f',v),just='center',col=if(abs(v)>.72)'white' else C$ink)
    }
    yl<-c('Epithelial\nbalance','ER–p63\nenrichment','Myeloid /\npresentation','Stromal /\nfibrovascular')
    for(i in 1:4)ktb_text(ax-2,ay+(i-.5)*sz,yl[i],just='right')
    for(i in 1:4)ktb_text(ax+(i-.5)*sz,ay+47,c('Balance','ER–p63','Myeloid','Stromal')[i],just='right',rot=45)
    for(i in 0:99)ktb_rect(x+78,ay+i*44/100,2,44/100,cmap(1-2*(i+.5)/100),NA,0)
    for(v in c(-1,0,1))ktb_text(x+81,ay+22*(1-v),as.character(v))
    ktb_text(x+50,y+2,'114 donors with all four states',just='center',col=C$gray)
    ktb_text(x+48,y+70,'Pearson correlation',just='center')
  }
  state_forest<-function(x,y) {
    q<-get('Fig2d');q<-q[match(K,q$state),];ax<-x+34;ay<-y+7;aw<-41;ah<-47;xl<-c(-1.2,1.2);yl<-c(-.5,3.5)
    ktb_axes(ax,ay,aw,ah,xl,yl,xticks=c(-1,0,1),yleft=FALSE,
      xlabel='Parous − nulliparous\n(state SD units)',xlab_offset=11);ktb_zero(ax,ay,aw,ah,xl)
    ktb_text(x+83,y+2,'Holm P',bold=TRUE,just='center')
    for(i in 1:4) {
      yy<-4-i;Y<-ay+ah*(yl[2]-yy)/diff(yl);sig<-q$p_holm[i]<.05;col<-if(sig)C$teal else C$gray
      ktb_ci(q$beta[i],q$lower[i],q$upper[i],yy,ax,ay,aw,ah,xl,yl,col,
        shape=if(sig)16 else 21,fill=if(sig)col else 'white',size=1.59)
      ktb_text(ax-4,Y,state_names[i],just='right');ktb_text(x+83,Y,qp(q$p_holm[i]),just='center')
    }
  }
  main<-function() {
    cover<-get('Fig2a_coverage');ns<-cover$donors[match(K,cover$state)]
    if(!identical(as.integer(ns),c(127L,120L,124L,124L)))stop('Unexpected state coverage; do not hardcode figure sample sizes.')
    ktb_title('a',2,6,'Four donor-level tissue states')
    ktb_text(181,6,'Assays integrated within each donor',just='right',col=C$gray)
    for(z in list(c(2,12,87,29),c(94,12,87,29),c(2,45,87,34),c(94,45,87,34))) {
      ktb_rect(z[1],z[2],z[3],z[4]);ktb_rect(z[1],z[2],1.1,z[4],C$teal,NA,0)
    }
    ktb_text(6,18,'ER–p63 epithelial balance',9,TRUE)
    cardlines(6,24,c('MxIF ER fraction − MxIF p63 fraction','+ IO360 ER-epithelial signature'))
    ktb_text(6,36,paste0('Three equally weighted components; n = ',ns[1]),col=C$gray)
    ktb_text(98,18,'ER–p63 neighborhood enrichment',9,TRUE)
    cardlines(98,24,c('MxIF mean lobule log-enrichment','Within-15-µm ER–p63 neighborhoods'))
    ktb_text(98,36,paste0('Single component; n = ',ns[2]),col=C$gray)
    cardlines(6,51,c('Myeloid/antigen-presentation','microenvironment'),9,TRUE,4)
    cardlines(6,61,c('½ DSP immune/myeloid block','+ ½ IO360 myeloid/antigen-presentation'))
    ktb_text(6,74,paste0('DSP block = mean(immune, myeloid); n = ',ns[3]),col=C$gray)
    ktb_text(98,51,'Stromal/fibrovascular state',9,TRUE)
    cardlines(98,61,c('½ DSP ECM/fibroblast signature','+ ½ IO360 stromal/vascular signature'))
    ktb_text(98,74,paste0('Equal assay-layer weights; n = ',ns[4]),col=C$gray)
    ktb_text(91.5,85,'Standardized components → donor-level combination → each state standardized (mean 0, SD 1)',just='center')
    ktb_title('b',2,94,'Regional protein signatures');ktb_text(181,94,'724 ROIs from 124 donors',just='right',col=C$gray)
    dsp_forest('Fig2b',0,98,'ROI','donor',c('ROI-level','Donor-level'))
    ktb_title('c',2,167,'State correlations');ktb_title('d',91,167,'Parity associations')
    heatmap(0,170);state_forest(89,170)
  }
  supp<-function() {
    ktb_title('a',2,6,'Regional sampling');ktb_title('b',96,6,'Epithelial content')
    ktb_rect(5,15,81,47)
    cardlines(45.5,32,c('GeoMx morphology image','600-µm whole ROI','Image pending'),9,FALSE,5,C$gray,'center')
    cardlines(5,70,c('Regions preselected on adjacent H&E','Full ROI collected; no compartment segmentation'),8.5,FALSE,4.8)
    d<-get('FigS5b');r<-get('FigS5b_annotation');stopifnot(nrow(r)==1L)
    d<-d[is.finite(d$mean_epi_score)&is.finite(d$ParousFlag),]
    if(any(d$mean_epi_score< -2|d$mean_epi_score>2))stop('Epithelial-content point outside approved axis.')
    n<-table(factor(d$ParousFlag,levels=0:1))
    ktb_text(146.5,16,paste0('Mean difference = ',ktb_num(r$beta,2,TRUE)),just='center')
    ktb_text(146.5,21,paste0('Welch P = ',sprintf('%.3f',r$p_raw)),just='center')
    p<-ktb_boxplot(d,'mean_epi_score',jitter,point_mm=.97)
    ktb_marks(p,114,28,65,43,c(-.5,1.5),c(-2,2))
    ktb_axes(114,28,65,43,c(-.5,1.5),c(-2,2),xticks=0:1,yticks=-2:2,
      xticklabels=c('',''),ylabel='Mean epithelial-content score',ylab_offset=12)
    for(i in 1:2)ktb_text(114+65*(i-.5)/2,75.5,
      paste0(c('Nulliparous','Parous')[i],'\nn = ',as.integer(n)[i]),just='center')
    ktb_title('c',2,99,'Sensitivity to additional isotype adjustment');ktb_text(181,99,'124 donors',just='right',col=C$gray)
    dsp_forest('FigS5c',0,105,'donor','isotype',c('Base model','+ Isotypes'))
    ktb_text(91.5,172,'Both models use nuclei-weighted donor signature summaries.',just='center',col=C$gray)
  }
  .graphics_state$audit<-list()
  dims<-rbind(export_ktb_page(main,'Figure2',183,246,out,crops=list(
    Fig2a=c(0,0,183,88),Fig2b=c(0,89,183,73),Fig2c=c(0,162,88,84),Fig2d=c(89,162,94,84))),
    export_ktb_page(supp,'FigureS5',183,174,out,crops=list(
      FigS5a=c(0,0,91,88),FigS5b=c(94,0,89,88),FigS5c=c(0,93,183,81))))
  write_csv(dims,file.path(out,'checks/figure_dimensions.csv'))
  audit<-do.call(rbind,.graphics_state$audit);write_csv(audit,file.path(out,'checks/text_bounds_review.csv'))
  if(any(!audit$approx_in_bounds))warning('Some text bounds require inspection; see text_bounds_review.csv. Do not shrink fonts to hide a layout error.')
  invisible(dims)
}
