# Shared fixed-physical-size graphics layer. Geometry is in millimetres; type in pt.
# ggplot2 draws statistical marks; grid assembles them without resizing the text.
# All statistical panels are drawn natively in R.

ktb_style <- list(font='Arial',letter=14,title=10,axis=9,text=8.5,
  blue='#377C99',rust='#A45137',teal='#237B7C',gray='#727B83',ink='#202428',
  other='#B8BEC2',line='#CBD0D3',pale='#F3F4F4',grid='#E5E8EA')
.graphics_state <- new.env(parent=emptyenv())
.graphics_state$width <- 183; .graphics_state$height <- 234
.graphics_state$audit <- list(); .graphics_state$figure <- ''

ktb_check_font <- function(font=config$font) {
  fonts <- systemfonts::system_fonts()
  hits <- fonts[tolower(fonts$family)==tolower(font),,drop=FALSE]
  if(!nrow(hits))stop('Arial is required. Install the Arial font locally; it is not bundled with this code.')
  invisible(hits)
}
ktb_p <- function(x) ifelse(!is.finite(x),'NA',ifelse(x<.0001,formatC(x,format='e',digits=1),ifelse(x<.01,sprintf('%.4f',x),sprintf('%.3f',x))))
ktb_num <- function(x,digits=1,signed=FALSE) sub('-', '\u2212',sprintf(paste0('%',if(signed)'+' else '', '.',digits,'f'),x),fixed=TRUE)
ktb_alpha <- function(col,a) grDevices::adjustcolor(col,alpha.f=a)
ktb_theme <- function() ggplot2::theme_void(base_family=config$font,base_size=9)+
  ggplot2::theme(plot.margin=ggplot2::margin(0,0,0,0),panel.background=ggplot2::element_rect(fill='white',colour=NA),
    plot.background=ggplot2::element_rect(fill='white',colour=NA),legend.position='none')

ktb_text <- function(x,y,label,size=8.5,bold=FALSE,just='left',col=ktb_style$ink,rot=0,vjust=.5) {
  h <- .graphics_state$height
  g <- grid::textGrob(label,x=grid::unit(x,'mm'),y=grid::unit(h-y,'mm'),
    just=just,vjust=vjust,rot=rot,gp=grid::gpar(fontfamily=config$font,fontsize=size,
      fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.05))
  grid::grid.draw(g)
  # Approximate text boxes in physical units, including rotated axis titles.
  w<-grid::convertWidth(grid::grobWidth(g),'mm',TRUE);hh<-grid::convertHeight(grid::grobHeight(g),'mm',TRUE)
  # grobWidth/grobHeight already include the text rotation.
  j<-if(is.numeric(just))just else switch(just,left=0,center=.5,centre=.5,right=1)
  .graphics_state$audit[[length(.graphics_state$audit)+1L]]<-data.frame(
    figure=.graphics_state$figure,text=label,x_mm=x,y_mm=y,width_mm=w,height_mm=hh,
    approx_in_bounds=(x-j*w>=-.5 && x+(1-j)*w<=.graphics_state$width+.5 && y-hh/2>=-.5 && y+hh/2<=h+.5))
  invisible(g)
}
ktb_title <- function(letter,x,y,title) {ktb_text(x,y,letter,14,TRUE);ktb_text(x+7,y,title,10,TRUE)}
ktb_rect <- function(x,y,w,h,fill=ktb_style$pale,col=ktb_style$line,lwd=.6) {
  grid::grid.rect(x=grid::unit(x,'mm'),y=grid::unit(.graphics_state$height-y,'mm'),
    width=grid::unit(w,'mm'),height=grid::unit(h,'mm'),just=c('left','top'),
    gp=grid::gpar(fill=fill,col=col,lwd=lwd*4/3))
}
ktb_line <- function(x0,y0,x1,y1,col=ktb_style$ink,lwd=.65,lty=1) {
  grid::grid.segments(x0=grid::unit(x0,'mm'),y0=grid::unit(.graphics_state$height-y0,'mm'),
    x1=grid::unit(x1,'mm'),y1=grid::unit(.graphics_state$height-y1,'mm'),
    gp=grid::gpar(col=col,lwd=lwd*4/3,lty=lty,lineend='butt'))
}
ktb_symbol <- function(x,y,shape=16,col=ktb_style$teal,fill=col,size=2,lwd=.7) {
  # R pch glyphs do not have the requested literal diameter. Draw the approved
  # geometry directly so size is a diameter/side length in mm and lwd is in pt.
  yy<-.graphics_state$height-y
  if(shape%in%c(16,21,1))grid::grid.circle(x=grid::unit(x,'mm'),y=grid::unit(yy,'mm'),
    r=grid::unit(size/2,'mm'),gp=grid::gpar(col=if(shape==16)NA else col,
      fill=if(shape==1)NA else fill,lwd=lwd*4/3))
  else if(shape==22)grid::grid.rect(x=grid::unit(x,'mm'),y=grid::unit(yy,'mm'),
    width=grid::unit(size,'mm'),height=grid::unit(size,'mm'),gp=grid::gpar(col=col,fill=fill,lwd=lwd*4/3))
  else if(shape==4)for(s in c(-1,1))grid::grid.segments(
    x0=grid::unit(x-size/2,'mm'),x1=grid::unit(x+size/2,'mm'),
    y0=grid::unit(yy-s*size/2,'mm'),y1=grid::unit(yy+s*size/2,'mm'),
    gp=grid::gpar(col=col,lwd=lwd*4/3))
  else stop('Unsupported physical symbol: ',shape)
}

# Draw only ggplot's plotting region into the exact supplied rectangle; axes and
# labels are placed separately, in mm. This avoids panel-specific font scaling.
ktb_marks <- function(p,x,y,w,h,xlim,ylim,ygrid=numeric()) {
  p<-p+ggplot2::coord_cartesian(xlim=xlim,ylim=ylim,expand=FALSE,clip='on')+ktb_theme()
  if(length(ygrid))p<-p+ggplot2::scale_y_continuous(breaks=ygrid)+
    ggplot2::theme(panel.grid.major.y=ggplot2::element_line(colour=ktb_style$grid,linewidth=.14))
  g<-ggplot2::ggplotGrob(p);ix<-which(g$layout$name=='panel')
  if(length(ix)!=1L)stop('One plotting region required; assemble facets explicitly.')
  grid::pushViewport(grid::viewport(x=grid::unit(x,'mm'),y=grid::unit(.graphics_state$height-y,'mm'),
    width=grid::unit(w,'mm'),height=grid::unit(h,'mm'),just=c('left','top'),clip='on'))
  grid::grid.draw(g$grobs[[ix]]);grid::popViewport()
}
ktb_axes <- function(x,y,w,h,xlim,ylim,xticks=numeric(),yticks=numeric(),xlabel=NULL,ylabel=NULL,
                     xticklabels=NULL,yticklabels=NULL,yleft=TRUE,bottom=TRUE,xlab_offset=9,ylab_offset=12,
                     grid_y=FALSE) {
  ix<-which(xticks>=xlim[1]&xticks<=xlim[2]);iy<-which(yticks>=ylim[1]&yticks<=ylim[2])
  xticks<-xticks[ix];yticks<-yticks[iy]
  if(!is.null(xticklabels))xticklabels<-xticklabels[ix]
  if(!is.null(yticklabels))yticklabels<-yticklabels[iy]
  X<-function(v)x+w*(v-xlim[1])/diff(xlim);Y<-function(v)y+h*(ylim[2]-v)/diff(ylim)
  if(grid_y)for(v in yticks)ktb_line(x,Y(v),x+w,Y(v),ktb_style$grid,.4)
  if(bottom)ktb_line(x,y+h,x+w,y+h)
  if(yleft)ktb_line(x,y,x,y+h)
  if(length(xticks))for(i in seq_along(xticks)){v<-xticks[i];ktb_line(X(v),y+h,X(v),y+h+.75,lwd=.6);ktb_text(X(v),y+h+2.8,
    if(is.null(xticklabels))format(v,trim=TRUE) else xticklabels[i],just='center')}
  if(length(yticks))for(i in seq_along(yticks)){v<-yticks[i];ktb_line(x-.75,Y(v),x,Y(v),lwd=.6);ktb_text(x-1.8,Y(v),
    if(is.null(yticklabels))format(v,trim=TRUE) else yticklabels[i],just='right')}
  if(!is.null(xlabel))ktb_text(x+w/2,y+h+xlab_offset,xlabel,9,just='center')
  if(!is.null(ylabel))ktb_text(x-ylab_offset,y+h/2,ylabel,9,just='center',rot=90)
  invisible(list(X=X,Y=Y))
}
ktb_ci <- function(b,lo,hi,yy,x,y,w,h,xlim,ylim,col=ktb_style$teal,shape=16,fill=col,size=1.85) {
  if(any(!is.finite(c(b,lo,hi)))||lo>b||hi<b)stop('Invalid confidence interval')
  if(lo<xlim[1]-1e-10||hi>xlim[2]+1e-10)stop('Confidence interval outside fixed axis; flag, do not clip silently.')
  X<-function(v)x+w*(v-xlim[1])/diff(xlim);Y<-y+h*(ylim[2]-yy)/diff(ylim)
  ktb_line(X(lo),Y,X(hi),Y,col,.9);ktb_line(X(lo),Y-.75,X(lo),Y+.75,col,.6);ktb_line(X(hi),Y-.75,X(hi),Y+.75,col,.6)
  ktb_symbol(X(b),Y,shape,col,fill,size,.75)
}
ktb_zero <- function(x,y,w,h,xlim) if(xlim[1]<=0&&xlim[2]>=0)ktb_line(x+w*(-xlim[1])/diff(xlim),y,x+w*(-xlim[1])/diff(xlim),y+h,'#99A1A6',.6,2)

ktb_box_stats <- function(x) {
  x<-x[is.finite(x)];qs<-quantile(x,c(.25,.5,.75),type=7,names=FALSE);iq<-qs[3]-qs[1]
  c(lo=min(x[x>=qs[1]-1.5*iq]),q1=qs[1],med=qs[2],q3=qs[3],hi=max(x[x<=qs[3]+1.5*iq]))
}
ktb_boxplot <- function(d,var,jitter, horizontal=FALSE, point_mm=1.3) {
  # ggplot's zero-stroke circular point diameter is 0.75 times its size argument.
  point_mm<-point_mm/.75
  d<-d[is.finite(d[[var]])&is.finite(d$ParousFlag),,drop=FALSE]
  d$.v<-d[[var]];d$.g<-d$ParousFlag
  d$.off<-jitter$offset[match(d$Barcode,jitter$Barcode)];if(anyNA(d$.off))stop('Missing display jitter key')
  d$.c<-ifelse(d$.g==0,ktb_style$blue,ktb_style$rust)
  stats<-do.call(rbind,lapply(0:1,function(g){v<-ktb_box_stats(d$.v[d$.g==g]);data.frame(g=g,t(v),c=if(g==0)ktb_style$blue else ktb_style$rust)}))
  if(!horizontal){
    p<-ggplot2::ggplot()+ggplot2::geom_segment(data=stats,ggplot2::aes(x=g,xend=g,y=lo,yend=hi,colour=c),linewidth=.25)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(x=g-.12,xend=g+.12,y=lo,yend=lo,colour=c),linewidth=.25)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(x=g-.12,xend=g+.12,y=hi,yend=hi,colour=c),linewidth=.25)+
      ggplot2::geom_rect(data=stats,ggplot2::aes(xmin=g-.24,xmax=g+.24,ymin=q1,ymax=q3,colour=c,fill=c),alpha=.12,linewidth=.25)+
      ggplot2::geom_point(data=d,ggplot2::aes(x=.g+.off,y=.v,colour=.c),size=point_mm,alpha=.68,stroke=0)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(x=g-.24,xend=g+.24,y=med,yend=med,colour=c),linewidth=.55)
  }else{
    # Nulliparous above parous, matching the approved horizontal S2e display.
    stats$g<-1-stats$g;d$.g<-1-d$.g
    p<-ggplot2::ggplot()+ggplot2::geom_segment(data=stats,ggplot2::aes(y=g,yend=g,x=lo,xend=hi,colour=c),linewidth=.22)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(y=g-.14,yend=g+.14,x=lo,xend=lo,colour=c),linewidth=.22)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(y=g-.14,yend=g+.14,x=hi,xend=hi,colour=c),linewidth=.22)+
      ggplot2::geom_rect(data=stats,ggplot2::aes(ymin=g-.275,ymax=g+.275,xmin=q1,xmax=q3,colour=c,fill=c),alpha=.12,linewidth=.22)+
      ggplot2::geom_point(data=d,ggplot2::aes(y=.g+.off,x=.v,colour=.c),size=point_mm,alpha=.65,stroke=0)+
      ggplot2::geom_segment(data=stats,ggplot2::aes(y=g-.275,yend=g+.275,x=med,xend=med,colour=c),linewidth=.45)
  }
  p+ggplot2::scale_colour_identity()+ggplot2::scale_fill_identity()
}

export_ktb_page <- function(draw,name,w,h,out,crops=list(),dpi=config$dpi) {
  ktb_check_font();dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(out,'panels'),recursive=TRUE,showWarnings=FALSE)
  ext<-c('pdf','png','svg')
  render <- function(path,type,crop=NULL) {
    ww<-if(is.null(crop))w else crop[3];hh<-if(is.null(crop))h else crop[4]
    if(type=='pdf') {
      if(!capabilities('aqua'))stop('This exporter uses the macOS Quartz PDF device. Source R/core/05_portable_exports.R for the cross-platform (Cairo) exporter.')
      grDevices::quartz(type='pdf',file=path,width=ww/25.4,height=hh/25.4,family=config$font)
    }
    if(type=='png')ragg::agg_png(path,width=ww,height=hh,units='mm',res=dpi,background='white',scaling=1)
    if(type=='svg')svglite::svglite(path,width=ww/25.4,height=hh/25.4,bg='white')
    on.exit(grDevices::dev.off(),add=TRUE)
    grid::grid.newpage();.graphics_state$width<-w;.graphics_state$height<-h;.graphics_state$figure<-name
    .graphics_state$cropped<-!is.null(crop)
    cx<-if(is.null(crop))0 else crop[1];cy<-if(is.null(crop))0 else crop[2]
    grid::pushViewport(grid::viewport(x=grid::unit(-cx,'mm'),y=grid::unit(hh+cy,'mm'),
      width=grid::unit(w,'mm'),height=grid::unit(h,'mm'),just=c('left','top'),clip='off'))
    draw();grid::popViewport()
  }
  for(e in ext)render(file.path(out,'figures',paste0(name,'.',e)),e)
  for(id in names(crops))for(e in ext)render(file.path(out,'panels',paste0(id,'.',e)),e,crops[[id]])
  data.frame(figure=name,width_mm=w,height_mm=h,font=config$font,pdf='Native R Quartz vector (exact dimensions)',svg='svglite vector',png_dpi=dpi)
}
