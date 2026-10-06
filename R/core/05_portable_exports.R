# Cross-platform graphics exporter. Overrides export_ktb_page() from 03_graphics.R;
# statistical helpers are unaffected. Uses the Cairo package device (not grDevices::cairo_pdf).
ktb_pdf_device <- function()getOption('ktb.pdf_device','cairo')
ktb_portable_packages <- function()c('Cairo','pdftools','officer','zip','xml2')
ktb_portable_require <- function() {
  p<-ktb_portable_packages();miss<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss))stop('Missing packages: ',paste(miss,collapse=', '),'. Run 00_INSTALL_PACKAGES.R --yes.')
}
# Cairo 1.7-0 adds 0.5 pt internally to grid text sizes. Measure the actual
# device instead of changing the approved style sizes or assuming a version.
# The correction affects only Cairo PDF calls; SVG and ragg keep their inputs.
ktb_cairo_type_calibration <- local({
  cached <- NULL
  function(out) {
    if(is.null(cached)) {
      requested<-c(8.5,9,10,14);labels<-c('general','axis','title','letter')
      probe<-function(offset) {
        path<-tempfile(fileext='.pdf');on.exit(unlink(path),add=TRUE)
        Cairo::Cairo(file=path,type='pdf',width=3,height=2,units='in',dpi=72,bg='white')
        dev<-grDevices::dev.cur()
        on.exit(if(dev%in%grDevices::dev.list())grDevices::dev.off(dev),add=TRUE)
        grid::grid.newpage()
        for(i in seq_along(requested))grid::grid.text(labels[i],x=.5,y=1-i/5,
          gp=grid::gpar(fontfamily=config$font,fontsize=requested[i]-offset))
        grDevices::dev.off(dev)
        d<-pdftools::pdf_data(path,font_info=TRUE)[[1]]
        stopifnot(all(labels%in%d$text))
        d$font_size[match(labels,d$text)]
      }
      initial<-probe(0);delta<-initial-requested
      if(any(!is.finite(delta))||diff(range(delta))>.001||abs(delta[1])>1)
        stop('Unexpected Cairo text-size behavior; inspect before exporting.')
      correction<-mean(delta);actual<-probe(correction)
      if(any(abs(actual-requested)>.01))stop('Cairo physical text-size calibration failed.')
      cached<-data.frame(requested_pt=requested,uncorrected_pt=initial,
        correction_pt=correction,measured_pt=actual,pass=abs(actual-requested)<=.01)
    }
    write_csv(cached,file.path(out,'checks/Cairo_text_size_calibration.csv'))
    cached$correction_pt[1]
  }
})
ktb_cairo_grob <- function(g, text_correction) {
  # Cairo's vector driver uses one-point lwd units (R's standard is 0.75 pt)
  # and adds 0.5 device point to circle radii. Correct the recorded R grobs,
  # including ggplot point grobs, before drawing; never rescale the page.
  if(!is.null(g$gp$lwd))g$gp$lwd<-g$gp$lwd*.75
  if(inherits(g,'text')&&!is.null(g$gp$fontsize))g$gp$fontsize<-g$gp$fontsize-text_correction
  if(inherits(g,'circle'))g$r<-g$r-grid::unit(.5,'bigpts')
  if(inherits(g,'points')) {
    circular<-g$pch%in%c(16,19)
    if(any(circular)) {
      if(!all(circular)||is.null(g$gp$fontsize))stop('Inspect mixed or inherited point sizes before Cairo export.')
      # R uses radius = 0.375 * fontsize for these filled circular pch values.
      g$gp$fontsize<-g$gp$fontsize-4/3
      if(any(g$gp$fontsize<=0))stop('Point too small for this Cairo radius correction.')
    }
  }
  if(length(g$children))g$children<-do.call(grid::gList,lapply(g$children,ktb_cairo_grob,text_correction=text_correction))
  if(length(g$grobs))g$grobs<-lapply(g$grobs,ktb_cairo_grob,text_correction=text_correction)
  if(inherits(g,'recordedGrob')&&length(g$list))g$list<-lapply(g$list,function(x)
    if(inherits(x,'grob'))ktb_cairo_grob(x,text_correction) else x)
  g
}
# Add a literal diamond in physical mm, preserving the tested symbols otherwise.
.ktb_base_symbol <- ktb_symbol
ktb_symbol <- function(x,y,shape=16,col=ktb_style$teal,fill=col,size=2,lwd=.7) {
  if(shape!=23)return(.ktb_base_symbol(x,y,shape,col,fill,size,lwd))
  for(i in seq_along(x))grid::grid.polygon(
    x=grid::unit(x[i]+c(0,size/2,0,-size/2),'mm'),
    y=grid::unit(.graphics_state$height-y[i]+c(size/2,0,-size/2,0),'mm'),
    gp=grid::gpar(col=col,fill=fill,lwd=lwd*4/3))
}

ktb_pdf_audit <- function(path,w,h,out) {
  sz<-pdftools::pdf_pagesize(path);if(nrow(sz)!=1L)stop('Expected a one-page PDF: ',path)
  expected<-c(w,h)*72/25.4;observed<-c(sz$width[1],sz$height[1])
  exact_size<-all(abs(expected-observed)<=.01)
  fonts<-pdftools::pdf_fonts(path)
  arial<-nrow(fonts)>0&&all(grepl('Arial',fonts$name,ignore.case=TRUE))
  embedded<-if('embedded'%in%names(fonts))all(as.logical(fonts$embedded)) else FALSE
  sizes<-sort(unique(pdftools::pdf_data(path,font_info=TRUE)[[1]]$font_size))
  allowed<-c(6.2,8,8.5,9,10,14) # Approved hierarchy, map labels and S2 superscript.
  type_ok<-length(sizes)>0L&&all(vapply(sizes,function(s)min(abs(s-allowed))<=.01,logical(1)))
  rec<-data.frame(file=basename(path),device=ktb_pdf_device(),width_mm=w,height_mm=h,
    pdf_width_pt=observed[1],pdf_height_pt=observed[2],dimensions_pass=exact_size,
    Arial_only=arial,all_fonts_embedded=embedded,
    actual_type_sizes_pt=paste(sizes,collapse=';'),type_sizes_pass=type_ok,stringsAsFactors=FALSE)
  dir.create(file.path(out,'checks/pdf_devices'),recursive=TRUE,showWarnings=FALSE)
  write_csv(rec,file.path(out,'checks/pdf_devices',paste0(basename(path),'.csv')))
  write_csv(fonts[,setdiff(names(fonts),'file'),drop=FALSE],file.path(out,'checks/pdf_devices',paste0(basename(path),'.fonts.csv')))
  if(!exact_size)stop('PDF physical dimensions failed the 0.01-point check: ',basename(path),'. Repair the device; do not rescale the artwork or relax the gate.')
  if(!arial||!embedded)stop('PDF font check failed: ',basename(path),'. Inspect the font report; no silent font substitution.')
  if(!type_ok)stop('PDF physical text-size check failed: ',basename(path),'. Keep the approved point sizes.')
  invisible(rec)
}

export_ktb_page <- function(draw,name,w,h,out,crops=list(),dpi=config$dpi) {
  ktb_check_font();ktb_portable_require()
  for(d in c('figures','panels','checks'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  render<-function(path,type,crop=NULL) {
    ww<-if(is.null(crop))w else crop[3];hh<-if(is.null(crop))h else crop[4]
    cairo_pdf<-type=='pdf'&&ktb_pdf_device()=='cairo'
    correction<-if(cairo_pdf)ktb_cairo_type_calibration(out) else 0
    if(type=='pdf') {
      device<-ktb_pdf_device()
      if(device=='cairo')Cairo::Cairo(file=path,type='pdf',width=ww/25.4,height=hh/25.4,
        units='in',dpi=72,pointsize=9,bg='white',canvas='white')
      else if(device=='quartz') {
        if(!capabilities('aqua'))stop('Quartz was requested but is unavailable (macOS only); use --pdf-device cairo.')
        grDevices::quartz(type='pdf',file=path,width=ww/25.4,height=hh/25.4,family=config$font)
      } else stop('Unknown explicit PDF device: ',device)
    }
    if(type=='svg')svglite::svglite(path,width=ww/25.4,height=hh/25.4,bg='white')
    if(type=='png')ragg::agg_png(path,width=ww,height=hh,units='mm',res=dpi,background='white',scaling=1)
    dev<-grDevices::dev.cur();on.exit(if(dev%in%grDevices::dev.list())grDevices::dev.off(dev),add=TRUE)
    scene<-function() {
      grid::grid.newpage();.graphics_state$width<-w;.graphics_state$height<-h;.graphics_state$figure<-name
      .graphics_state$cropped<-!is.null(crop)
      cx<-if(is.null(crop))0 else crop[1];cy<-if(is.null(crop))0 else crop[2]
      grid::pushViewport(grid::viewport(x=grid::unit(-cx,'mm'),y=grid::unit(hh+cy,'mm'),
        width=grid::unit(w,'mm'),height=grid::unit(h,'mm'),just=c('left','top'),clip='off'))
      draw();grid::popViewport()
    }
    if(cairo_pdf) {
      # Record on the portable SVG device to keep R's standard measurements.
      # Only the native grobs are replayed into Cairo, not the SVG artwork.
      scratch<-tempfile(fileext='.svg');on.exit(unlink(scratch),add=TRUE)
      g<-grid::grid.grabExpr(scene(),width=ww/25.4,height=hh/25.4,
        device=function(width,height)svglite::svglite(scratch,width=width,height=height),
        wrap=TRUE,wrap.grobs=TRUE)
      grid::grid.newpage();grid::grid.draw(ktb_cairo_grob(g,correction))
    } else scene()
    grDevices::dev.off(dev)
    if(type=='pdf')ktb_pdf_audit(path,ww,hh,out)
  }
  for(e in c('pdf','png','svg'))render(file.path(out,'figures',paste0(name,'.',e)),e)
  for(id in names(crops))for(e in c('pdf','png','svg'))render(file.path(out,'panels',paste0(id,'.',e)),e,crops[[id]])
  manifest<-do.call(rbind,lapply(names(crops),function(id) {
    z<-crops[[id]];data.frame(figure=name,panel=id,left_mm=z[1],top_mm=z[2],width_mm=z[3],height_mm=z[4],
      figure_width_mm=w,figure_height_mm=h,stringsAsFactors=FALSE)
  }))
  if(nrow(manifest))write_csv(manifest,file.path(out,'checks',paste0(name,'_panel_layout.csv')))
  data.frame(figure=name,width_mm=w,height_mm=h,font=config$font,
    pdf=paste('Native R vector; device =',ktb_pdf_device()),svg='svglite vector',png_dpi=dpi)
}

ktb_export_powerpoint <- function(root,out) {
  # Each panel is one independent image object. This is intentionally NOT a claim
  # that individual statistical dots/labels are editable PowerPoint chart objects.
  # SVGs remain available separately for vector insertion; PNG is the reliable
  # compatibility default in the assembly deck, with no font substitution on open.
  files<-list.files(file.path(out,'checks'),pattern='_panel_layout\\.csv$',full.names=TRUE)
  if(!length(files))stop('No per-panel layout manifests available for PowerPoint.')
  dir.create(file.path(out,'powerpoint'),showWarnings=FALSE)
  layout<-do.call(rbind,lapply(files,read_csv))
  annotations<-list()
  for(name in unique(layout$figure)) {
    q<-layout[layout$figure==name,];w<-q$figure_width_mm[1];h<-q$figure_height_mm[1]
    panel_count<-nrow(q)
    template<-file.path(root,'templates',paste0('blank_',format(w,trim=TRUE),'x',format(h,trim=TRUE),'mm.pptx'))
    if(!file.exists(template))stop('Missing exact-size PPTX template: ',template)
    ppt<-officer::read_pptx(template);sz<-officer::slide_size(ppt)
    if(any(abs(c(sz$width,sz$height)*25.4-c(w,h))>1e-5))stop('PowerPoint template size mismatch.')
    if(length(ppt)>0L)for(i in rev(seq_len(length(ppt))))ppt<-officer::remove_slide(ppt,index=i)
    ls<-officer::layout_summary(ppt);i<-which(ls$layout=='Blank')[1]
    if(is.na(i))stop('No Blank PowerPoint layout.')
    ppt<-officer::add_slide(ppt,layout=ls$layout[i],master=ls$master[i])
    for(i in seq_len(nrow(q))) {
      v<-q[i,];image<-file.path(out,'panels',paste0(v$panel,'.png'));if(!file.exists(image))stop('Missing panel PNG: ',image)
      ppt<-officer::ph_with(ppt,officer::external_img(image,width=v$width_mm,height=v$height_mm,unit='mm',
        alt=paste(v$panel,'; exported from this run; see results and legend')),
        location=officer::ph_location(left=v$left_mm/25.4,top=v$top_mm/25.4,width=v$width_mm/25.4,
          height=v$height_mm/25.4,newlabel=v$panel),use_loc_size=TRUE)
    }
    # The Figure 1 standalone panel crops intentionally omit shared notes.
    # Restore those notes separately in the complete PowerPoint assembly,
    # taking the actual text/position from the full-page R drawing audit.
    if(name%in%c('FigureS3','FigureS4')) {
      a<-read_csv(file.path(out,'checks/text_bounds_review.csv'))
      yy<-if(name=='FigureS3')84 else 183.8
      note<-unique(a[a$figure==name&abs(a$y_mm-yy)<1e-8,c('text','x_mm','y_mm')])
      if(nrow(note)!=1L)stop('Missing/nonunique shared Figure1 note: ',name)
      id<-paste0(name,'_shared_note');top<-yy-2.5
      dir.create(file.path(out,'annotations'),showWarnings=FALSE)
      png<-file.path(out,'annotations',paste0(id,'.png'))
      ragg::agg_png(png,width=w,height=5,units='mm',res=config$dpi,background='white')
      grid::grid.newpage();grid::grid.text(note$text,x=grid::unit(note$x_mm,'mm'),
        y=grid::unit(2.5,'mm'),just='left',vjust=.5,
        gp=grid::gpar(fontfamily=config$font,fontsize=8.5,col=ktb_style$gray,lineheight=1.05))
      grDevices::dev.off()
      ppt<-officer::ph_with(ppt,officer::external_img(png,width=w,height=5,unit='mm',alt=note$text),
        location=officer::ph_location(left=0,top=top/25.4,width=w/25.4,height=5/25.4,newlabel=id),use_loc_size=TRUE)
      v<-data.frame(figure=name,panel=id,left_mm=0,top_mm=top,width_mm=w,height_mm=5,
        figure_width_mm=w,figure_height_mm=h,stringsAsFactors=FALSE)
      q<-rbind(q,v);annotations[[length(annotations)+1L]]<-v
    }
    dest<-file.path(out,'powerpoint',paste0(name,'_assembly.pptx'));print(ppt,target=dest)
    # Reopen and inspect shapes; no flattened page may replace independent panels.
    contents<-utils::unzip(dest,list=TRUE)$Name
    stopifnot('ppt/slides/slide1.xml'%in%contents)
    td<-tempfile('ppt_audit');dir.create(td);utils::unzip(dest,files=c('ppt/slides/slide1.xml','ppt/presentation.xml'),exdir=td)
    doc<-xml2::read_xml(file.path(td,'ppt/slides/slide1.xml'))
    pics<-xml2::xml_find_all(doc,'//*[local-name()="pic"]')
    if(length(pics)!=nrow(q))stop('PPTX does not contain one separate image per panel: ',name)
    for(i in seq_len(nrow(q))) {
      xf<-xml2::xml_find_first(pics[[i]],'.//*[local-name()="xfrm"]')
      off<-xml2::xml_find_first(xf,'./*[local-name()="off"]');ex<-xml2::xml_find_first(xf,'./*[local-name()="ext"]')
      got<-as.numeric(c(xml2::xml_attr(off,'x'),xml2::xml_attr(off,'y'),xml2::xml_attr(ex,'cx'),xml2::xml_attr(ex,'cy')))/36000
      expected<-as.numeric(unlist(q[i,c('left_mm','top_mm','width_mm','height_mm')],use.names=FALSE))
      if(any(abs(got-expected)>1e-4))stop('PowerPoint panel placement/size mismatch: ',q$panel[i])
      lab<-xml2::xml_attr(xml2::xml_find_first(pics[[i]],'.//*[local-name()="cNvPr"]'),'name')
      if(!identical(lab,q$panel[i]))stop('PowerPoint Selection Pane label mismatch: ',q$panel[i])
    }
    write_csv(data.frame(figure=name,expected_panels=panel_count,picture_objects=length(pics),
      shared_note_objects=nrow(q)-panel_count,
      independent_panels=TRUE,artwork='600-dpi panel PNGs; individual SVGs supplied separately'),
      file.path(out,'checks',paste0(name,'_PowerPoint_audit.csv')))
    unlink(td,recursive=TRUE)
  }
  write_csv(layout,file.path(out,'powerpoint/PANEL_DIMENSIONS_MM.csv'))
  if(length(annotations))write_csv(do.call(rbind,annotations),file.path(out,'powerpoint/ANNOTATION_DIMENSIONS_MM.csv'))
  writeLines(c('Each assembly file contains one figure at publication dimensions, with every panel as an independent movable picture.',
    'Insert individual panels from panels/*.svg for vector artwork in SVG-capable PowerPoint, or use panels/*.png for broad compatibility.',
    'The PNGs are 600 dpi at the recorded physical size. Keep Lock aspect ratio selected and use the original mm sizes for manuscript assembly.',
    'Do not shrink the whole figure and then expect the lettering to retain its publication size.',
    'The deck does not turn statistical data into editable Office charts. R and the source tables remain the place to change data or statistical labels.'),
    file.path(out,'powerpoint/README.txt'))
  invisible(layout)
}
