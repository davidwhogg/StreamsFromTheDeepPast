;+
; NAME:
;   bhb_pm.pro
; BUGS:
;   - Not written.
;   - Header not written.
;-
pro bhb_pm

pm= mrdfits('../../data/big*.fits',1)
pm= pm[where(pm.i lt 18.0)]

ugr= (((pm.u-pm.g) gt 0.8) and ((pm.u-pm.g) lt 1.35) and $
      ((pm.g-pm.r) gt (-0.4)) and ((pm.g-pm.r) lt 0.0))
ugrindx= where(ugr)
notindx= where(ugr EQ 0)

set_plot, 'ps'
hogg_plot_defaults
xsize= 7.5 & ysize= 7.5
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8

xrange= [150.0,140.0]
yrange= [5.0,30.0]
hogg_scatterplot, pm[notindx].ra,pm[notindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=notgrid,/nocontours,title='not in ugr box'
hogg_scatterplot, pm[ugrindx].ra,pm[ugrindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=ugrgrid,/nocontours,title='in ugr box'

ratio= total(ugrgrid)/total(notgrid)
splog, 'ratio',ratio
diffgrid= ugrgrid-notgrid*ratio
;diffgrid= diffgrid-min(diffgrid)

hogg_scatterplot, pm[ugrindx].ra,pm[ugrindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=diffgrid,/usegrid,/nocontours,title='in - not'

return
end
