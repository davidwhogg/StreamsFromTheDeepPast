;+
; NAME:
;   switch_pm.pro
; BUGS:
;   - Not written.
;   - Header not written.
;-
pro switch_pm

pm= mrdfits('../../data/pm_monoceros.fits',1)
pm= pm[where(pm.i LT 19.5)]

out= (pm.ra gt 130.0)
in= (pm.ra lt 130.0)
outindx= where(out)
inindx= where(in)
gi= (((pm.g-pm.i) gt 0.7) and ((pm.g-pm.i) lt 0.8) and $
     (pm.i gt 17.0) and (pm.i lt 19.5))
giindx= where(gi)
notindx= where(gi EQ 0)

set_plot, 'ps'
hogg_plot_defaults
xsize= 7.5 & ysize= 7.5
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8

xrange= [0.0,1.0]
yrange= [19.5,17.0]
hogg_scatterplot, (pm[outindx].g-pm[outindx].i),pm[outindx].i, $
  xnpix=xnpix,xrange=xrange,xtitle='(g-i)', $
  ynpix=ynpix,yrange=yrange,ytitle='i', $
  grid=outgrid,/nocontours,title='out of stream'
hogg_scatterplot, (pm[inindx].g-pm[inindx].i),pm[inindx].i, $
  xnpix=xnpix,xrange=xrange,xtitle='(g-i)', $
  ynpix=ynpix,yrange=yrange,ytitle='i', $
  grid=ingrid,/nocontours,title='in stream'

zero= ((pm.i lt 19.0) AND ((pm.g-pm.i) lt 0.6))
ratio= total(in AND zero)/total(out AND zero)
splog, 'ratio',ratio
diffgrid= ingrid-outgrid*ratio
;diffgrid= diffgrid-min(diffgrid)

hogg_scatterplot, (pm[outindx].g-pm[outindx].i),pm[outindx].i, $
  xnpix=xnpix,xrange=xrange,xtitle='(g-i)', $
  ynpix=ynpix,yrange=yrange,ytitle='i', $
  grid=diffgrid,/usegrid,/nocontours,title='in - out'

xrange= [150.0,110.0]
yrange= [-5.0,30.0]
hogg_scatterplot, pm[notindx].ra,pm[notindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=notgrid,/nocontours,title='out of gi box'
hogg_scatterplot, pm[giindx].ra,pm[giindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=gigrid,/nocontours,title='in gi box'

ratio= total(gigrid)/total(notgrid)
splog, 'ratio',ratio
diffgrid= gigrid-notgrid*ratio
;diffgrid= diffgrid-min(diffgrid)

hogg_scatterplot, pm[giindx].ra,pm[giindx].dec, $
  xnpix=xnpix,xrange=xrange,xtitle='RA', $
  ynpix=ynpix,yrange=yrange,ytitle='Dec', $
  grid=diffgrid,/usegrid,/nocontours,title='in - not'

return
end
