;+
; NAME:
;   reduced_pm.pro
; BUGS:
;   - Not written.
;   - Header not written.
;-
pro reduced_pm

pm= mrdfits('../../data/big*.fits',1)
pm= pm[where(pm.pmdecerr LT 5.0)]

totalpm= sqrt((pm.pmra*cos(pm.dec*!DPI/1.8D2))^2+(pm.pmdec)^2)
dm= 5.0*alog10(totalpm)
pm.u= pm.u+dm
pm.g= pm.g+dm
pm.r= pm.r+dm
pm.i= pm.i+dm
pm.z= pm.z+dm

set_plot, 'ps'
hogg_plot_defaults
xsize= 7.5 & ysize= 7.5
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8

xrange= [0.0,1.0]
yrange= [26.0,15.0]
hogg_scatterplot, pm.g-pm.i,pm.i, $
  xnpix=xnpix,xrange=xrange,xtitle='g - i', $
  ynpix=ynpix,yrange=yrange,ytitle='i + 5 log!d10!n !7l!3', $
  grid=grid,/nocontours,title='reduced proper motions'

return
end
