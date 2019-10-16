;+
; NAME:
;   ngp_pm.pro
; BUGS:
;   - Do I really need to multiply by cos Dec?
;   - Header not written.
;-
pro ngp_pm

prefix= 'ngp'
seed= -1L
pm= mrdfits('../../data/'+prefix+'*.fit*',1)
pm= pm[where(pm.pmDecErr lt 10.0)]
nstar= n_elements(pm)
pm.pmRA= pm.pmRA+randomu(seed,nstar)
pm.pmDec= pm.pmDec+randomu(seed,nstar)

set_plot, 'ps'
hogg_plot_defaults
xsize= 7.5 & ysize= 7.5
psfilename= prefix+'.ps'
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
!P.MULTI= [0,2,2]
!X.MARGIN= !X.OMARGIN
!X.OMARGIN= [0,0]
!Y.MARGIN= !Y.OMARGIN
!Y.OMARGIN= [0,0]

xtitle= 'RA'
ytitle= 'Dec'
hogg_scatterplot, pm.ra,pm.dec, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix

xtitle= '[g - i]'
xrange= [0.0,1.0]
ytitle= 'i'
yrange= [19.5,15.5]
hogg_scatterplot, (pm.g-pm.i),pm.i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix

xtitle= '!7l!da!n!3'
xrange= [-1,1]*25.0
hogg_scatterplot, pm.pmra,pm.i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix
xtitle= '!7l!dd!n!3'
hogg_scatterplot, pm.pmdec,pm.i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix

thickdisk= where(((pm.g-pm.i) gt 0.45) and ((pm.g-pm.i) lt 0.60))
xtitle= '!7l!da!n!3'
hogg_scatterplot, pm[thickdisk].pmra,pm[thickdisk].i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' thick disk'
xtitle= '!7l!dd!n!3'
hogg_scatterplot, pm[thickdisk].pmdec,pm[thickdisk].i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' thick disk'

halo= where(((pm.g-pm.i) gt 0.30) and ((pm.g-pm.i) lt 0.45))
xtitle= '!7l!da!n!3'
hogg_scatterplot, pm[halo].pmra,pm[halo].i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' halo'
xtitle= '!7l!dd!n!3'
hogg_scatterplot, pm[halo].pmdec,pm[halo].i, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' halo'

absmagi= 2.5
absmagdist= 0.01 ; kpc
distance= 1D1^(0.2*(pm.i-absmagi))*absmagdist ; kpc
mpykps= 4.73 ; 1 mas/yr at 1 kpc in km/s
vdec= pm.pmra*cos(pm.dec*!DPI/1.8D2)*distance*mpykps ; km/s
vra= pm.pmdec*distance*mpykps ; km/s
ytitle= 'height (kpc)'
yrange= [0.0,10.0]
xtitle= 'v!d!7a!3!n'
xrange= [-1,1]*200.0
hogg_scatterplot, vra[halo],distance[halo], $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' halo'
plot, vra[halo],distance[halo],psym=1, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  title=prefix+' halo'
xtitle= 'v!d!7d!3!n'
hogg_scatterplot, vdec[halo],distance[halo], $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  /nocontours,title=prefix+' halo'
plot, vdec[halo],distance[halo],psym=1, $
  xrange=xrange,xtitle=xtitle, $
  yrange=yrange,ytitle=ytitle, $
  title=prefix+' halo'

return
end
