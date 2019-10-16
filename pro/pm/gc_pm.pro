;+
; NAME:
;   gc_pm.pro
; BUGS:
;   - Do I really need to multiply by cos Dec?
;   - Header not written.
;-
pro gc_pm

if 0 then begin
    prefix= 'ngc6205m13'
    ra0= 250.423                ; deg J2000
    dec0= 36.460                ; deg J2000
    rad1= 0.2                   ; deg
endif

if 1 then begin
    prefix= 'pal5'
    ra0= 229.022                ; deg J2000
    dec0= -0.111                ; deg J2000
    rad1= 3.0/60.0              ; deg
endif

rad2= 2.0
seed= -1L

pm= mrdfits('../../data/pm_'+prefix+'*.fit*',1)
pm= pm[where((pm.pmDecErr lt 10.0) AND (pm.i gt 16.5))]
nstar= n_elements(pm)
pm.pmRA= pm.pmRA+randomu(seed,nstar)
pm.pmDec= pm.pmDec+randomu(seed,nstar)

adist= djs_diff_angle(pm.ra,pm.dec,ra0,dec0)
near= where(adist LT rad1)
far= where((adist GT rad1) AND (adist LT rad2))
small= !DPI*rad1^2
big= !DPI*rad2^2-small
ratio= small/big
splog, 'ratio',ratio

set_plot, 'ps'
hogg_plot_defaults
xsize= 7.5 & ysize= 7.5
psfilename= prefix+'.ps'
device, file=psfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,bits=8
!P.MULTI= [0,3,3]
!X.MARGIN= !X.OMARGIN
!X.OMARGIN= [0,0]
!Y.MARGIN= !Y.OMARGIN
!Y.OMARGIN= [0,0]

for pp=0,2 do begin

    if (pp EQ 0) then begin
        xx= pm.ra
        xrange= 0
        xtitle= 'RA (deg) J2000'
        xnpix= 0
        yy= pm.dec
        yrange= 0
        ytitle= 'Dec (deg) J2000'
        ynpix= 0
    endif else if (pp EQ 1) then begin
        xx= pm.g-pm.i
        xrange= [0.0,1.0]
        xtitle= 'g - i'
        xnpix= 100
        yy= pm.i
        yrange= [19.5,16.5]
        ytitle= 'i'
        ynpix= xnpix
    endif else begin
        xx= pm.pmRA*cos(pm.Dec*!DPI/1.8D2)
        xrange= [-1,1]*15.0
        xtitle= '!7l!da!n!3 cos !7d!3'
        xnpix= 30
        yy= pm.pmDec
        yrange= xrange
        ytitle= '!7l!dd!n!3'
        ynpix= xnpix
    endelse

    hogg_scatterplot, xx[far],yy[far], $
      xnpix=xnpix,xrange=xrange,xtitle=xtitle,xvec=xvec, $
      ynpix=ynpix,yrange=yrange,ytitle=ytitle,yvec=yvec, $
      grid=fargrid,/nocontours,title=prefix+': far stars'
    hogg_scatterplot, xx[near],yy[near], $
      xnpix=xnpix,xrange=xrange,xtitle=xtitle, $
      ynpix=ynpix,yrange=yrange,ytitle=ytitle, $
      grid=neargrid,/nocontours,title=prefix+': near stars'

    diffgrid= neargrid-fargrid*ratio

    hogg_scatterplot, xx,yy, $
      xnpix=xnpix,xrange=xrange,xtitle=xtitle, $
      ynpix=ynpix,yrange=yrange,ytitle=ytitle, $
      grid=diffgrid,/usegrid,/nocontours,title=prefix+': near - far'

endfor 

; oplot, [1,1]*(-0.90-1.00),!Y.CRANGE,psym=0,color=djs_icolor('red')
; oplot, [1,1]*(-0.90+1.00),!Y.CRANGE,psym=0,color=djs_icolor('red')
; oplot, !X.CRANGE,[1,1]*( 5.50-1.00),psym=0,color=djs_icolor('red')
; oplot, !X.CRANGE,[1,1]*( 5.50+1.00),psym=0,color=djs_icolor('red')

weight= fltarr(nstar)+1.0
weight[far]= -ratio
oplot, [1,1]*weighted_quantile(xx,weight),!Y.CRANGE,psym=0
oplot, !X.CRANGE,[1,1]*weighted_quantile(yy,weight),psym=0
oplot, [1,1]*weighted_quantile(xx[near]),!Y.CRANGE,psym=0,linestyle=1
oplot, !X.CRANGE,[1,1]*weighted_quantile(yy[near]),psym=0,linestyle=1

pmq= mrdfits('../../data/pm_'+prefix+'*qso*.fit*',1)
pmq= pmq[where(pmq.pmdec lt 10.0)]
nqso= n_elements(pmq)
pmq.pmRA= pmq.pmRA+randomu(seed,nqso)
pmq.pmDec= pmq.pmDec+randomu(seed,nqso)

xx= pmq.pmRA*cos(pmq.Dec*!DPI/1.8D2)
yy= pmq.pmDec

hogg_scatterplot, xx,yy, $
  xnpix=xnpix,xrange=xrange,xtitle=xtitle,xvec=xvec, $
  ynpix=ynpix,yrange=yrange,ytitle=ytitle,yvec=yvec, $
  grid=qsogrid,/nocontours,title=prefix+': QSOs'
oplot, [1,1]*weighted_quantile(xx),!Y.CRANGE,psym=0
oplot, !X.CRANGE,[1,1]*weighted_quantile(yy),psym=0

return
end
