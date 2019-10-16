;+
; NAME:
;  hogg_bell_analysis
; PURPOSE:
;  Perform Bell-like analysis on MW sims.
; BUGS:
;  - Magic number of 8.0 kpc hard-coded.
;-
pro hogg_bell_analysis
nside= 32
if (not keyword_set(infile)) then $
  infile='/global/data/scr/adi2/obsmw1_newamiga.fits'
obs= mrdfits(infile,1)
set_plot, 'ps'
xsize= 7.5 & ysize= xsize
device, /inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, bits=8
hogg_plot_defaults
shiftvec= [[obs.x-8.0],[obs.y],[obs.z]]
help, shiftvec
d2= total(shiftvec*shiftvec,2)
dd= sqrt(d2)
logd= alog10(dd)
theta= acos(shiftvec[*,2]/dd)
Dec= (0.5*!DPI-theta)*!radeg
phi= atan(shiftvec[*,1],shiftvec[*,0])
RA= phi*!radeg
negative= where(RA LT 0.0,nneg)
if (nneg GT 0) then RA[negative]= RA[negative]+360.0
splog, minmax(RA),minmax(Dec)
delta= 0.5
for inner=0.0,2.0,delta do begin
    innerstr= string(inner,format='(F3.1)')
    outerstr= string(inner+delta,format='(F3.1)')
    !P.TITLE= innerstr+' < log!d10!n(!8d!3 / kpc) < '+outerstr
    slice= where((logd GT float(innerstr)) AND (logd LT float(outerstr)))
    !X.RANGE=[-1,1]*80.0
    !Y.RANGE=!X.RANGE
    plot, obs[slice].x,obs[slice].z,psym=3, $
      xtitle='x (kpc)', $
      ytitle='z (kpc)'
    ang2pix_ring, nside,theta[slice],phi[slice],ind
    map=fltarr(12L*nside*nside)
    map[ind]++
    hogg_healpix_greyscale, map
endfor
device,/close
return
end
