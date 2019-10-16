;+
; NAME:
;  faintpm_coadd_plot
; PURPOSE:
;  Plot proper motions determined by faintpm_coadd_catalog.
; BUGS:
;  - Need galactic coordinates?
; LICENSE:
;  Copyright 2008 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_coadd_plot
prefix= 'fermi_coadd_catalog'
fitsfilename= prefix+'.fits'
doterrlim= 50.0                 ; mas/yr
band= 4

; read files
fcc= mrdfits(prefix+'.fits*',1)
ncc= n_elements(fcc)
hogg= mrdfits(prefix+'.fits*',2)
pmmag= sqrt(hogg.radot^2+hogg.decdot^2)

; classify objects badness
badness= bytarr(ncc)
bad= where((hogg.fluxerr LE 0.0) OR $
           (hogg.raerr LE 0.0) OR $
           (hogg.decerr LE 0.0) OR $
           (hogg.radoterr LE 0.0) OR $
           (hogg.decdoterr LE 0.0),nbad)
if (nbad GT 0) then badness[bad]= badness[bad]+4
boring= where((hogg.radoterr GT doterrlim) OR $
              (hogg.decdoterr GT doterrlim),nboring)
if (nboring GT 0) then badness[boring]= badness[boring]+1
crappy= where((hogg.flux LT (3.0*hogg.fluxerr)),ncrappy)
if (ncrappy GT 0) then badness[crappy]= badness[crappy]+2
good= where(badness EQ 0,ngood)

; output BD candidates for SJ
bdcandidate= where((badness EQ 0) AND $
                   (pmmag GT 250.0),nbd)
if (nbd GT 0) then begin
    bdfile= 'hogg_for_sj.fits'
    mwrfits, fcc[bdcandidate],bdfile,/create
    mwrfits, hogg[bdcandidate],bdfile
endif

; match to known BDs and QSOs
readcol, '../../data/faintpm/Jiang_LTdata.dat',foo,bar,format='A,A',comment='#'
sdss= where((foo EQ 'SDSS'))
bdname= foo[sdss]+' '+bar[sdss]
nbd= n_elements(bdname)
bdra= dblarr(nbd)
bddec= dblarr(nbd)
for jj=0L,nbd-1 do begin
    hogg_name2radec, bdname[jj],tra,tdec
    bdra[jj]= tra
    bddec[jj]= tdec
endfor
readcol, '../../data/faintpm/Jiang_*quasar*.dat',foo,format='A',comment='#'
qsoname= 'SDSS '+foo
nqso= n_elements(qsoname)
qsora= dblarr(nqso)
qsodec= dblarr(nqso)
for jj=0L,nqso-1 do begin
    hogg_name2radec, qsoname[jj],tra,tdec
    qsora[jj]= tra
    qsodec[jj]= tdec
endfor
spherematch, bdra,bddec,fcc[good].ra,fcc[good].dec,2.0/3600.0, $
  match1,bdindx,dist12
nbd= 0
if (match1[0] NE -1) then begin
    nbd= n_elements(match1)
    bdindx= good[bdindx]
    bdname= bdname[match1]
    bdra= bdra[match1]
    bddec= bddec[match1]
endif
spherematch, qsora,qsodec,fcc[good].ra,fcc[good].dec,2.0/3600.0, $
  match1,qsoindx,dist12
nqso= 0
if (match1[0] NE -1) then begin
    nqso= n_elements(match1)
    qsoindx= good[qsoindx]
    qsoname= qsoname[match1]
    qsora= qsora[match1]
    qsodec= qsodec[match1]
endif
help, nbd,nqso

; setup plot
set_plot,'ps'
hogg_plot_defaults
xsize= 6.0
ysize= xsize
deltamargin= [2,-2]
!X.OMARGIN= !X.OMARGIN+deltamargin
!P.TITLE= 'all well-measured i-dropouts in the Fermi Coadd Catalog'

; make pm plot
device, filename=prefix+'_pm_pm.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches,/color
pmmax= 1000.0
pmrange= [-1,1]*pmmax
plot, hogg[good].radot,hogg[good].decdot,psym=6,symsize=0.1, $
  xrange=pmrange,xtitle='d!8RA!3 / d!8t!3 (mas / yr)', $
  yrange=pmrange,ytitle='d!8Dec!3 / d!8t!3 (mas / yr)'
hoggscale= 1.0
hogg_usersym, 5,/stellar,scale=hoggscale,/fill
oplot, [hogg[bdindx].radot],[hogg[bdindx].decdot],psym=8
hogg_usersym, 12,scale=hoggscale,/fill
oplot, [hogg[qsoindx].radot],[hogg[qsoindx].decdot],psym=8
device,/close

; make pm vs mag plot
device, filename=prefix+'_pm_mag.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches,/color
pmrange= [0,1]*pmmax
zrange= [20,21]
plot, fcc[good].psfcounts[band],pmmag[good],psym=6,symsize=0.1, $
  xrange=zrange,xtitle='!8z!3 (mag)', $
  yrange=pmrange,ytitle='proper motion magnitude (mas / yr)'
hogg_usersym, 5,/stellar,scale=hoggscale,/fill
oplot, [fcc[bdindx].psfcounts[band]],[pmmag[bdindx]],psym=8
hogg_usersym, 12,scale=hoggscale,/fill
oplot, [fcc[qsoindx].psfcounts[band]],[pmmag[qsoindx]],psym=8
device,/close

return
end
