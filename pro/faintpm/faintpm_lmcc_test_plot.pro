;+
; BUGS:
;  - No proper header.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_lmcc_test_plot
ddmy= 365.25*3.6D6 ; conversion factor
prefix= 'faintpm_lmcc'

; read structures
filename= prefix+'.fits'
lmcc= mrdfits(filename,1)
hogg= mrdfits(filename,2)
good= where((hogg.fluxerr GT 0.0) AND $
            (hogg.raerr GT 0.0) AND $
            (hogg.decerr GT 0.0) AND $
            (hogg.radoterr GT 0.0) AND $
            (hogg.decdoterr GT 0.0) AND $
            (hogg.radoterr LT 100.0) AND $
            (hogg.decdoterr LT 100.0))
lmcc= lmcc[good]
hogg= hogg[good]

; setup plot
set_plot, 'ps'
hogg_plot_defaults
xsize= 6.0
ysize= xsize
device, filename=prefix+'.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches,/color
size= 2.0*!P.SYMSIZE
!X.OMARGIN= !X.OMARGIN+[1,-1]
!X.RANGE= [-350,350]
!X.TITLE= 'LMCC proper motion (mas/yr)'
!Y.RANGE= !X.RANGE
!Y.TITLE= 'Hogg proper motion (mas/yr)'
charsize= 0.75*!P.CHARSIZE

!P.TITLE= 'proper motions along RA direction'
plot, !X.RANGE,!X.RANGE
hogg_usersym, 12,/fill,size=size
oplot, lmcc.ra_pm*ddmy,hogg.radot,psym=8
djs_oploterr, lmcc.ra_pm*ddmy,hogg.radot, $
  xerr=lmcc.ra_pm_err*ddmy, $
  yerr=hogg.radoterr
bad= where(abs(lmcc.ra_pm*ddmy-hogg.radot) GT 150,nbad)
for bb=0L,nbad-1L do begin
    xyouts, lmcc[bad[bb]].ra_pm*ddmy,hogg[bad[bb]].radot, $
      ' '+hogg_iau_name(lmcc[bad[bb]].ra_mean,lmcc[bad[bb]].dec_mean)+' ', $
      align= float(lmcc[bad[bb]].ra_pm*ddmy LT hogg[bad[bb]].radot), $
      charsize=charsize
endfor

!P.TITLE= 'proper motions along Dec direction'
plot, !X.RANGE,!X.RANGE
oplot, lmcc.dec_pm*ddmy,hogg.decdot,psym=8
djs_oploterr, lmcc.dec_pm*ddmy,hogg.decdot, $
  xerr=lmcc.dec_pm_err*ddmy, $
  yerr=hogg.decdoterr
bad= where(abs(lmcc.dec_pm*ddmy-hogg.decdot) GT 150,nbad)
for bb=0L,nbad-1L do begin
    xyouts, lmcc[bad[bb]].dec_pm*ddmy,hogg[bad[bb]].decdot, $
      ' '+hogg_iau_name(lmcc[bad[bb]].dec_mean,lmcc[bad[bb]].dec_mean)+' ', $
      align= float(lmcc[bad[bb]].dec_pm*ddmy LT hogg[bad[bb]].decdot), $
      charsize=charsize
endfor

device, /close
return
end
