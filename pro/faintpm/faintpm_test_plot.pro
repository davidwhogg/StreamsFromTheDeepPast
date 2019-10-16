;+
; BUGS:
;  - No proper header.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_test_plot
prefix= 'faintpm_test_plot'

; read and arrange data
tstr= mrdfits('faintpm_test.fits*',1)
ntstr= n_elements(tstr)
badness= bytarr(ntstr)
theorydoterr= tstr.meanfwhm/(tstr.totalsnr*tstr.sigmat)*3.6E6 ; mas
raokay= where(tstr.radoterr^2 GT theorydoterr^2,nraokay)
decokay= where(tstr.decdoterr^2 GT theorydoterr^2,ndecokay)
bad= where((tstr.decdoterr^2 GT (2.0*theorydoterr)^2) OR $
           (tstr.decdoterr^2 GT (2.0*theorydoterr)^2),nbad)
terrible= where((tstr.radotunc LT 0.0) OR $
                (tstr.decdotunc LT 0.0),nterrible)
if (nraokay GT 0) then badness[raokay]= badness[raokay]+1
if (ndecokay GT 0) then badness[decokay]= badness[decokay]+2
if (nbad GT 0) then badness[bad]= badness[bad]+4
if (nterrible GT 0) then badness[terrible]= badness[terrible]+8
splog, 'RA > 1 sigma:',float(nraokay)/float(ntstr)
splog, 'Dec > 1 sigma:',float(ndecokay)/float(ntstr)
splog, 'something > 2 sigma:',float(nbad)/float(ntstr)
splog, 'something hit chisq edge:',float(nterrible)/float(ntstr)

set_plot,'ps'
hogg_plot_defaults
xsize= 6.0
ysize= xsize
device, filename=prefix+'.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches,/color
!X.TITLE= 'total !8s/n!3 in all epochs combined'
!X.RANGE= [2,30]
!Y.TITLE= 'number !8N!3 of epochs'
!Y.RANGE= [4,51]
nside= 12
size= 2.0*!P.SYMSIZE
plot, tstr.totalsnr,tstr.nepoch,/nodata

; plot 3-sigma line
nepochline= findgen(100)+1
snrline= 3.0*sqrt(nepochline)
oplot, snrline,nepochline,thick=3.0*!P.THICK

; plot circles
this= where(badness GE 4,nthis)
if (nthis GT 0) then begin
    hogg_usersym, nside,/fill,color=djs_icolor('black'),size=size
    oplot, tstr[this].totalsnr,tstr[this].nepoch,psym=8
endif
this= where((badness GT 1) AND (badness LT 4),nthis)
if (nthis GT 0) then begin
    hogg_usersym, nside,/fill,color=djs_icolor('grey'),size=size
    oplot, tstr[this].totalsnr,tstr[this].nepoch,psym=8
endif
this= where((badness EQ 0),nthis)
if (nthis GT 0) then begin
    hogg_usersym, nside,/fill,color=djs_icolor('white'),size=size
    oplot, tstr[this].totalsnr,tstr[this].nepoch,psym=8
endif
hogg_usersym, nside,thick=!P.THICK,size=size
oplot, tstr.totalsnr,tstr.nepoch,psym=8
device,/close

cmd= 'pstopnm --stdout --portrait --dpi=200 --xborder=0 --yborder=0 '+prefix+'.ps | pamscale 0.5 | pnmtojpeg --quality=90 > '+prefix+'.jpg'
splog, cmd
spawn, cmd
return
end
