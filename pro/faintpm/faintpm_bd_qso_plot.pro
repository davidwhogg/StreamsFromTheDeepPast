;+
; NAME:
;  faintpm_bd_qso_plot
; PURPOSE:
;  Compare PM measurements for BDs and QSOs.
; BUGS:
;  - Not written.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_bd_qso_plot
prefix= 'faintpm_bd_qso'
set_plot,'ps'
hogg_plot_defaults
xsize= 6.0
ysize= xsize
device, filename=prefix+'.ps',xsize=xsize,xoffset=0.5*(8.5-xsize), $
  ysize=ysize,yoffset=0.5*(11.0-ysize), $
  /inches,/color
!X.TITLE= 'approximate !8z!3-band magnitude'
!X.RANGE= [18,21]
!Y.TITLE= 'proper motion (mas/yr)'
!Y.RANGE= [0,350]
size= 2.0*!P.SYMSIZE

zmag= !X.RANGE[0]+(!X.RANGE[1]-!X.RANGE[0])*findgen(101)/100.0
pm= 40.0*10.0^(-0.4*(20.5-zmag)) ; rough estimate of 1-sigma error
plot, zmag,pm
oplot, zmag,3.0*pm,thick=3.0*!P.THICK
pm= 200*10.0^(0.2*(20.5-zmag)) ; proper-motion scaling
oplot, zmag,pm,linestyle=2
maglimit= 20.5-2.5*alog10(1.5) ; approximation to 3-sigma SDSS limit
oplot, [1,1]*maglimit,!Y.CRANGE,thick=3.0*!P.THICK

filename= ['faintpm_bd.fits','faintpm_highz.fits']
label= ['BD','QSO']
nside= [5,12]
stellar= [1,0]
nfile= n_elements(filename)

for ii=0,nfile-1 do begin
    qso= mrdfits(filename[ii],1)
    good= where((qso.flux GT 0.0) AND $
                (qso.fluxerr GT 0.0) AND $
                (qso.raerr GT 0.0) AND $
                (qso.decerr GT 0.0) AND $
                (qso.radoterr GT 0.0) AND $
                (qso.radoterr LT 50.0) AND $
                (qso.decdoterr GT 0.0) AND $
                (qso.decdoterr LT 50.0) AND $
                (qso.nepoch GE 6))
    help, good
    qso= qso[good]
    zmag= 22.5-2.5*alog10(qso.flux)-0.2 ; aperture mag hack
    cosdec= cos(!PI*qso.dec/180.0)
    pm= sqrt((qso.radot/cosdec)^2+qso.decdot^2)
    splog, minmax(pm)
    hogg_usersym, nside[ii],size=size,stellar=stellar[ii],/fill
    oplot, zmag,pm,psym=8
    xlabel= !X.CRANGE[0]+(!X.CRANGE[1]-!X.CRANGE[0])*0.05
    dylabel= (!Y.CRANGE[1]-!Y.CRANGE[0])*0.025
    ylabel= !Y.CRANGE[0]+dylabel*(38-1.5*ii)
    oplot, [xlabel],[ylabel+0.2*dylabel],psym=8
    xyouts, xlabel,ylabel,' '+label[ii]
endfor

device,/close
cmd= 'pstopnm --stdout --portrait --dpi=200 --xborder=0 --yborder=0 '+prefix+'.ps | pamscale 0.5 | pnmtojpeg --quality=90 > '+prefix+'.jpg'
splog, cmd
spawn, cmd
return
end
