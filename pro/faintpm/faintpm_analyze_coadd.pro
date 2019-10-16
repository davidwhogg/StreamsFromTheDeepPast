;+
; PURPOSE:
;  Run faintpm_analyze on coadd results.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_analyze_coadd
doterrlim= 60.0 ; mas/yr
prefix= 'fermi_coadd_catalog'
orig= mrdfits(prefix+'.fits',1)
cc= mrdfits(prefix+'.fits',2)
ncc= n_elements(cc)
pmmag=sqrt(cc.radot^2+cc.decdot^2)
for ii=0,ncc-1L do begin
    thiscc= cc[ii]
    if ((thiscc.fluxerr GT 0.0) AND $
        (thiscc.raerr GT 0.0) AND $
        (thiscc.decerr GT 0.0) AND $
        (thiscc.radoterr GT 0.0) AND $
        (thiscc.radoterr LT doterrlim) AND $
        (thiscc.decdoterr GT 0.0) AND $
        (thiscc.decdoterr LT doterrlim) AND $
        (pmmag[ii] GT 300.0) AND $
        (thiscc.flux GT (5.0*thiscc.fluxerr))) then begin
        splog, 'going to look at '+thiscc.filename
        if file_test(thiscc.filename) then begin
            faintpm_analyze, thiscc.filename,thiscc.radot,thiscc.decdot
        endif
    endif
endfor
return
end
