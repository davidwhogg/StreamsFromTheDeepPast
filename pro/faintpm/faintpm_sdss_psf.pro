;+
; NAME:
;  faintpm_sdss_psf
; PURPOSE:
;  measure the SDSS PSF FWHM for a field
; INPUTS:
;  flist   - flist structure (one element only, please)
;  x,y     - position of source in that field
; OUTPUTS:
;  [return value]  - PSF FWHM in pixels
; BUGS:
;  - Hard-coded upper and lower bounds on the PSF.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_sdss_psf, flist,xx,yy
sqrt8ln2= sqrt(8.0*alog(2.0))
fwhmlist= [2.0,3.5,5.0]          ; min, default, max
psfieldname= sdss_name('psField',flist.run,flist.camcol, $
                       flist.field,rerun=flist.rerun)
psfield= mrdfits(psfieldname,flist.filter+1,/silent)
npsfhalf= 2
npsf= 2*npsfhalf+1
xpsf= findgen(npsf)#replicate(1.0,npsf)-npsfhalf
ypsf= replicate(1.0,npsf)#findgen(npsf)-npsfhalf
psfimage= sdss_psf_recon(psfield,xx,yy,trimdim=[npsf,npsf])
psfimage= psfimage-mean(psfimage) ; remove mean level like in verify step
psfchisq= fltarr(3)
for ff=0,2 do begin
    psfvar= (fwhmlist[ff]/sqrt8ln2)^2
    model= exp(-0.5*(xpsf*xpsf+ypsf*ypsf)/psfvar)
    model= model-mean(model) ; remove mean level like in verify step
; perhaps we should really marginalize over amp?
    amp= total(psfimage*model,/double)/total(model*model,/double)
    dmodel= psfimage-amp*model
    psfchisq[ff]= total(dmodel*dmodel,/double)
endfor
invfwhmlist= 1.0/fwhmlist
foo= faintpm_parabola(invfwhmlist,psfchisq,invfwhmlist,invfwhm,bar,sigmafwhm)
invfwhm= invfwhm > min(invfwhmlist)
invfwhm= invfwhm < max(invfwhmlist)
fwhm= 1.0/invfwhm
if (sigmafwhm LE 0.0) then begin
    splog, 'something went wrong with PSF FWHM determination for ' $
      +psfieldname
    splog, fwhm,sigmafwhm
    splog, psfimage
    splog, fwhmlist
    splog, psfchisq
    splog, 'using min value'
    fwhm= min(fwhmlist) ; typically failures are images with small PSF
endif
splog, 'fwhm',fwhm,' for '+psfieldname
return, fwhm
end
