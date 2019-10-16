;+
; NAME:
;  faintpm_verify
; PURPOSE:
;  Check a point-source hypothesis against an image.
; COMMENTS:
; INPUTS:
;  image   - 2-d image
;  invvar  - inverse variance estimates, matched to image
;  xx,yy   - [N] vectors of positions for a putative source
;  flux    - [N] vector of fluxes
;  fwhm    - PSF FWHM in pixels
; OUTPUTS:
;  [returns] - vector of delta-chi-squareds, matched to xx,yy
; BUGS:
;  - Could be sped up with knowledge of which calls differ only in
;    flux.
;  - Could be sped up by removing all divisions at the expense of
;    making the code significantly harder to understand.
; REVISION HISTORY:
;  2007-07-14  started by Hogg (NYU).
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_verify, image,invvar,xx,yy,flux,fwhm
tinvvar= total(invvar,/double)
if (tinvvar LE 0.0) then tinvvar= 1.0
timage= image-total(image*invvar,/double)/tinvvar ; subtract mean
sqrt8ln2= sqrt(8.0*alog(2.0))
foo= size(image,/dimens)
nx= foo[0]
ny= foo[1]
npos= n_elements(xx)
xpix= findgen(nx)#replicate(1.0,ny)
ypix= replicate(1.0,nx)#findgen(ny)
psfvar= (fwhm/sqrt8ln2)^2
deltachisq= 0.0*xx
for pp=0L,npos-1L do begin
    deltax= xpix-xx[pp]
    deltay= ypix-yy[pp]
    model= 1.0/(2*!PI*psfvar)*exp(-0.5*(deltax*deltax+deltay*deltay)/psfvar)
    model= model-total(model*invvar,/double)/tinvvar ; subtract mean
    deltachisq[pp]= total((flux[pp]*model-2.0*timage)*flux[pp]*model*invvar, $
                          /double)
endfor
return, deltachisq
end
