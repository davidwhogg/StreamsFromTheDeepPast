;+
; NAME:
;  faintpm_verify_grid
; PURPOSE:
;  Run proper-motion verification on a grid, to measure proper
;  motions.
; INPUTS:
;  ra     - [nra] vector of RAs to try (deg)
;  dec    - [ndec] vector of Decs to try (deg)
;  radot  - [nradot] vector of dRA/dts to try (mas/yr)
;  decdot - [ndecdot] vector of dDec/dts to try (mas/yr)
;  flux   - [nflux] vector of fluxes to consider (assumes images are
;           from the same band and calibrated the same)
;  epoch  - time of validity of (ra,dec) in same units as tname
;  tname  - name of the mean time of observation in FITS headers
;  file   - multi-EXTEN FITS file (calibrated and with astrometry
;           and mean times of observation); read the source for
;           documentation on this file.
; OPTIONAL INPUTS:
;  tunit  - multiplier to get tname time into years (eg, 1.0/365.25 if
;           tname time is in days); default 1.0
;  priorvec - information for a set of gaussian priors; read the
;             source for format; default no prior
; OUTPUTS:
;  deltachisq  - [nra,ndec,nradot,ndecdot,nflux] grid of
;                delta-chi-squared values
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
function faintpm_verify_grid, ra,dec,radot,decdot,flux,epoch,tname,file, $
                              tunit=tunit,priorvec=priorvec

; get dimensions
nra= n_elements(ra)
ndec= n_elements(dec)
nradot= n_elements(radot)
ndecdot= n_elements(decdot)
nflux= n_elements(flux)
deltachisq= fltarr(nra,ndec,nradot,ndecdot,nflux)
hdr= headfits(file)
nepoch= sxpar(hdr,'NEPOCH')
splog, 'found',nepoch,' epochs'

; make useful matrices
ralist= reform(ra#replicate(1D0,ndec*nradot*ndecdot*nflux), $
               nra,ndec,nradot,ndecdot,nflux)
declist= reform(reform(replicate(1D0,nra)#dec, $
                       nra*ndec) $
                #replicate(1.0,nradot*ndecdot*nflux), $
                nra,ndec,nradot,ndecdot,nflux)
radotlist= reform(reform(replicate(1D0,nra*ndec)#radot, $
                         nra*ndec*nradot) $
                  #replicate(1.0,ndecdot*nflux), $
                  nra,ndec,nradot,ndecdot,nflux)
decdotlist= reform(reform(replicate(1.0,nra*ndec*nradot)#decdot, $
                          nra*ndec*nradot*ndecdot) $
                   #replicate(1.0,nflux), $
                   nra,ndec,nradot,ndecdot,nflux)
fluxlist= reform(replicate(1.0,nra*ndec*nradot*ndecdot)#flux, $
                 nra,ndec,nradot,ndecdot,nflux)

; add in priors
if keyword_set(priorvec) then begin
    deltachisq= deltachisq+((ralist-priorvec[0])/priorvec[1])^2
    deltachisq= deltachisq+((declist-priorvec[2])/priorvec[3])^2
    deltachisq= deltachisq+((radotlist-priorvec[4])/priorvec[5])^2
    deltachisq= deltachisq+((decdotlist-priorvec[6])/priorvec[7])^2
    deltachisq= deltachisq+((fluxlist-priorvec[8])/priorvec[9])^2
endif

; loop over images
for ff=0,nepoch-1 do begin
    exten= 2*ff+1
    image= mrdfits(file,exten,hdr,/silent)
    invvar= mrdfits(file,exten+1,/silent)
    extast, hdr,astr
    tt= sxpar(hdr,tname)
    fwhm= sxpar(hdr,'FWHM')

; deal with coordinates
    thisra= ralist+radotlist*(tt-epoch)*tunit/3.6D6 ; mas -> deg
    thisdec= declist+decdotlist*(tt-epoch)*tunit/3.6D6 ; mas -> deg
    ad2xy, thisra,thisdec,astr,thisxx,thisyy

; verify
    thischisq= faintpm_verify(image,invvar,thisxx,thisyy,fluxlist,fwhm)
    if min(finite(thischisq)) LT 1 then stop

; end loop over images
    deltachisq= deltachisq+thischisq
endfor
return, deltachisq
end
