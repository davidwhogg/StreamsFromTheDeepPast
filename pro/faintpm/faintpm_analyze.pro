;+
; NAME:
;  faintpm_analyze
; PURPOSE:
;  Make plots for the purpose of analyzing a proper-motion
;    determination.
; INPUTS:
;  fitsfile     - multi-epoch fits file like that made by
;                 faintpm_coadd_catalog
;  radot,decdot - proper motion to test / analyze / compare to zero (mas/yr)
; BUGS:
;  - Won't work on sets of images with different sizes or pixel scales.
;  - Won't work near the celestial poles.
;  - Requires various header keywords (read the source, Luke).
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_analyze, fitsfile,radot,decdot
prefix= strmid(fitsfile,0,strpos(fitsfile,'.fits'))
splog, prefix
smosaic_make_jpg_scales, scales,nonlinearity
scales= 3*[scales[1],scales[1],scales[1]]
quality= 95

; read the files and figure out the epoch
hdr= headfits(fitsfile)
nepoch= sxpar(hdr,'NEPOCH')
mjd= fltarr(nepoch)
for ii=0,nepoch-1 do begin
    exten= 2*ii+1
    thishdr= headfits(fitsfile,exten=exten)
    extast, thishdr,thisastr
    if (not keyword_set(astr)) then begin
        astr= replicate(thisastr,nepoch)
    endif
    astr[ii]= thisastr
    mjd[ii]= sxpar(thishdr,'MJD')
    if (not keyword_set(image)) then begin
        nx= sxpar(thishdr,'NAXIS1')
        ny= sxpar(thishdr,'NAXIS2')
        image= fltarr(nx,ny,nepoch)
    endif
    image[*,*,ii]= mrdfits(fitsfile,exten)
endfor
epoch= mean(mjd)

; find the image closest to the epoch and grab the WCS
foo= min(abs(mjd-epoch),ii)
bigastr= astr[ii]
margin= 2 ; pixels
bigastr= struct_addtags(bigastr,{naxis:[nx+2*margin,ny+2*margin]})
bigastr.crpix= bigastr.crpix+[margin,margin]

; stack all the images at zero proper motion and with the proper
; motion compensated
; NB: you need to subtract the motion from the RA,Dec of the image (CRVAL)
for kk=0,1 do begin
    numerator= fltarr(bigastr.naxis)
    denominator= numerator+1e-9 ; tiny
    if (kk EQ 0) then medianblock= fltarr([bigastr.naxis,nepoch])
    for ii=0,nepoch-1 do begin
        thisimage= reform(image[*,*,ii],nx,ny)
        weight= fltarr(nx,ny)+1.0
        thisastr= astr[ii]
        if (kk EQ 1) then thisastr.crval= thisastr.crval $
          -[radot,decdot]*(mjd[ii]-epoch) $
          *(0.001/365.25/3600.0) ; mas/yr conversion
        smosaic_remap, thisimage,thisastr,bigastr,outimage,offset, $
          weight=weight,outweight=outweight,reflimits=reflim, $
          outlimits=outlim
        numerator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]]= $
          numerator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]] $
          +outweight[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]] $
          *outimage[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]]
        denominator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]]= $
          denominator[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1]] $
          +outweight[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]]
        if (kk EQ 0) then $
          medianblock[reflim[0,0]:reflim[0,1],reflim[1,0]:reflim[1,1],ii]= $
          outimage[outlim[0,0]:outlim[0,1],outlim[1,0]:outlim[1,1]]
    endfor
    thisimage= numerator/denominator
    thishdr= hdr
    putast, thishdr,bigastr
    thisfilename= prefix+'_zero.fits'
    if (kk EQ 1) then thisfilename= prefix+'_pm.fits'
    mwrfits, thisimage,thisfilename,thishdr,/create
    mwrfits, denominator,thisfilename
    thisjpgname= prefix+'_zero.jpg'
    if (kk EQ 1) then thisjpgname= prefix+'_pm.jpg'
    nw_rgb_make, thisfilename,thisfilename,thisfilename, $
      name=thisjpgname,scales=scales,nonlinearity=nonlinearity, $
      quality=quality,rebinfactor=0.25
    if (kk EQ 0) then begin
        tnx= bigastr.naxis[0]
        tny= bigastr.naxis[1]
        thisimage= median(medianblock,dimension=3)
        help, thisimage
        thisfilename= prefix+'_median.fits'
        mwrfits, thisimage,thisfilename,thishdr,/create
        thisjpgname= prefix+'_median.jpg'
        nw_rgb_make, thisfilename,thisfilename,thisfilename, $
          name=thisjpgname,scales=scales,nonlinearity=nonlinearity, $
          quality=quality,rebinfactor=0.25
        stop
    endif
endfor

return
end
