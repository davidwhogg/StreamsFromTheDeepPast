;+
; INPUTS:
;  RA,Dec       - position at epoch (deg)
;  radot,decdot - time derivatives (mas/yr)
;  epoch        - time of applicability of RA,Dec (MJD)
;  filename     - name for the fits file for the stack
; OPTIONAL INPUTS:
;  run          - list of runs to consider (otherwise consider all
;                 photometric runs)
; BUGS:
;  - No proper header.
;  - Cutout size hard-coded.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_sdss_cutout, ra,dec,radot,decdot,epoch,filename,run=run
nhalf= 7

; get SDSS information
flist= faintpm_sdss_flist(ra,dec,run=run)
nfield= n_elements(flist)
if tag_exist(flist,'RUN') then begin
    mjd= sdss_run2mjd(flist.run)
    splog, 'mean mjd for this source:',mean(mjd)
    nepoch= 0
endif else begin
    nfield= 0
    nepoch= -1
endelse

; prepare material for fits file
pixscale= 0.396/3600.0          ; deg
dra= pixscale*(2.0*nhalf+1.0)   ; deg
defaultastr= hogg_make_astr(ra,dec,dra,dra,pixscale=pixscale,npixround=1)
dummy= fltarr(2,2)
mwrfits, dummy,filename,/create

oldrun= -1
for rr=0,nfield-1 do if (flist[rr].run NE oldrun) then begin
    thisdt= (mjd[rr]-epoch)/365.25 ; yr
    thisra= ra+radot*thisdt/3.6D6 ; mas -> deg
    thisdec= dec+decdot*thisdt/3.6D6 ; mas -> deg

; compute x,y position in field
    thisimage= smosaic_cacheimg(flist[rr],invvar=thisinvvar, $
                                /globalsky,ncache=1,/dontcrash)
    foo= size(thisimage,/dimens)
    thisnx= foo[0]
    thisny= foo[1]
    thisgsa= sdss_astrom(flist[rr].run,flist[rr].camcol, $
                         flist[rr].field,rerun=flist[rr].rerun, $
                         filter=flist[rr].filter)
    astrom_adxy, thisgsa,thisra,thisdec,xpix=thisx,ypix=thisy
    thisx= thisx[0]
    thisy= thisy[0]
    xint= round(thisx)
    yint= round(thisy)

; are we inside the image?
    if ((xint GE (nhalf)) AND (xint LT (thisnx-nhalf-1)) AND $
        (yint GE (nhalf)) AND (yint LT (thisny-nhalf-1))) then begin
        splog, flist[rr].run
        oldrun= flist[rr].run

; compute astrometric derivatives
        astrom_xyad, thisgsa,double(xint),double(yint),ra=crval0,dec=crval1
        thisastr= defaultastr
        thisastr.crval= [crval0,crval1]
        dx= 1D0
        astrom_xyad, thisgsa,double(xint)+dx,double(yint),ra=rax,dec=decx
        dradx= (rax-crval0)/dx
        ddecdx= (decx-crval1)/dx
        dy= dx
        astrom_xyad, thisgsa,double(xint),double(yint)+dy,ra=ray,dec=decy
        drady= (ray-crval0)/dy
        ddecdy= (decy-crval1)/dy
        thisastr.cd= transpose([[dradx,drady],[ddecdx,ddecdy]])

; make cutout
        thisimage= thisimage[xint-nhalf:xint+nhalf,yint-nhalf:yint+nhalf]
        thisinvvar= thisinvvar[xint-nhalf:xint+nhalf,yint-nhalf:yint+nhalf]

; write images
        mkhdr, thishdr,thisimage
        sxaddpar, thishdr,'MJD',mjd[rr]
        sxaddpar, thishdr,'FWHM',faintpm_sdss_psf(flist[rr],thisx,thisy)
        putast, thishdr,thisastr
        mwrfits, thisimage,filename,thishdr
        mwrfits, thisinvvar,filename
        nepoch= nepoch+1

; test FITS file
        imageerror= thisimage-mrdfits(filename,nepoch*2-1,/silent)
        invvarerror= thisinvvar-mrdfits(filename,nepoch*2,/silent)
        if ((max(abs(imageerror)) GT 0.0) OR $
            (max(abs(invvarerror)) GT 0.0)) then begin
            splog, 'ERROR: FITS file corrupt'
            stop
        endif
    endif
endif

; write nepoch into fits file and finish
splog, 'wrote',nepoch,' epochs'
mkhdr, hdr,dummy
sxaddpar, hdr,'NEPOCH',nepoch
djs_modfits, filename,dummy,hdr
return
end
