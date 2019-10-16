;+
; NAME:
;  faintpm_makefake
; PURPOSE:
;  Make fake data for faint-source proper motion studies.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_makefake, nreal

; set defaults
if (not keyword_set(nreal)) then nreal= 2
sqrt8ln2= sqrt(8.0*alog(2.0))
seed= -1L
raref= 180.0
decref= 0.0
pixscale= 0.25/3600.0 ; deg
dra= pixscale*16.0
ddec= pixscale*16.0
invvar= 1.0 ; inverse variance for each pixel

; begin loop over realizations
for rr=0,nreal-1 do begin
    rrstr= string(rr,format='(I3.3)')
    prefix= 'fake_'+rrstr

; make epochs and seeing
    nepoch= floor(5.0+46.0*randomu(seed))
    tt= 5.0*randomu(seed,nepoch)
    meant= mean(tt)
    fwhm= 3.0+1.0*randomu(seed,nepoch) ; pix
    snr= fltarr(nepoch)

; make astrometric headers
    for ii=0,nepoch-1 do begin
        thisorientation=   (randomu(seed)*2.0-1.0)*45.0
        thisraref=  raref+ (randomu(seed)*2.0-1.0)*2.0*pixscale
        thisdecref= decref+(randomu(seed)*2.0-1.0)*2.0*pixscale
        thisastr= hogg_make_astr(thisraref,thisdecref,dra,ddec, $
                                 pixscale=pixscale, $
                                 orientation=thisorientation)
        if (ii EQ 0) then astr= replicate(thisastr,nepoch)
        astr[ii]= thisastr
    endfor

; make random position and proper motion
    dramean=  (randomu(seed)*2.0-1.0)*1.0*pixscale ; deg
    ddecmean= (randomu(seed)*2.0-1.0)*1.0*pixscale ; deg
    radot=    (randomu(seed)*2.0-1.0)*1.0*pixscale*3.6D6 ; mas/yr
    decdot=   (randomu(seed)*2.0-1.0)*1.0*pixscale*3.6D6 ; mas/yr

; make noise-free images
    nx= astr[0].naxis[0]
    ny= astr[0].naxis[1]
    image= fltarr(nx,ny,nepoch)
    xx= findgen(nx)#replicate(1.0,ny)
    yy= replicate(1.0,nx)#findgen(ny)
    flux= (10.0+150.0*randomu(seed))/sqrt(nepoch) ; units of per-pixel sigma
    for ii=0,nepoch-1 do begin
        ad2xy, raref +dramean +radot *(tt[ii]-meant)/3.6D6, $
               decref+ddecmean+decdot*(tt[ii]-meant)/3.6D6, $
          astr[ii],tmpx,tmpy
        thisvar= (fwhm[ii]/sqrt8ln2)^2 ; check this!
        thisx= xx-tmpx
        thisy= yy-tmpy
        model= 1.0/(2.0*!DPI*thisvar) $
          *exp(-0.5*(thisx*thisx+thisy*thisy)/thisvar)
        image[*,*,ii]= flux*model
        foo= model/total(model)
        snr[ii]= flux*sqrt(total(foo*foo)*invvar)
    endfor

; add noise
    image= image+randomn(seed,nx,ny,nepoch)/sqrt(invvar)

; make fits files
    weight= fltarr(nx,ny)+invvar
    mkhdr, hdr,0
    sxaddpar, hdr,'NEPOCH',nepoch
    sxaddpar, hdr,'MEANT',meant
    sxaddpar, hdr,'FLUX',flux
    sxaddpar, hdr,'RAREF',raref
    sxaddpar, hdr,'DECREF',decref
    sxaddpar, hdr,'DRAMEAN',dramean
    sxaddpar, hdr,'DDECMEAN',ddecmean
    sxaddpar, hdr,'RADOT',radot ; mas/yr
    sxaddpar, hdr,'DECDOT',decdot ; mas/yr
    scales= [1,1,1]*2.0
    filename= prefix+'.fits'
    mwrfits, 0,filename,hdr,/create
    for ii=0,nepoch-1 do begin
        mkhdr, thishdr,weight
        putast, thishdr,astr[ii]
        sxaddpar, thishdr,'TT',tt[ii]
        sxaddpar, thishdr,'FWHM',fwhm[ii]
        sxaddpar, thishdr,'SNR',snr[ii]
        thisimage= reform(image[*,*,ii],nx,ny)
        mwrfits, thisimage,filename,thishdr
        mwrfits, weight,filename
    endfor
    spawn, 'gzip -v --best '+filename

; end loop over realizations
endfor
return
end
