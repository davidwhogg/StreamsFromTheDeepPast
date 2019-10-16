;+
; BUGS:
;  - No proper header.
;  - Hard-coded to z band.
;  - Epoch hacked.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_lmcc_test
band= 4
ddmy= 365.25*3.6D6 ; conversion factor
prefix= 'faintpm_lmcc'

; read structure and make mine
lmcc= mrdfits('../../data/faintpm/run1_LMCC.fits*',1)
good= where((lmcc.mean_psfmag[band] GT 18.0) AND $
            (lmcc.mean_psfmag[band] LT 20.0) AND $
            (lmcc.n_good_epochs[band] GE 10) AND $
            ((lmcc.ra_pm*ddmy) LT 400.0) AND $
            ((lmcc.dec_pm*ddmy) LT 400.0),nlmcc)
lmcc= lmcc[good]
hogg= faintpm_struct(nlmcc)

; loop over good objects
for ii=0L,nlmcc-1L do begin
    lmccii= lmcc[ii]
    splog, ii,'; z =',lmccii.mean_psfmag[band],' mag;', $
      lmccii.n_good_epochs[band],' epochs'
    flux= 1D1^(0.4*(22.5-lmccii.mean_psfmag[band])) ; nMgy
    ra= lmccii.ra_mean          ; deg
    dec= lmccii.dec_mean        ; deg
    radot= lmccii.ra_pm*ddmy    ; mas/yr ?
    decdot= lmccii.dec_pm*ddmy ; mas/yr ?
    fluxerr= lmccii.rms_psfmag[band]/sqrt(lmccii.n_good_epochs[band])
    raerr= lmccii.ra_mean_err
    decerr= lmccii.dec_mean_err
    radoterr= lmccii.ra_pm_err*ddmy
    decdoterr= lmccii.dec_pm_err*ddmy
    hogg[ii].epoch= 0.5*(lmccii.ra_t0+lmccii.dec_t0) ; mjd ; HACK!
    faintpm_sdss_verify_iterate, flux,ra,dec,radot,decdot,hogg[ii].epoch, $
      fluxerr,raerr,decerr,radoterr,decdoterr,nepoch=nepoch, $
      meant=meant,deltat=deltat,sigmat=sigmat
    hogg[ii].nepoch= nepoch
    hogg[ii].meant= meant
    hogg[ii].deltat= deltat
    hogg[ii].sigmat= sigmat
    hogg[ii].flux= flux
    hogg[ii].ra= ra
    hogg[ii].dec= dec
    hogg[ii].radot= radot
    hogg[ii].decdot= decdot
    hogg[ii].fluxerr= fluxerr
    hogg[ii].raerr= raerr
    hogg[ii].decerr= decerr
    hogg[ii].radoterr= radoterr
    hogg[ii].decdoterr= decdoterr

    if (((ii MOD 10) EQ 0) OR $
        (ii EQ (nlmcc-1L))) then begin
        filename= prefix+'.fits'
        mwrfits, lmcc,filename,/create
        mwrfits, hogg,filename
    endif
endfor
return
end
