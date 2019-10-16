;+
; BUGS:
;  - No proper header.
;  - Hard-coded to z band.
;  - Epoch hacked.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_coadd_catalog, clobber=clobber
prefix= 'fermi_coadd_catalog'
outfilename= prefix+'.fits'
savefilename= prefix+'.sav'
band= 4

; read or make structures
if (file_test(outfilename) and (not keyword_set(clobber))) then begin
    cc= mrdfits(outfilename,1)
    ncc= n_elements(cc)
    hogg= mrdfits(outfilename,2)
endif else begin
    cc= mrdfits('../../data/faintpm/FermiCoadd_allzmag_izgt2_rzgt2_all.fits.gz',1)
    good= where(((cc.status AND 8192) EQ 8192) AND $
                (cc.psfcounts[band] LT 21.0),ncc)
    cc= cc[good]
    hogg= faintpm_struct(ncc)
endelse
splog, 'found',ncc,' good sources'

; loop over good objects
for ii=0L,ncc-1L do begin
    if (hogg[ii].nepoch EQ 0) then begin
        splog, 'working on',ii
        ccii= cc[ii]
        splog, ii,'; z =',ccii.psfcounts[band],' mag'
        flux= 1D1^(0.4*(22.5-ccii.psfcounts[band])) ; nMgy
        ra= ccii.ra             ; deg
        dec= ccii.dec           ; deg
        radot= 0.0              ; mas/yr
        decdot= 0.0             ; mas/yr
        fluxerr= 0.1*flux
        raerr= 1.0e-5           ; deg
        decerr= 1.0e-5          ; deg
        radoterr= 200.0         ; mas/yr
        decdoterr= 200.0        ; mas/yr
        hogg[ii].epoch= 53000   ; mjd ; HACK!
        faintpm_sdss_verify_iterate, flux,ra,dec,radot,decdot,hogg[ii].epoch, $
          fluxerr,raerr,decerr,radoterr,decdoterr,filename=thisfilename, $
          nepoch=nepoch,meant=meant,deltat=deltat,sigmat=sigmat
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
        hogg[ii].filename= thisfilename
        thisfilename= 0
    endif else begin
        splog, 'skipping',ii
    endelse

    if (((ii MOD 100) EQ 0) OR $
        (ii EQ (ncc-1L))) then begin
        save, filename=savefilename
        splog, 'writing FITS file dont kill me!'
        mwrfits, cc,outfilename,/create
        mwrfits, hogg,outfilename
    endif
endfor
return
end
