;+
; NAME:
;  faintpm_sdss_verify_iterate
; PURPOSE:
;  Sensibly iterate the faintpm_sdss_verify code.
; INPUTS:
;  flux         - first-guess flux (nMgy)
;  RA,Dec       - first-guess position of the source at the epoch (deg)
;  radot,decdot - first-guess time derivatives of position at the
;                 epoch (mas/yr)
;  epoch        - time of applicability of RA,Dec (MJD)
;  fluxerr,raerr,decerr,[etc] - first guesses at uncertainties
; OPTIONAL INPUTS:
;  [extra]     - all extras passed to faintpm_sdss_verify
; OUTPUTS:
;  RA,Dec       - measured positions at the epoch (deg)
;  radot,decdot - measured time derivatives at the epoch (mas/yr)
;  fluxerr,raerr,decerr,[etc] - uncertainties
; OPTIONAL OUTPUTS:
;  [extra]     - all extras passed from faintpm_sdss_verify
; BUGS:
;  - See all bugs in faintpm_sdss_verify.  There are MANY.
;  - Doesn't re-build data cutout set even if RA,Dec change substantially.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_sdss_verify_iterate, flux,ra,dec,radot,decdot,epoch, $
                                 fluxerr,raerr,decerr,radoterr,decdoterr, $
                                 run=run,filename=filename,nepoch=nepoch, $
                                 meant=meant,deltat=deltat,sigmat=sigmat

influx= flux
inra= ra
indec= dec
inradot= radot
indecdot= decdot
influxerr= fluxerr
inraerr= raerr
indecerr= decerr
inradoterr= radoterr
indecdoterr= decdoterr

pospriorsigma= 0.2/3600.0       ; deg
dotpriorsigma= 200.0            ; mas/yr
fluxprior= 10.0
fluxpriorsigma= 10.0
priorvec=[inra,pospriorsigma,indec,pospriorsigma, $
          inradot,dotpriorsigma,indecdot,dotpriorsigma, $
          influx,fluxpriorsigma]

for iter=0,9 do begin
    if (finite(flux) NE 1) then flux= influx
    if (finite(ra) NE 1) then ra= inra
    while (ra LT 0D0) do ra= ra+3.6D2
    while (ra GT 3.6D2) do ra= ra-3.6D2
    if (finite(dec) NE 1) then dec= indec
    if (finite(radot) NE 1) then radot= inradot
    if (finite(decdot) NE 1) then decdot= indecdot
    if ((fluxerr LE 0.0) OR (flux LE 0.0)) then begin
        flux= influx
        fluxerr= influxerr*1.1^(iter)
        ra= inra
        dec= indec
        radot= inradot
        decdot= indecdot
        raerr= inraerr
        decerr= indecerr
        radoterr= inradoterr
        decdoterr= indecdoterr
    endif
    if ((raerr LE 0.0) or (raerr GT 3d-4)) then begin
        raerr= inraerr*1.1^(iter)
        radot= inradot
        decdot= indecdot
        radoterr= inradoterr
        decdoterr= indecdoterr
    endif
    if ((decerr LE 0.0) or (decerr GT 3d-4)) then begin
        decerr= indecerr*1.1^(iter)
        radot= inradot
        decdot= indecdot
        radoterr= inradoterr
        decdoterr= indecdoterr
    endif
    if ((radoterr LE 0.0) or (radoterr GT 1000.0)) then begin
        radoterr= inradoterr*1.1^(iter)
        radot= inradot
        decdot= indecdot
    endif
    if ((decdoterr LE 0.0) or (decdoterr GT 1000.0)) then begin
        decdoterr= indecdoterr*1.1^(iter)
        radot= inradot
        decdot= indecdot
    endif
    faintpm_sdss_verify, flux,ra,dec,radot,decdot,epoch, $
      fluxerr,raerr,decerr,radoterr,decdoterr, $
      run=run,filename=filename,nepoch=nepoch, $
      meant=meant,deltat=deltat,sigmat=sigmat, $
      priorvec=priorvec
endfor
return
end
