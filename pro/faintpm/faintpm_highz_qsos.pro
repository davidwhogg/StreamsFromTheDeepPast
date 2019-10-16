;+
; BUGS:
;  - No proper header.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_highz_qsos,bd=bd,clobber=clobber
epoch= 53000.0

if keyword_set(bd) then begin
    readcol, '../../data/faintpm/Jiang_LTdata.dat',foo,bar,format='A,A'
    sdss= where((foo EQ 'SDSS'))
    bdlist= foo[sdss]+' '+bar[sdss]
    outfilename= 'faintpm_bd.fits'
endif else begin
    bdlist= ['SDSS J000552.34-000655.8', $
             'SDSS J020332.39+001229.3', $
             'SDSS J030331.40-001912.9', $
             'SDSS J035349.72+010404.4', $
             'SDSS J205406.49-000514.8', $
             'SDSS J231546.57-002358.1']
    outfilename= 'faintpm_highz.fits'
endelse
nbd= n_elements(bdlist)
hogg= faintpm_struct(nbd)

for qq=0,nbd-1 do begin
    flux= 3.5                   ; nMgy
    hogg_name2radec, bdlist[qq],ra,dec
    radot= 0.0
    decdot= 0.0
    fluxerr= 0.2                ; nMgy
    raerr= 1D-5                 ; deg
    decerr= raerr
    radoterr= 30.0              ; mas/yr
    decdoterr= radoterr
    filename= strjoin(strsplit(bdlist[qq],' ',/extract),'_')+'.fits'
    if keyword_set(clobber) then spawn, '\rm -vf '+filename+'*'
    hogg[qq].epoch= epoch
    faintpm_sdss_verify_iterate, flux,ra,dec,radot,decdot,hogg[qq].epoch, $
      fluxerr,raerr,decerr,radoterr,decdoterr,filename=filename, $
      nepoch=nepoch,meant=meant,deltat=deltat,sigmat=sigmat
    hogg[qq].nepoch= nepoch
    hogg[qq].meant= meant
    hogg[qq].deltat= deltat
    hogg[qq].sigmat= sigmat
    hogg[qq].flux= flux
    hogg[qq].ra= ra
    hogg[qq].dec= dec
    hogg[qq].radot= radot
    hogg[qq].decdot= decdot
    hogg[qq].fluxerr= fluxerr
    hogg[qq].raerr= raerr
    hogg[qq].decerr= decerr
    hogg[qq].radoterr= radoterr
    hogg[qq].decdoterr= decdoterr
endfor

mwrfits, hogg,outfilename,/create
return
end
