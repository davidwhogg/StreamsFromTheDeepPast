;+
; NAME:
;  faintpm_sdss_verify
; PURPOSE:
;  Try to fit a moving point source to multi-epoch SDSS data to verify
;  a proposed proper motion.
; INPUTS:
;  flux         - first-guess flux (nMgy)
;  RA,Dec       - first-guess position of the source at the epoch (deg)
;  radot,decdot - first-guess time derivatives of position at the
;                 epoch (mas/yr)
;  epoch        - time of applicability of RA,Dec (MJD)
;  fluxerr,raerr,decerr,[etc] - first guesses at uncertainties
; OPTIONAL INPUTS:
;  run          - if set, array of SDSS runs to use
;  filename     - use this filename for data cutouts
;  priorvec     - information for a gaussian prior (see faintpm_verify
;                 code for an explanation)
; OUTPUTS:
;  RA,Dec       - measured positions at the epoch (deg)
;  radot,decdot - measured time derivatives at the epoch (mas/yr)
;  fluxerr,raerr,decerr,[etc] - uncertainties
; OPTIONAL OUTPUTS:
;  nepoch        - number of epochs
;  meant,deltat,sigmat - mean time (MJD), time span (years) and
;                        root-variance of times (years)
; BUGS:
;  - Runs only on SDSS z-band images.
;  - Relies on invvar output of smosaic_cacheimg(), which is (a)
;    unsupported and (b) doesn't work when returning data from cache.
;    For this reason, smosaic_cacheimg() is run with ncache=1 to
;    disable the cache.  This needs to be replaced with direct calls
;    to the data -- don't know how to do this.
;  - Makes cutout at zero proper motion to avoid getting large-pm
;    runaway solutions.  That's a hack.  Ought to have a capability
;    for taking a prior of some kind.
;  - Will fail for very large input raerr,decerr,radoterr,decdoterr
;    because of image-size issues.
; LICENSE:
;  Copyright 2007 David W. Hogg (NYU).  All rights reserved.
;-
pro faintpm_sdss_verify, flux,ra,dec,radot,decdot,epoch, $
                         fluxerr,raerr,decerr,radoterr,decdoterr, $
                         run=run,filename=filename,nepoch=nepoch, $
                         meant=meant,deltat=deltat,sigmat=sigmat, $
                         priorvec=priorvec

; set defaults
ngridhalf= 2 ; (make a 5-d cube in parameter space of side 2*ngrid+1)
if (NOT keyword_set(filename)) then begin
    name= hogg_iau_name(ra,dec)
    filename= strjoin(strsplit(name,' ',/extract),'_')+'.fits'
    filedir= strmid(filename,6,2)
    spawn, 'mkdir -p '+filedir
    filename= filedir+'/'+filename
endif
nepoch= 0
meant= 0.0
deltat= 0.0
sigmat= 0.0

; get sdss data
clobber= 0
if file_test(filename) then begin
    if (sxpar(headfits(filename),'NEPOCH') EQ 0) then clobber= 1
endif
if ((NOT file_test(filename)) OR clobber) then begin
    splog, 'making '+filename
    faintpm_sdss_cutout, ra,dec,0.0,0.0,epoch,filename
endif
splog, 'using cutouts in '+filename
nepoch= sxpar(headfits(filename),'NEPOCH')
nmin= 6
if (nepoch LT nmin) then begin
    splog, 'fewer than',nmin,' epochs!'
    splog, 'returning without change'
    return
endif

; get time statistics
tlist= fltarr(nepoch)
for ee=0,nepoch-1 do tlist[ee]= sxpar(headfits(filename,exten=2*ee+1),'MJD')
meant= total(tlist)/float(nepoch)
deltat= (max(tlist)-min(tlist))/365.25 ; yr
sigmat= sqrt(total((tlist-meant)^2)/float(nepoch))/365.25 ; yr

; run delta-chisq on a grid
ngrid= 2*ngridhalf+1
nsigma= double(ngridhalf)
gridvec= nsigma*(dindgen(ngrid)-ngridhalf)/double(ngridhalf)
ragrid= ra+raerr*gridvec
decgrid= dec+decerr*gridvec
radotgrid= radot+radoterr*gridvec
decdotgrid= decdot+decdoterr*gridvec
fluxgrid= flux+fluxerr*gridvec
dcs= faintpm_verify_grid(ragrid,decgrid,radotgrid,decdotgrid,fluxgrid, $
                         epoch,'MJD',filename,tunit=1D0/3.6525D2, $
                         priorvec=priorvec)

; marginalize and optimize for each parameter
tmpdcs= faintpm_marginalize(dcs,5)
three= (sort(tmpdcs))[0:2]
tmpxx= fluxgrid[three]
tmpyy= faintpm_parabola(tmpxx,tmpdcs[three],tmpxx,flux,ymin,fluxerr)
splog, 'flux:',flux,fluxerr

tmpdcs= faintpm_marginalize(dcs,1)
three= (sort(tmpdcs))[0:2]
tmpxx= ragrid[three]
tmpyy= faintpm_parabola(tmpxx,tmpdcs[three],tmpxx,ra,ymin,raerr)
splog, 'ra:',ra,raerr,' at mjd',epoch

tmpdcs= faintpm_marginalize(dcs,2)
three= (sort(tmpdcs))[0:2]
tmpxx= decgrid[three]
tmpyy= faintpm_parabola(tmpxx,tmpdcs[three],tmpxx,dec,ymin,decerr)
splog, 'dec:',dec,decerr,' at mjd',epoch

tmpdcs= faintpm_marginalize(dcs,3)
three= (sort(tmpdcs))[0:2]
tmpxx= radotgrid[three]
tmpyy= faintpm_parabola(tmpxx,tmpdcs[three],tmpxx,radot,ymin,radoterr)
splog, 'dra/dt:',radot,radoterr

tmpdcs= faintpm_marginalize(dcs,4)
three= (sort(tmpdcs))[0:2]
tmpxx= decdotgrid[three]
tmpyy= faintpm_parabola(tmpxx,tmpdcs[three],tmpxx,decdot,ymin,decdoterr)
splog, 'ddec/dt:',decdot,decdoterr

return
end
