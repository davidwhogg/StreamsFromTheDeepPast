;+
; NAME:
;   projected_gauss_mixtures_c
; PURPOSE:
;   iterate on projected gaussian mixtures using the C program
; INPUTS:
;   ngauss     - number of desired gaussians to fit
;   ydata      - [ndimy,ndata] observed velocities
;   ycovar     - [ndimy,ndimy,ndata] observed velocities' errors
;   projection - [ndimx, ndimy, ndata] non-square matrices
;                implementing the projection
;   amp        - [ngauss] list of relative amplitudes
;   xmean      - [ndimx,ngauss] list of initial guesses for the means
;   xcovar     - [ndimx,ndimx,ngauss] list of initial guesses for the
;                 covar
;
; OPTIONAL INPUTS:
;   fixamp     - [ngauss] list of integers: 0 to update amp, 1 not to
;   fixmean    - [ngauss] list of integers: 0 to update xmean, 1 not to
;   fixcovar   - [ngauss] list of integers: 0 to update xcovar, 1 not
;                to
;   tol        - tolerance (convergence iff difference in avgloglike <
;                tol)
;   maxiter    - maximum number of iterations
;   w          - constrain the any variance to be larger than this
;                value
;   logfile    - append messages to this logfile
;
; KEYWORDS:
;   likeonly   - only compute and return the avgloglikedata; do not
;                update
;   quiet      - don't print messages
;
; OUTPUTS:
;   avgloglikedata - Average log-likelihood of the data after
;                   convergence
;   +updata xmean, xcovar and amp...
;
; REVISION HISTORY:
;   2003-02-18  first implementation of projected_gauss_mixtures
;               written in IDL - Blanton and Roweis
;   2008-06-10  patched up + docs - Bovy
;   2008-07-17  added 'w' input - Bovy
;   2008-09-15  rewritten the algorithm in C, written this wrapper for
;               it - Bovy
;-
pro projected_gauss_mixtures_c, ngauss, ydata, ycovar, projection, $
                                amp, xmean, xcovar, $
                                fixamp=fixamp,fixmean=fixmean, $
                                fixcovar=fixcovar, $
                                avgloglikedata=avgloglikedata, $
                                tol=tol, maxiter=maxiter, $
                                likeonly=likeonly, w=w, $
                                quiet=quiet, logfile=logfile, $
                                splitnmerge=splitnmerge

;ON_ERROR, 2

;;Perform some massive checking of the inputs here, return if (some
;;critical) parameters aren't set
;;We don't want fixamp, ... to be changed by the C-code
if(NOT keyword_set(fixamp)) then fixamp_temp= bytarr(ngauss) else fixamp_temp = byte(fixamp)
if(NOT keyword_set(fixmean)) then fixmean_temp= bytarr(ngauss) else fixmean_temp = byte(fixmean)
if(NOT keyword_set(fixcovar)) then fixcovar_temp= bytarr(ngauss) else fixcovar_temp = byte(fixcovar)
if(NOT keyword_set(tol)) then tol=1.D-6
if(NOT keyword_set(maxiter)) then maxiter=1000000000L
if ~keyword_set(splitnmerge) then splitnmerge = 1L else splitnmerge = long(splitnmerge)
;;amplitudes have to add up to 1
if not (total(amp) EQ 1D) then amp = amp/total(amp)

;;For now this assumes that the dimensions of all datapoints are equal
ndimy=n_elements(ycovar)/n_elements(ydata)
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(projection)/ndata/ndimy

avgloglikedata=0.

;likeonly
If ~keyword_set(likeonly) THEN likeonly = 0B
likeonly = byte(likeonly)

;logfile
IF ~keyword_set(logfile) THEN clog = '' ELSE clog = logfile + '_c.log'
n_clog = strlen(clog)



clog = byte(clog)
value= bytarr(21)
value[3]=1
value[4]=1
value[8]=1
value[9]=1
value[14]=1
value[15]=1
value[16]=1
value[17]=1
value[19]=1
value[20]=1


IF ~keyword_set(quiet) THEN splog, "Loading C-implementation to run the projected_gauss_mixtures algorithm..."


result = CALL_EXTERNAL("../EM.so","EM_IDL", $
                       ydata, ycovar, projection, ndata, ndimy, $
                       amp, xmean, xcovar, ndimx, ngauss, $
                       fixamp_temp, fixmean_temp, fixcovar_temp, $
                       avgloglikedata, tol, maxiter, likeonly, $
                       w,clog,n_clog,splitnmerge, $
                       /CDECL,/AUTO_GLUE,VALUE=value)

IF ~keyword_set(quiet) THEN splog, "Successfully ran the C-code"


stop

end
