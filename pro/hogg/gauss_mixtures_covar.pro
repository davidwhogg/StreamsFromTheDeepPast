;+
; NAME:
;   gauss_mixtures_covar
; PURPOSE:
;   iterate on gaussian mixtures (including data covariances)
; INPUTS:
;   ngauss - number of desired gaussians to fit
;   ydata - [ndim,ndata] observed velocities
;   ycovar - [ndim,ndim,ndata] observed velocities' errors
; OUTPUTS:
; BUGS:
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro gauss_mixtures_covar_step, data, covar,gauss_logamp,gauss_mean, $
                               gauss_covar, loglikedata=loglikedata, $
                               quiet=quiet, debug=debug

ngauss=n_elements(logamp)
ndim=n_elements(covar)/n_elements(data)
ndata=n_elements(data)/ndim
twopiterm=0.5*double(ndim)*alog(2.*!DPI)

loglikedata=0.D
logposterior=dblarr(ngauss,ndata)
for ii=0L, ndata-1L do begin
;   calculate log likelihood for each gaussian for this data point
    loglike=dblarr(ngauss)
    for jj=0L, ngauss-1L do begin
        tinv=special_invert(gauss_covar[*,*,jj]+covar[*,*,ii])
        delta=data[*,ii]-gauss_mean[*,jj]
        loglike[jj]=gauss_logamp[jj]+0.5*alog(determ(tinv,/double))- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    endfor

;   calculate posterior probability to be in each gaussian
    curr_loglikedata=logsum(loglike,/double)
    loglikedata=loglikedata+curr_loglikedata
    logposterior[*,ii]=loglike-curr_loglikedata
endfor

if(NOT keyword_set(quiet)) then $
  splog,'loglikedata= '+strtrim(string(loglikedata,format='(f40.9)'),2)

for jj=0L, ngauss-1L do begin

;   accumulate stuff for this gaussian
    curr_gauss_mean=dblarr(ndim)
    curr_gauss_covar1=dblarr(ndim,ndim)
    curr_gauss_covar2=dblarr(ndim,ndim)
    
    for ii=0L, ndata-1L do begin
        tinv=special_invert(gauss_covar[*,*,jj]+covar[*,*,ii])
        delta=ydata[*,ii]-gauss_mean[*,jj]
        posterior=exp(logposterior[jj,ii])
        curr_gauss_guess=gauss_mean[*,jj]+gauss_covar[*,*,jj]#tinv#delta
        curr_gauss_mean=curr_gauss_mean+posterior*curr_gauss_guess
        curr_gauss_covar1=curr_gauss_covar1+ $
          posterior*(curr_gauss_guess##curr_gauss_guess)
;       fix this if you fix the below 
        curr_gauss_covar2=curr_gauss_covar2+posterior*tinv

        if(n_elements(xcovarguess) gt 0) then begin
;           fix this if you fix the above 
            gauss_covarguess[*,*,ii,jj]=gauss_covar[*,*,jj]- $
              gauss_covar[*,*,jj]#tinv#transpose(gauss_covar[*,*,jj])
        endif
    endfor

;   now update parameters 
    logamp[jj]=logsum(logposterior[jj,*])
    curr_ndata=exp(logamp[jj])
    logamp[jj]=logamp[jj]-alog(double(ndata))
    gauss_mean[*,jj]=curr_gauss_mean/curr_ndata
    gauss_covar[*,*,jj]=gauss_covar[*,*,jj]- $
      gauss_covar[*,*,jj]#curr_gauss_covar2#transpose(gauss_covar[*,*,jj])/ $
      curr_ndata+curr_gauss_covar1/curr_ndata- $
      gauss_mean[*,jj]##gauss_mean[*,jj]

;   check for non-symmetry is NOT clean
    nonsym=max(abs(gauss_covar[*,*,jj]-transpose(gauss_covar[*,*,jj]))) $
      gt 1.D-10
    if(nonsym) then stop

    if(keyword_set(debug)) then $
      splog,'jj= '+string(jj)+' ; curr_ndata= '+string(curr_ndata)
endfor

end
;
pro gauss_mixtures_covar, ngauss, data, covar, $
                          gauss_amp, gauss_mean, gauss_covar, $
                          loglikedata=loglikedata, seed=seed, $
                          quiet=quiet, tol=tol, maxiter=maxiter

if(n_params() lt 3) then begin
    print,'Usage: gauss_mixtures_covar,ngauss,data,covar [,gauss_amp, gauss_mean, $'
    print,'          gauss_covar, loglikedata=, seed=, /quiet, tol=, maxiter=]'
endif

ndim=n_elements(covar)/n_elements(data)
ndata=n_elements(data)/ndim

if(NOT keyword_set(tol)) then tol=1.D-6
if(NOT keyword_set(maxiter)) then maxiter=1000000000L

; set initial conditions if necessary
if(n_elements(gauss_amp) ne ngauss or $
   n_elements(gauss_mean) eq 0 or $
   n_elements(gauss_covar)) then begin
    kmeans_covar, ngauss, data, covar, igauss, gauss_mean, gauss_covar, $
      seed=seed
    for i=0L, ngauss-1L do begin
        gauss_indx=where(igauss eq i, gauss_count)
        if(gauss_count eq 0) then $
          gauss_amp=1./double(ndata) $
        else $
          gauss_amp=double(gauss_count)/double(ndata)
    endif
    sumamp=total(gauss_amp,/double)
    gauss_amp=gauss_amp/sumamp
endif
  
diff=2.*tol
niter=0
while(niter lt maxiter and diff gt tol) do begin
    
    gauss_logamp=alog(gauss_amp)
    gauss_mixtures_covar_step, data, covar,  $
      gauss_logamp, gauss_mean, gauss_covar, loglikedata=loglikedata, $
      quiet=quiet, debug=debug
    gauss_amp=exp(gauss_logamp)
      
    if(niter gt 0) then $
      diff=(loglikedata-oldloglikedata)
    if(diff lt 0.) then $
        message,'FATAL ERROR: log likelihood decreased.'
    oldloglikedata=loglikedata
    niter=niter+1L
endwhile

if(n_elements(xguess) gt 0) then begin
    gauss_mixtures_covar_step, data, covar, $
      gauss_logamp, gauss_mean, gauss_covar, loglikedata=loglikedata, $
      quiet=quiet, debug=debug
endif

end
