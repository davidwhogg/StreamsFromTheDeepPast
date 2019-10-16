;+
; NAME:
;   gauss_mixtures
; PURPOSE:
;   iterate on gaussian mixtures 
; CALLING SEQUENCE:
;   gauss_mixtures, ngauss, data, amp, mean, covar [, weight=, 
;     loglikedata=, seed=, /quiet, tol=, maxiter=, peg_dims= ]
; INPUTS:
;   ngauss      - number of desired gaussians to fit
;   data        - [ndim,ndata] observations
; OPTIONAL INPUTS:
;   weight      - [ndata] weights for data
;   peg_dims    - { read the source, Luke }
;   maxiter     - { read the source, Luke }
;   tol         - { read the source, Luke }
; KEYWORDS:
;   quiet       - minimize output
;   debug       - maximize output
; OUTPUTS:
;   gauss_amp   - [ngauss]
;   gauss_mean  - [ndim,ngauss]
;   gauss_covar - [ndim,ndim,ngauss]
; OPTIONAL OUTPUTS:
;   loglikedata
; BUGS:
;   Doesn't work with non-trivial weights -- we get singular and
;     negative variances!!
;   IDL sucks, so covariance matrices go non-symmetric; this code just
;     "re-symmetrizes" them, which is VERY BAD.
;   Weights are not passed to kmeans.
;   Hogg changed the restriction that "peg_dims" doesn't allow
;     non-zero covariances of the pegged dimensions with anything
;     else to no covariances with anything unpegged; ie, no
;     peg-unpeg covariances.
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;   2003-03-30  added weights - Hogg
;-
pro gauss_mixtures_step, data, gauss_logamp, gauss_mean, gauss_covar, $
                         weight=weight, loglikedata=loglikedata, $
                         quiet=quiet, debug=debug, peg_dims=peg_dims

ngauss=n_elements(gauss_logamp)
ndim=n_elements(gauss_mean)/ngauss
ndata=n_elements(data)/ndim

data=reform(data, ndim, ndata)
gauss_logamp=reform(gauss_logamp, ngauss)
gauss_mean=reform(gauss_mean, ndim, ngauss)
gauss_covar=reform(gauss_covar, ndim, ndim, ngauss)

twopiterm=0.5*double(ndim)*alog(2.*!DPI)

; deal with pegging dimensions
notpegged_mean=lindgen(ndim)
notpegged_covar=lindgen(ndim*ndim)
notpegged_count=ndim
if(n_elements(peg_dims) gt 0) then begin
    notpegged=1+lonarr(ndim)
    notpegged[peg_dims]=0
    notpegged_mean=where(notpegged gt 0, notpegged_count)
    if(notpegged_count gt 0) then begin
        notpegged_covar=reform(notpegged_covar,ndim,ndim)
        notpegged_covar[*,peg_dims]=-1
        notpegged_covar[peg_dims,*]=-1
    endif
endif
indx=where(notpegged_covar ge 0)
notpegged_covar_1=notpegged_covar[indx]/ndim
notpegged_covar_2=notpegged_covar[indx] mod ndim

; invert matrices
tinvlist= gauss_covar
dtinvlist= dblarr(ngauss)
for jj=0L, ngauss-1L do begin
    tinvlist[*,*,jj]=special_invert(gauss_covar[*,*,jj])
    if(ndim gt 1) then $
      dtinvlist[jj]= determ(tinvlist[*,*,jj],/double) $
    else $
      dtinvlist[jj]= tinvlist[0,0,jj]
endfor

loglikedata=0.D
logposterior=dblarr(ngauss,ndata)

for ii=0L, ndata-1L do begin
;   calculate log likelihood for each gaussian for this data point
    loglike=dblarr(ngauss)
    for jj=0L, ngauss-1L do begin
        tinv= tinvlist[*,*,jj]
        delta=data[*,ii]-gauss_mean[*,jj]
        loglike[jj]=gauss_logamp[jj]+0.5*alog(dtinvlist[jj])- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    endfor

;   calculate posterior probability to be in each gaussian
    curr_loglikedata=logsum(loglike)
    
    loglikedata=loglikedata+curr_loglikedata*weight[ii]
    logposterior[*,ii]=loglike-curr_loglikedata
endfor

if(NOT keyword_set(quiet)) then $
  splog,'loglikedata= '+strtrim(string(loglikedata,format='(e20.9)'),2)
if(NOT keyword_set(quiet)) then $
  splog, gauss_logamp
if(NOT keyword_set(quiet)) then $
  splog, gauss_mean
if(NOT keyword_set(quiet)) then $
  splog, gauss_covar
if(NOT keyword_set(quiet)) then $
  splog, dtinvlist

for jj=0L, ngauss-1L do begin

;   accumulate stuff for this gaussian
    curr_gauss_mean=dblarr(ndim)
    curr_gauss_covar1=dblarr(ndim,ndim)
    curr_gauss_covar2=dblarr(ndim,ndim)
    tinv=tinvlist[*,*,jj]
 
    for ii=0L, ndata-1L do begin
        delta=data[*,ii]-gauss_mean[*,jj]
        posterior=exp(logposterior[jj,ii])
        curr_gauss_guess=gauss_mean[*,jj]+gauss_covar[*,*,jj]#tinv#delta
        curr_gauss_mean=curr_gauss_mean $
          +(posterior*curr_gauss_guess)*weight[ii]
        curr_gauss_covar1=curr_gauss_covar1+ $
          posterior*(curr_gauss_guess##curr_gauss_guess)*weight[ii]
;       fix this if you fix the below 
        curr_gauss_covar2=curr_gauss_covar2+posterior*tinv*weight[ii]
    endfor

;   now update parameters 
    gauss_logamp[jj]=logsum(logposterior[jj,*]+alog(weight))
    curr_ndata=exp(gauss_logamp[jj])
    gauss_logamp[jj]=gauss_logamp[jj]-alog(total(weight,/double))
    if(n_elements(peg_dims) eq 0) then begin
        gauss_mean[*,jj]=curr_gauss_mean/curr_ndata
        upcovar= -gauss_covar[*,*,jj]#(curr_gauss_covar2/curr_ndata)# $
          transpose(gauss_covar[*,*,jj])+(curr_gauss_covar1/curr_ndata)- $
          gauss_mean[*,jj]##gauss_mean[*,jj]
        gauss_covar[*,*,jj]= gauss_covar[*,*,jj]+upcovar[*,*]
    endif else if(notpegged_count gt 0) then begin
        gauss_mean[notpegged_mean,jj]=curr_gauss_mean[notpegged_mean]/ $
          curr_ndata
        upcovar= -gauss_covar[*,*,jj]#(curr_gauss_covar2/curr_ndata)# $
          transpose(gauss_covar[*,*,jj]) $
          +(curr_gauss_covar1/curr_ndata)- $
          gauss_mean[*,jj]##gauss_mean[*,jj]
        gauss_covar[notpegged_mean,notpegged_mean,jj]= $
          gauss_covar[notpegged_mean,notpegged_mean,jj]+ $
          upcovar[notpegged_covar_1,notpegged_covar_2]
    endif

; re-symmetrize
; TOTAL HACK
    gauss_covar[*,*,jj]= $
      5D-1*(gauss_covar[*,*,jj]+transpose(gauss_covar[*,*,jj]))

; talk
    if(keyword_set(debug)) then $
      splog,'jj= '+string(jj)+' ; curr_ndata= '+string(curr_ndata)
endfor

end
;
pro gauss_mixtures, ngauss, data, $
                    gauss_amp, gauss_mean, gauss_covar, $
                    weight=weight, $
                    loglikedata=loglikedata, seed=seed, $
                    quiet=quiet, tol=tol, maxiter=maxiter, $
                    peg_dims=peg_dims

if(n_params() lt 2) then begin
    print,'Usage: gauss_mixtures,ngauss,data [,gauss_amp, gauss_mean, $'
    print,'          gauss_covar, weight=, loglikedata=, seed=, /quiet, tol=, maxiter=]'
endif

if(n_elements(size(data,/dimensions)) eq 1) then $
  ndim=1 $
else $
  ndim=(size(data,/dimensions))[0]
ndata=n_elements(data)/ndim
splog, 'fitting ',ndata,' data points of ',ndim,' dimensions with ',ngauss,' gaussians'
flush, -1

if(NOT keyword_set(tol)) then tol=1.D-6
if(NOT keyword_set(maxiter)) then maxiter=1000000000L
if(NOT keyword_set(weight)) then weight= replicate(1D0,ndata)

; set initial conditions if necessary
if(n_elements(gauss_amp) ne ngauss or $
   n_elements(gauss_mean) eq 0 or $
   n_elements(gauss_covar) eq 0) then begin
    kmeans_streams, ngauss, data, igauss, group_mean=gauss_mean, $
      group_covar=gauss_covar, seed=seed, quiet=quiet
    gauss_amp=dblarr(ngauss)
    for i=0L, ngauss-1L do begin
        gauss_indx=where(igauss eq i, gauss_count)
        if(gauss_count eq 0) then $
          gauss_amp[i]=1./double(ndata) $
        else $
          gauss_amp[i]=double(gauss_count)/double(ndata)
    endfor
    sumamp=total(gauss_amp,/double)
    gauss_amp=gauss_amp/sumamp
endif

; if the dimension is pegged, we don't let it have covariances with
; any unpegged dimension
npeg= n_elements(peg_dims)
if (npeg GT 0) then begin
    for kk=0L,ndim-1 do begin
        blah= where(peg_dims EQ kk,nkk)
        if nkk EQ 0 then begin
            for jj=0,ngauss-1 do begin
                gauss_covar[kk,peg_dims,jj]= 0D0
                gauss_covar[peg_dims,kk,jj]= 0D0
            endfor
        endif
    endfor
endif
  
diff=2.*tol
niter=0
while(niter lt maxiter and diff gt tol) do begin
    
    gauss_logamp=alog(gauss_amp)
    gauss_mixtures_step, data, $
      gauss_logamp, gauss_mean, gauss_covar, $
      weight=weight, loglikedata=loglikedata, $
      quiet=quiet, debug=debug, peg_dims=peg_dims
    gauss_amp=exp(gauss_logamp)
      
    if(niter gt 0) then $
      diff=(loglikedata-oldloglikedata)
    if(diff lt 0.) then begin
        splog,'NOT FATAL ERROR: log likelihood decreased.'
        splog,'NOT FATAL ERROR: Either Blanton or Roweis or Hogg or RSI suck'
    endif
    oldloglikedata=loglikedata
    niter=niter+1L
endwhile

end
