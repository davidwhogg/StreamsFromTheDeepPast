;+
; NAME:
;   uber
; PURPOSE:
;   Execute 1-d, 1-gaussian uber algorithm
; INPUTS:
;   weight    - vector of data-point weights
;   point     - vector of data-point values
;   covar     - vector of data-point measurement variances
;   nvec      - number of eigenvectors to fit
; OPTIONAL INPUTS:
;   nboot     - number of bootstrap resamples to do; default zero
;   seed      - random seed for bootstrap resamples
; OUTPUTS:
;   mean      - mean of the external distribution
;   var       - variance of the external distribution
; OPTIONAL OUTPUTS:
;   bootmean  - vector of bootstrap means
;   bootvar   - vector of bootstrap variances
; COMMENTS:
;   Works better with a prior variance which is way larger than the
;     true value than a prior variance which is way smaller.
; DEPENDENCIES:
;   idlutils
; BUGS:
;   doesn't initialize to PCA solution
; REVISION HISTORY:
;   2002-04-09  (re)written - Hogg
;-
function hogg_determ,a
  if n_elements(a) EQ 1 then return, double(a[0]) else return, determ(a,/double)
end

function uber_log_like,npoint,ndimen,nvec,weight,point,covar,mean,vec

; build external variance matrix
  extvar= transpose(vec)##vec

; start with orthogonalization constraint
  sum= 0D0
;  for ii=0L,nvec-1L do begin
;    for jj=0L,ii-1L do sum= sum-alog((vec[*,ii]##transpose(vec[*,jj]))[0])
;  endfor

; loop over points
  for ii= 0,npoint-1 do begin
    thispoint= point[*,ii]-mean
    bigvar= covar[*,*,ii]+extvar
    invvar= invert(bigvar,/double)
    sum= sum+weight[ii]*5d-1*(alog(hogg_determ(invvar)) $
      -thispoint##invvar##transpose(thispoint))
  endfor

  return, sum
end

function uber_log_like_wrapper,pars,npoint=npoint,ndimen=ndimen,nvec=nvec,weight=weight,point=point,covar=covar
  mean= double(pars[0:ndimen-1L])
  vec= dblarr(ndimen,nvec)
  for ii= 0L,nvec-1L do vec[*,ii]= pars[(ii+1L)*ndimen:(ii+2L)*ndimen-1L]
  return, uber_log_like(npoint,ndimen,nvec,weight,point,covar,mean,vec)
end

pro uber,weight,point,covar,nvec,mean,vec, $
         nboot=nboot,seed=seed,bootmean=bootmean,bootvec=bootvec

; re-form and re-cast inputs into the common block
  npoint= n_elements(weight)
  ndimen= n_elements(point)/npoint
  if (nvec GT npoint) OR (nvec GT ndimen) then begin
    splog, "ERROR: nvec not compatible with npoint or ndimen"
    return
  endif
  weight= reform(double([weight]),npoint)
  point= reform(double([point]),ndimen,npoint)
  covar= reform(double([covar]),ndimen,ndimen,npoint)

; initialize for iteration
  mean= total(point,2)/double(npoint)
  par= [mean]
  for ii=0L,nvec-1L do par= [par,mean^ii]

; maximize
  par= tnmin('uber_log_like_wrapper',par,epsabs=1D-6, $
    /maximize,/autoderivative,/quiet,status=status, $
    functargs={npoint:npoint,ndimen:ndimen,nvec:nvec,weight:weight, $
               point:point,covar:covar})
  splog, 'tnmin completed with status',status
  mean= double(par[0:ndimen-1])
  vec= dblarr(ndimen,nvec)
  for ii=0L,nvec-1L do vec[*,ii]= par[(ii+1L)*ndimen:(ii+2L)*ndimen-1L]

; recursively bootstrap, if asked to do so
  if keyword_set(nboot) then if nboot GT 0 then begin
    bootmean= dblarr(ndimen,nboot)
    bootvec= dblarr(ndimen,nvec,nboot)
    for ii=0,nboot-1 do begin
      index= long(randomu(seed,npoint)*double(npoint))
      uber, weight[index],point[index],covar[index],nvec,tmean,tvec
      bootmean[*,ii]= tmean
      bootvec[*,*,ii]= tvec
    endfor
  endif
  return
end
