;+
; NAME:
;   projected_gauss_analyze
; PURPOSE:
;   analyze results of projected_gauss_mixtures (likelihood and predictions)
; INPUTS:
;   ydata - [ndimy,ndata] observed velocities
;   ycovar - [ndimy,ndimy,ndata] observed velocities' errors
;   projection - [ndimx, ndimy, ndata] non-square matrices
;                implementing the projection
; OUTPUTS:
; BUGS:
;   this should be linked to projected_gauss_step
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro projected_gauss_analyze, ydata, ycovar, projection, $
                             amp, xmean, xcovar, $
                             avgloglikedata=avgloglikedata, $
                             quiet=quiet

ngauss=n_elements(amp)
ndimy=n_elements(ycovar)/n_elements(ydata)
ndimx=n_elements(xcovar)/n_elements(xmean)
ndata=n_elements(ydata)/ndimy
twopiterm=0.5*double(ndimy)*alog(2.*!DPI)
logamp=alog(amp)

loglikedata=0.D
logposterior=dblarr(ngauss,ndata)
for ii=0L, ndata-1L do begin
;   calculate log likelihood for each gaussian for this data point
    loglike=dblarr(ngauss)
    for jj=0L, ngauss-1L do begin
;       may want to hack 2D invert case to be faster
        tinv=invert(transpose(projection[*,*,ii])#xcovar[*,*,jj]# $
                    projection[*,*,ii]+ycovar[*,*,ii])
        delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
        loglike[jj]=logamp[jj]+0.5*alog(determ(tinv))- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    endfor

;   calculate posterior probability to be in each gaussian
    curr_loglikedata=logsum(loglike,/double)
    loglikedata=loglikedata+curr_loglikedata
    logposterior[*,ii]=loglike-curr_loglikedata
endfor

avgloglikedata=loglikedata/double(ndata)
if(NOT keyword_set(quiet)) then $
  splog,'avgloglikedata= '+strtrim(string(avgloglikedata,format='(f40.9)'),2)

end
