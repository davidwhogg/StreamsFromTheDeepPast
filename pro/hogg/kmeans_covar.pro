;+
; NAME:
;   kmeans 
; PURPOSE:
;   K-means a distribution
; CALLING SEQUENCE:
;   kmeans, ngroup, data, covar, group [,group_mean=,group_covar=,seed=]
; INPUTS:
;   ngroup - number of groups
;   data - [M,N] data points 
;   covar - [M,M,N] data covariances 
; OPTIONAL INPUTS:
;   seed - seed for initial conditions
; OUTPUTS:
;   group - [N] group assignments for each data point
; OPTIONAL OUTPUTS:
;   group_mean - [M,ngroup] mean of each group
;   group_covar - [M,M,ngroup] covar of each group
; BUGS:
; REVISION HISTORY:
;   2003-03-03  written - Blanton Roweis and Hogg
;-
pro kmeans_covar, ngroup, data, covar, group, group_mean=group_mean, $
                  group_covar=group_covar, seed=seed

epsilon=1.e-2
ndim=n_elements(covar)/n_elements(data)
ndata=n_elements(data)/ndim

; set initial conditions if necessary (means at ngroup random points)
if(n_elements(group_mean) ne ngroup*ndim) then begin
    rearrange=shuffle_indx(ndata,num_sub=ngroup,seed=seed) 
    group_mean=data[*,rearrange]
endif

group=lonarr(ndata)
converged=0

while(converged eq 0) do begin
;   assigned groups
    total_error=0.D
    converged=1
    minerror=dblarr(ndata)
    for ii=0L, ndata-1L do begin
        error_group=dblarr(ngroup)
        invcovar=invert(covar[*,*,ii])
        for jj=0L, ngroup-1L do begin
            delta=data[*,ii]-group_mean[*,jj]
            error_group[jj]=transpose(delta)#invcovar#delta
        endfor
        minerror[ii]=min(error_group,jjminerror)
        total_error=total_error+minerror[ii]
        if(group[ii] ne jjminerror) then begin
            converged=0
            group[ii]=jjminerror
        endif
    endfor
    splog,'total_error= '+string(total_error)

;   update group centers
    if(converged eq 0) then begin
        for jj=0L, ngroup-1L do begin
            indx_group=where(group eq jj, indx_count)
            if(indx_count gt 0) then begin
                denom=dblarr(ndim,ndim)
                numer=dblarr(ndim,ndim)
                for kk=0L, n_elements(indx_group)-1L do begin
                    ii=indx_group[kk]
                    invcovar=invert(covar[*,*,ii])
                    denom=denom+invcovar
                    numer=numer+invcovar#data[*,ii]
                endfor
                group_mean[*,jj]=invert(denom)#numer
            endif else begin
                maxminerror=max(minerror,iiminerror)
                group_mean[*,jj]=data[*,ii]
                minerror[iiminerror]=0.D
            endelse
        endfor
    endif
endwhile

if(arg_present(group_covar)) then begin
    group_covar=dblarr(ndim,ndim,ngroup)
    for jj=0L, ngroup-1L do begin
        indx_group=where(group eq jj, indx_count)
        group_covar[*,*,jj]=transpose(data[*,indx_group])## $
          data[*,indx_group]/double(indx_count)
        if(indx_count le ndim) then $
          group_covar[*,*,jj]=group_covar[*,*,jj]+epsilon*identity(ndim)
    endfor
endif

end
