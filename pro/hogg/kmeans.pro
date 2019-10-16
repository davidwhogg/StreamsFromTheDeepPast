;+
; NAME:
;   kmeans 
; PURPOSE:
;   K-means a distribution
; CALLING SEQUENCE:
;   kmeans, ngroup, data, group [, group_mean=, seed= ]
; INPUTS:
;   ngroup - number of groups
;   data - [M,N] data points 
; OPTIONAL INPUTS:
;   seed - seed for initial conditions
; KEYWORDS:
;   quiet - shut up
; OUTPUTS:
;   group - [N] group assignments for each data point
; OPTIONAL OUTPUTS:
;   group_mean - [M,ngroup] mean of each group
; BUGS:
; REVISION HISTORY:
;   2003-03-03  written - Blanton Roweis and Hogg
;-
pro kmeans, ngroup, data, group, group_mean=group_mean, $
            group_covar=group_covar, seed=seed,quiet=quiet

if(n_elements(size(data,/dimensions)) eq 1) then $
  ndim=1 $
else $
  ndim=(size(data,/dimensions))[0]
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
        for jj=0L, ngroup-1L do begin
            delta=data[*,ii]-group_mean[*,jj]
            error_group[jj]=transpose(delta)#delta
        endfor
        minerror[ii]=min(error_group,jjminerror)
        total_error=total_error+minerror[ii]
        if(group[ii] ne jjminerror) then begin
            converged=0
            group[ii]=jjminerror
        endif
    endfor
    if NOT keyword_set(quiet) then splog,'total_error= '+string(total_error)

;   update group centers
    if(converged eq 0) then begin
        for jj=0L, ngroup-1L do begin
            indx_group=where(group eq jj, indx_count)
            if(indx_count gt 0) then begin
                group_mean[*,jj]=total(data[*,indx_group],2,/double)/ $
                  double(indx_count)
            endif else begin
                maxminerror=max(minerror,iiminerror)
                group_mean[*,jj]=data[*,ii]
                minerror[iiminerror]=0.D
            endelse
        endfor
    endif
endwhile

epsilon=1.e-2
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
