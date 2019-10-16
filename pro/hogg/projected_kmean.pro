;+
; NAME:
;   projected_kmean
; PURPOSE:
;   K-means in projected space 
; INPUTS:
;   num - number of groups
;   y - [2,N] data points 
; OUTPUTS:
; BUGS:
;   inputs and outputs not documented; comment header incomplete.
;   weights not used or passed in 
;   we are not sure that the xcovars are calculated correctly
; REVISION HISTORY:
;   2003-02-18  written - Blanton Roweis and Hogg
;-
pro projected_kmean, ngroup, ydata, ycovar, projection, amp, $
                     xmean, xcovar, quiet=quiet

epsilon=1.d
ndimy=n_elements(ycovar)/n_elements(ydata)
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(projection)/ndata/ndimy

; set initial conditions if necessary
if(n_elements(amp) ne ngroup or n_elements(xmean) eq 0 or n_elements(xcovar)) $
  then begin
    projected_initialize, ngroup, ydata, ycovar, projection, amp, $
      xmean, xcovar
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
        invycovar=invert(ycovar[*,*,ii])
        for jj=0L, ngroup-1L do begin
            delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,jj]
            error_group[jj]=transpose(delta)#invycovar#delta
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
                denom=dblarr(ndimx,ndimx)
                numer=dblarr(ndimx,ndimx)
                for kk=0L, n_elements(indx_group)-1L do begin
                    ii=indx_group[kk]
                    projinvycovar=projection[*,*,ii]#invert(ycovar[*,*,ii])
                    denom=denom+projinvycovar#transpose(projection[*,*,ii])
                    numer=numer+projinvycovar#ydata[*,ii]
                endfor
                xmean[*,jj]=invert(denom)#numer
            endif else begin
                maxminerror=max(minerror,iiminerror)
                xmean[*,jj]=projection[*,*,iiminerror]#ydata[*,ii]
                minerror[iiminerror]=0.D
            endelse
            amp[jj]=double(indx_count)/double(ndata)
        endfor
    endif
endwhile

for jj=0L, ngroup-1L do begin
    indx_group=where(group eq jj, indx_count)
    xdata=dblarr(ndimx,indx_count)
    for ii=0L, indx_count-1L do $
      xdata[*,ii]=projection[*,*,indx_group[ii]]# $
      ydata[*,indx_group[ii]]-xmean[*,jj]
    xcovar[*,*,jj]=transpose(xdata)##xdata/double(indx_count)
    if(indx_count le ndimx) then $
        xcovar[*,*,jj]=xcovar[*,*,jj]+epsilon*identity(ndimx)
endfor
        
end
