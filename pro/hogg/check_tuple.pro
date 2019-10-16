;+
; NAME:
;   check_tuple
; PURPOSE:
;   check that an N-tuple of 3-d gaussians has a best point within nsigma of
;     all N gaussians
; INPUTS:
;   viv      - vector of (mean^T . covar^(-1} . mean)
;   iv       - vector of (covar^(-1} . mean)
;   icov     - vector of (covar^(-1))
;   tuple    - list of indices for the N-tuple
;   nsigma   - threshold
; OUTPUTS:
;   accept   - binary yes/no
; REVISION HISTORY:
;   2003-02-18  written - Roweis and Hogg
;-
function check_tuple, viv,iv,icov,tuple,nsigma

ntup= n_elements(tuple)
dimen= (size(icov,/dimensions))[0]
cbar= dblarr(dimen,dimen)
vbar= dblarr(dimen)
for ii=0L,ntup-1 do begin
    cbar= cbar+icov[*,*,tuple[ii]]
    vbar= vbar+iv[*,tuple[ii]]
endfor
vbar= invert(cbar) # vbar
accept= 1B
ii= 0L
while accept EQ 1 AND ii LT ntup do begin
    if (viv[tuple[ii]]-2D0*transpose(vbar)#iv[*,tuple[ii]]+ $
      transpose(vbar)#icov[*,*,tuple[ii]]#vbar) GT nsigma^2 then accept=0B
    ii= ii+1
endwhile
return, accept
end
