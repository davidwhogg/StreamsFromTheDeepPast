;+
; NAME:
;   find_all_tuples
; PURPOSE:
;   find all N-tuples consistent with a single velocity vector
; INPUTS:
;   hip       - hipparcos data structure (from read_hipparcos)
;   ntup      - N
;   nsigma    - threshold for "consistency"
; OPTIONAL INPUTS:
;   ntup_list - list of starting point M-tuple index lists
; OUTPUTS:
;   ntup_list - list of N-tuple index lists
; REVISION HISTORY:
;   2003-02-18  written - Roweis and Hogg
;-
function find_all_tuples,hip,ntup,nsigma,ntup_list=ntup_list

; setup parameters
nhip= n_elements(hip)
if not keyword_set(ntup_list) then ntup_list= lindgen(nhip)
dims= size(ntup_list,/dimensions)
if n_elements(dims) EQ 1 then start_ntup= 2L else start_ntup=dims[0]

; allocate arrays
viv_uvw= dblarr(nhip)
v_uvw= dblarr(3,nhip)
iv_uvw= dblarr(3,nhip)
icov_uvw= dblarr(3,3,nhip)

; set the prior radial velocity and its prior variance
vrad_assumed= 0D0
sigma_vrad_assumed= 3D2

; loop over data to get things into coord system
for ii=0L,nhip-1 do begin

; make square matrices from the non-square matrices by adding stuff
; the operation "uvw = RR # rlb" changes r,l,b vector to U,V,W
    ll= hip[ii].l*!DPI/1.80D2
    bb= hip[ii].b*!DPI/1.80D2
    RR= dblarr(3,3)
    RR[*,1:2]= hip[ii].nsm
    RR[0,0]= cos(bb)*cos(ll)
    RR[1,0]= cos(bb)*sin(ll)
    RR[2,0]= sin(bb)

; build mean and covariance for each star in UVW
    v_uvw[*,ii]= RR # [vrad_assumed, hip[ii].vl, hip[ii].b]
    icov_rlb= dblarr(3,3)
    icov_rlb[0,0]= sigma_vrad_assumed^2
    icov_rlb[1:2,1:2]= hip[ii].vlvbc
    icov_rlb= invert(icov_rlb)
    icov_uvw[*,*,ii]= RR # icov_rlb # transpose(RR)

; pre-cache icovs times means etc
    iv_uvw[*,ii]= icov_uvw[*,*,ii] # v_uvw[*,ii]
    viv_uvw[ii]= transpose(v_uvw[*,ii]) # iv_uvw[*,ii]
endfor

; start master loop
for tt=start_ntup,ntup do begin

; take all the (ntup-1)-tuples and find the ntup-tuples
    ntup_list= increment_tuples(viv_uvw,iv_uvw,icov_uvw,ntup_list,nsigma, $
                                failed=failed)
    if failed then splog, 'failed at ntup = ',tt
endfor

; done
return, ntup_list
end
