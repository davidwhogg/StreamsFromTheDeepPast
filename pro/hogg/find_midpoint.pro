;+
; NAME:
;   v= find_midpoints(hip1,hip2,[dv=dv])
; PURPOSE:
;   find the best-guess 3-space velocity given two stars, assuming that they
;   are at the exact same velocity
; INPUTS:
;   hip1,hip2   - hipparcos structures for the two stars
; OUTPUTS:
;   v           - [u,v,w] "best" velocity
; OPTIONAL OUTPUTS:
;   dv          - magnitude of the difference between the best velocity and the
;                   closest allowed velocity of each star
; BUGS:
; REVISION HISTORY:
;   2003-02-17  written - Hogg
;-
function find_midpoint,hip1,hip2,dv=dv

; define tangential velocity three-vectors in UVW
c1= hip1.nsm # [hip1.vl,hip1.vb]
c2= hip2.nsm # [hip2.vl,hip2.vb]

; define radial velocity unit-three-vectors in UVW
ll1= hip1.l*!DPI/1.80D2
bb1= hip1.b*!DPI/1.80D2
a1= dblarr(3)
a1[0]=  cos(bb1)*cos(ll1)
a1[1]=  cos(bb1)*sin(ll1)
a1[2]=  sin(bb1)
ll2= hip2.l*!DPI/1.80D2
bb2= hip2.b*!DPI/1.80D2
a2= dblarr(3)
a2[0]=  cos(bb2)*cos(ll2)
a2[1]=  cos(bb2)*sin(ll2)
a2[2]=  sin(bb2)

; get dot products, etc
a1a2= (transpose(a1)#a2)[0]
a1c2= (transpose(a1)#c2)[0]
a2c1= (transpose(a2)#c1)[0]
determ= 1D0/(1D0-a1a2^2)

; make velocities and midpoint
v1= determ*(a1c2+a1a2*a2c1)*a1+c1
v2= determ*(a1a2*a1c2+a2c1)*a2+c2
dv= 0.5*(v1-v2)
dv= sqrt((transpose(dv)#dv)[0])
v= 0.5*(v1+v2)
return, v
end
