;+
;
;-
pro blob_mixtures_oned, ww,yy,SS

; initialize
phi= 0.0
totww= total(ww)
mm= total(ww*yy)/totww
VV= total(ww*(yy-mm)^2)/totww

iteration= 0
repeat begin
    iteration= iteration+1

; compute likelihood
    TT= VV+SS
    resid= yy-mm
    oldphi= phi
    phi= 0.5*total(ww*alog(TT))-0.5*total(ww*resid^2/TT)
    dphi= (phi-oldphi)/phi
    splog, dphi,phi,mm,VV

; e-step -- see STR documentation
    bb= mm+VV*resid/TT
    BBm= VV-VV^2/TT

; m-step
    mm= total(ww*bb)/totww
    VV= total(ww*((mm-bb)^2+BBm))/totww

endrep until (dphi LT 1D-8) and (iteration GT 2)
end
