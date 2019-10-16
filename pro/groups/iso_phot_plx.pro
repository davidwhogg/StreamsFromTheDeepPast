;+
;   NAME:
;      is_phot_plx
;   PURPOSE:
;      calculate the log probability of a parallax given an isochrone
;   INPUT:
;   INPUT:
;      plx - parallax to evaluate the probability of (can be an array)
;      bvcolor - B-V color of the plx star
;      vmag - apparent magnitude of the plx star
;      sig_plx - uncertainty on the plx
;      iso - isochrone [B,V] dblarr(2,N)
;      alpha - alpha parameter
;      back_loglike - background model's likelihood
;      sspspread - spread to add to the SSP prediction
;   OUTPUT:
;      p(plx|foreground model)
;   REVISION HISTORY:
;      2009-11-10 - Written - Bovy (NYU)
;-
FUNCTION ISO_PHOT_PLX, plx, bvcolor,vmag,sig_plx, iso,$
                       alpha=alpha, $
                       back_loglike=back_loglike,sspspread=sspspread
modelplx= phot_plx(bvcolor,vmag,iso,/ms)
IF modelplx EQ -1. THEN BEGIN
    RETURN, back_loglike+alog(alpha)
ENDIF ELSE BEGIN
    thisspread= sig_plx^2.+sspspread^2.*plx^2.
    RETURN, alog((1.-alpha)/sqrt(2.*!DPI*thisspread)*exp(-0.5*(plx-modelplx)^2./thisspread)+alpha*exp(back_loglike))
ENDELSE
END
