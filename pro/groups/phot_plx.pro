;+
;   NAME:
;      phot_plx
;   PURPOSE:
;      given an isochrone and photometry, predict the plx
;   CALLING SEQUENCE:
;      plx= phot_plx(BV,V,iso)
;   INPUT:
;      BV  - color
;      V   - apparent mag
;      iso - isochrone [B,V] dblarr(2,N)
;   KEYWORDS:
;      ms - isolate the main-sequence explicitly
;   OUTPUT:
;      plx - photometric plx
;      absmag - the absolute magnitude the star should have had
;   REVISION HISTORY:
;      2009-09-11 - Written - Bovy (NYU)
;-
FUNCTION PHOT_PLX, BV, V, iso, absmag=absmag, ms=ms
;;isochrone color=
isoBV= iso[0,*]-iso[1,*]
;;Find where BV is on the isochrone
foundBV= 0
end_ms_is_near= 0 ;;What follows are some pretty ad-hoc criteria to signal the end of the MS, hopefully these work
ii=0
jj=0
WHILE ~foundBV AND ii LT (n_elements(isoBV[0,*])-1) AND jj LT 15 DO BEGIN
    IF isoBV[ii] GE BV AND isoBV[ii+1] LE BV THEN BEGIN
        foundBV= 1
        BREAK
    ENDIF
    IF end_ms_is_near EQ 0 AND keyword_set(ms) AND isoBV[ii] LT isoBV[ii+5] THEN BEGIN
        end_ms_is_near= 1
    ENDIF
    ii+= 1
    IF end_ms_is_near EQ 1 THEN jj+= 1
ENDWHILE
IF ii EQ (n_elements(isoBV[0,*])-1) OR foundBV EQ 0 THEN RETURN, -1
;;Now interpolate between the BV to find absV
IF isoBV[ii+1] EQ isoBV[ii] THEN BEGIN
    absV= iso[1,ii]
ENDIF ELSE BEGIN
    absV= (iso[1,ii]-iso[1,ii+1])/(isoBV[ii]-isoBV[ii+1])*(BV-isoBV[ii])+iso[1,ii]
ENDELSE
absmag= absV
;;And calculate the plx
plx= 10.^(-(V-absV)/5.+2.)
RETURN, plx
END
