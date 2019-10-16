;+
;   NAME:
;      bovy_dec_to_damas
;   PURPOSE:
;      convert dec in deg to degree, am, as
;   CALLING SEQUENCE:
;      as= bovy_dec_to_damas(dec,/arcsec)
;   INPUT:
;      dec  - dec in deg
;   KEYWORDS:
;      arcmin - if set, return am
;      arcsec - if set, return as
;   OUTPUT:
;      d (or am or as)
;   REVISION HISTORY:
;      2009-01-18 - Written Bovy (NYU)
;-
FUNCTION bovy_dec_to_damas, dec, arcmin=arcmin, arcsec=arcsec
int_dec= dec

;;Calculate d am as for int_int_dec given in deg, returns am if /arcmin, as if /arcsec
IF bovy_sign(int_dec) EQ -1 THEN BEGIN
    deg= ceil(int_dec)
    int_dec= deg-int_dec
ENDIF ELSE BEGIN
    deg= floor(int_dec)
    int_dec= int_dec-deg
ENDELSE
int_dec*= 60.
am= floor(int_dec)
int_dec-= am
as= int_dec*60

IF keyword_set(arcsec) THEN BEGIN
    return_value= as
ENDIF ELSE IF keyword_set(arcmin) THEN BEGIN
    return_value= am
ENDIF ELSE Begin
    return_value= deg
ENDELSE

RETURN, return_value
END
