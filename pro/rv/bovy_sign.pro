;+
;   NAME:
;      bovy_sign
;
;   PURPOSE:
;      sign-function
;
;   CALLING SEQUENCE:
;      s = sign(x)
;
;   INPUT:
;      x    : number
;
;   OUTPUT:
;     s     : +1 if positive, -1 if negative, +1 if 0
;
;   REVISION HISTORY:
;      07-10-08 - Written Bovy
;-
FUNCTION BOVY_SIGN, x

RETURN, 1 - (x LT 0)*2

END
