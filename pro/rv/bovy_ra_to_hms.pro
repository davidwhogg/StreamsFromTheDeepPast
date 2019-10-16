;+
;   NAME:
;      bovy_ra_to_hms
;   PURPOSE:
;      convert ra in deg to hms
;   CALLING SEQUENCE:
;      s= bovy_ra_to_hms(ra,/sec)
;   INPUT:
;      ra  - ra in deg
;   KEYWORDS:
;      min - if set, return m
;      sec - if set, return s
;   OUTPUT:
;      h (or m or s)
;   REVISION HISTORY:
;      2009-01-18 - Written Bovy (NYU)
;-
FUNCTION bovy_ra_to_hms, ra, min=min, sec=sec
int_ra= ra

;;Calculate h m s for int_ra given in deg, returns m if /m, s if /s
h= floor(int_ra/15.)
int_ra-= 15*h
m= floor(int_ra)
int_ra-= m
m*= 4
int_ra*= 60.
m_tmp= floor(int_ra/15.)
int_ra-= m_tmp*15
m+= m_tmp
s= floor(int_ra)
int_ra-= s
s*= 4.
int_ra*= 60.
s+= int_ra/15.

IF keyword_set(sec) THEN BEGIN
    return_value= s
ENDIF ELSE IF keyword_set(min) THEN BEGIN
    return_value= m
ENDIF ELSE Begin
    return_value= h
ENDELSE

RETURN, return_value
END
