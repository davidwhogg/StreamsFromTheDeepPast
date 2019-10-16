;+
;   NAME:
;      bovy_in_range_periodic
;   PURPOSE:
;      determine whether a value is in a given range, when the values
;      are part of a circle [i.e., the values wrap around]
;   CALLING SEQUENCE:
;      bovy_in_range_periodic(value,range,orange)
;   INPUT:
;      value    - value which is in the range or not (in orange)
;      range    - range (one of these should be in orange)
;      orange - overall range, like [0,360.]
;   KEYWORDS:
;   OUTPUT:
;      1 if the value's in the range, 0 otherwise
;   REVISION HISTORY:
;      2009-01-20 - Written Bovy (NYU)
;-
FUNCTION BOVY_IN_RANGE_PERIODIC, value, range, orange

inrange= 0B

IF (range[0] LT orange[0] AND range[1] GT orange[1]) THEN RETURN, 1B

IF (range[0] GE orange[0] AND range[1] LE orange[1]) THEN BEGIN
    IF (value GE range[0] AND value LE range[1]) THEN inrange= 1B ELSE inrange= 0B
ENDIF ELSE IF (range[0] LT orange[0]) THEN BEGIN
    IF ((value GE orange[0] AND value LE range[1]) OR $
        (value GE (range[0]+(orange[1]-orange[0])) AND value LE orange[1])) THEN $
        inrange= 1B ELSE inrange= 0B
ENDIF ELSE IF (range[1] GT orange[1]) THEN BEGIN
    IF ((value GE range[0] AND value LE orange[1]) OR $
        (value GE orange[0] AND (value LE (range[1]-(orange[1]-orange[0]))))) THEN $
      inrange= 1B ELSE inrange= 0B
ENDIF


RETURN, inrange

END
