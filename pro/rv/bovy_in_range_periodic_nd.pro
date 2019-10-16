;+
;   NAME:
;      bovy_in_range_periodic_nd
;   PURPOSE:
;      determine whether a value is in a given N-d range, when the
;      values are part of a sphere (that is, they wrap around)
;   CALLING SEQUENCE:
;   INPUT:
;      value    - value which is in the range or not (in orange) dblarr(N)
;      range    - range (one of these should be in orange) dblarr(2,N)
;      orange - overall range, like [[0,360.],[-90.,90.]] dblarr(2,N)
;   KEYWORDS:
;   OUTPUT:
;      1 if the value's in the range, 0 otherwise
;   REVISION HISTORY:
;      2009-01-21 - Written Bovy (NYU)
;-
FUNCTION BOVY_IN_RANGE_PERIODIC_ND,value, range, orange

d= n_elements(value)
in_1d_range= dblarr(d)

FOR ii= 0L, d-1 DO $
  in_1d_range[ii]= bovy_in_range_periodic(value[ii],range[*,ii],orange[*,ii])

return_value= 1B
FOR ii= 0L, d-1 DO return_value*= in_1d_range[ii]

RETURN, byte(return_value)

END
