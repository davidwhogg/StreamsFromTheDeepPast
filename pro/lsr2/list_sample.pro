;+
;   NAME:
;      list_sample
;   PURPOSE:
;      write a list of hip numbers to send to Aumer for comparison
;   CALLING SEQUENCE:
;   INPUT:
;      filename - 
;   OUTPUT:
;   KEYWORDS:
;   REVISION HISTORY:
;      2009-08-04 - Written - Bovy (NYU)
;-
PRO LIST_SAMPLE, filename=filename

IF ~keyword_set(filename) THEN filename='BovySample.dat'

;;Restore sample
restore, filename='hip2-aumer.sav'

hipn= hip.hip
;;Sort the hipn
hipn= hipn[sort(hipn)]

;;Write to file
OPENW, wlun, filename, /GET_LUN
FOR ii= 0L, n_elements(hip.hip)-1 DO BEGIN
    PRINTF, wlun, strtrim(string(hipn[ii]),2)
ENDFOR
FREE_LUN, wlun

END
