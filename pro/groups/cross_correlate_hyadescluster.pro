;+
;   NAME:
;      cross_correlate_hyadescluster.pro
;   PURPOSE:
;      see how many of the stars in our sample are Hyades members
;   INPUT:
;   OUTPUT:
;   HISTORY:
;      2010-03-05 - Written - Bovy (NYU)
;-
PRO CROSS_CORRELATE_HYADESCLUSTER, hyadessavefilename=hyadessavefilename, $
                                   samplefilename=samplefilename
IF ~keyword_set(hyadessavefilename) THEN hyadessavefilename= 'hyadesMembers_full.sav'
IF ~keyword_set(samplefilename) THEN samplefilename= '../lsr2/hip2-aumer.sav'

IF ~file_test(hyadessavefilename) THEN BEGIN
    hyades= read_hyades(savefilename=hyadessavefilename,/single,/tenpc)
ENDIF ELSE BEGIN
    restore, hyadessavefilename
ENDELSE

restore, samplefilename

nhyades=n_elements(hyades.hip)
insample= bytarr(nhyades)
FOR ii=0L, nhyades-1 DO BEGIN
    match= where(hyades[ii].hip EQ hip.hip)
    IF match[0] EQ -1 THEN insample[ii]= 0B ELSE insample[ii]= 1B
ENDFOR
print, n_elements(where(insample EQ 1B))
print, n_elements(hip.hip)
END
