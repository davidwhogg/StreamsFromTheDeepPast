;+
;   NAME:
;      test_phot_plx
;   PURPOSE:
;      test the function 'phot_plx'
;   CALLING SEQUENCE:
;   INPUT:
;   OUTPUT:
;      prints some stuff
;   REVISION HISTORY:
;      2009-11-09 - Written - Bovy (NYU)
;-
PRO TEST_PHOT_PLX

;;Read an isochrone, I don't care which one
datafilename='pleiades_isochrone.dat'
;;Read the isochrone
nfield= 16L
; most data types are double
fieldtypes= lonarr(nfield)+5
; some are long
fieldtypes[[12]]= 3
; give dumb names
fieldnames= 'H'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: string(9B), missingvalue: 0.0, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
isochrone= read_ascii(datafilename,num_records=nline,template=template)
iso= dblarr(2,40)
iso[0,*]= isochrone.h08[0:39]
iso[1,*]= isochrone.h09[0:39]
;;Now predict all of the Hipparcos plxs
restore, filename='../lsr2/hip2-aumer.sav'
nhip= n_elements(hip.hip)
plxs= dblarr(nhip)
FOR ii=0L, nhip-1 DO Begin
    BV= hip[ii].bt-hip[ii].vt
    V= hip[ii].vt
;    IF BV GE 1.5 OR BV LE .2 THEN CONTINUE
    plxs[ii]= phot_plx(BV,V,iso)
ENDFOR
stop

END
