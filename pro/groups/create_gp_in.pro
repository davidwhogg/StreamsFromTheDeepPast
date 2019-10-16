;+
;
;-
PRO CREATE_GP_IN, filename=filename

restore, filename='../lsr2/hip2-aumer.sav'
nhip= n_elements(hip.hip)

out= create_struct('bvcolor',dblarr(nhip),'MV',dblarr(nhip),'e_MV',dblarr(nhip))
FOR ii=0L, nhip-1 DO BEGIN
    out.bvcolor[ii]= hip[ii].bvcolor
    out.MV[ii]= hip[ii].vmag-5D0*alog10(100/hip[ii].plx)
    out.e_MV[ii]= 5D0/alog(10D)*hip[ii].e_plx/hip[ii].plx
ENDFOR

mwrfits, out, filename

END
