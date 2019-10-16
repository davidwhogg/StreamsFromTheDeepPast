;+
;   NAME:
;      loo_kernel
;   PURPOSE:
;      calculate the leave-one-out log likelihood of a certain kernel width
;   INPUT:
;      width - width parameter of the kernel function
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      silent - be silent
;   OUTPUT:
;      log pseudo-likelihood
;   HISTORY:
;      2009-11-05 - Written - Bovy (NYU)
;-
FUNCTION loo_kernel, width, kernel=kernel, silent=silent
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN datafilename='../lsr2/hip2-aumer.sav'
restore, filename=datafilename
nhip= n_elements(hip.hip)
logls= dblarr(nhip)
FOR ii=0L, nhip-1 DO BEGIN
    ;;Create loo data set
    indx= where(hip.hip NE hip[ii].hip)
    data= dblarr(nhip-1,2)
    data[*,0]= hip[indx].bvcolor
    data[*,1]= hip[indx].vmag-5.*alog10(100./hip[indx].plx)
    logls[ii]= alog(kernel_phot_plx(hip[ii].plx,hip[ii].bvcolor,$
                                    hip[ii].vmag,$
                                    hip[ii].e_plx,data,$
                                    width=width,kernel=kernel,silent=silent))
ENDFOR
nonzero= where(logls NE -(machar(/double)).xmax AND FINITE(logls))
minnonzero= min(logls[nonzero])
zero=where(logls EQ -(machar(/double)).xmax OR FINITE(logls,/NaN) OR FINITE(logls,/INFINITY))
IF zero[0] NE -1 THEN logls[zero]= minnonzero
RETURN, total(logls)
END
