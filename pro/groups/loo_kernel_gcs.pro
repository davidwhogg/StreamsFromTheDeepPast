;+
;   NAME:
;      loo_kernel_gcs
;   PURPOSE:
;      calculate the leave-one-out log likelihood of a certain kernel
;      width for the GCS sample
;   INPUT:
;      width - width parameter of the kernel function
;      zcut - metallicity at which to cut the sample in half
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      silent - be silent
;      highz - use high metallicity subsample
;      lowz - use low metallicity subsample
;   OUTPUT:
;      log pseudo-likelihood
;   HISTORY:
;      2009-11-05 - Adapted from loo_kernel - Bovy (NYU)
;-
FUNCTION loo_kernel_gcs, width, kernel=kernel, silent=silent,highz=highz,zcut=zcut, lowz=lowz
IF ~keyword_set(zcut) THEN zcut= -0.132429
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN datafilename='phot_gcs_0.35_BV_0.95.sav'
restore, filename=datafilename
IF keyword_set(highz) THEN BEGIN
    indx= where(gcs.feh GT zcut)
    gcs= gcs[indx]
ENDIF ELSE IF keyword_set(lowz) THEN BEGIN
    indx= where(gcs.feh LT zcut)
    gcs= gcs[indx]
ENDIF
ngcs= n_elements(gcs.hip)
logls= dblarr(ngcs)
FOR ii=0L, ngcs-1 DO BEGIN
    ;;Create loo data set
    indx= where(gcs.hip NE gcs[ii].hip)
    data= dblarr(ngcs-1,2)
    data[*,0]= gcs[indx].bvcolor
    data[*,1]= gcs[indx].vmag-5.*alog10(100./gcs[indx].plx)
    logls[ii]= alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                                    gcs[ii].vmag,$
                                    gcs[ii].e_plx,data,$
                                    width=width,kernel=kernel,silent=silent))
ENDFOR
nonzero= where(logls NE -(machar(/double)).xmax AND FINITE(logls))
minnonzero= min(logls[nonzero])
zero=where(logls EQ -(machar(/double)).xmax OR FINITE(logls,/NaN) OR FINITE(logls,/INFINITY))
IF zero[0] NE -1 THEN logls[zero]= minnonzero
RETURN, total(logls)
END
