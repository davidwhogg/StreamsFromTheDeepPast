;+
;   NAME:
;      correlate_gcs_hip
;
;   PURPOSE:
;      correlate the gcs and dehnen's hipparcos sample
;
;   CALLING SEQUENCE:
;      hip_gcs=correlate_gcs_hip(indx=indx,/compare_pm)
;
;   OUTPUT:
;      hip_gcs  - Array of HIP numbers that are in both samples
;
;   OPTIONAL OUTPUT:
;      indx     - list with correlated indices (dblarr(2,ncorr))
;
;   REVISION HISTORY:
;      2008-09-03 - Written Bovy
;-
FUNCTION CORRELATE_GCS_HIP, indx=indx, compare_pm=compare_pm

;;restore hip-dehnen
restore, '/global/data/hipparcos/hip-dehnen.sav'
restore, '/global/data/rv/gcs/gcs.sav'

nhip= n_elements(hip.hip)
hip_gcs= dblarr(nhip)
FOR ii= 0L, nhip-1 DO IF (where(gcs.hip EQ hip[ii].hip))[0] NE -1 THEN hip_gcs[ii]= hip[ii].hip

hip_gcs= hip_gcs[where(hip_gcs NE 0)]

IF arg_present(indx) OR keyword_set(compare_pm) THEN BEGIN
    ncorr= n_elements(hip_gcs)
    indx= dblarr(2,ncorr)
    FOR ii= 0L, ncorr-1 DO indx[*,ii]= [(where(hip.hip EQ hip_gcs[ii]))[0], (where(gcs.hip EQ hip_gcs[ii]))[0]]
    IF keyword_set(compare_pm) THEN BEGIN
        k_print, filename='hip_gcs_pm.ps'
        hogg_plothist, sqrt(hip.pmra^2+hip.pmdec^2), xrange=[0,1000], $
          /totalweight, xtitle='!4I!7l!4I!x'+textoidl(' [mas yr^{-1}]'), $
          ytitle='number of stars', yrange=[0,(3000)],npix=40
        hogg_plothist, sqrt(hip[indx[0,*]].pmra^2+hip[indx[0,*]].pmdec^2), $
          /overplot, /totalweight, linestyle=2, xrange=[0,1000],npix=40
;;make inset
        bangP= !P
        !P.POSITION=[0.4,0.4,0.8,0.8]
        !P.NOERASE=1
        hogg_plothist, sqrt(hip.pmra^2+hip.pmdec^2), xrange=[400,1000], $
          /totalweight, yrange=[0,(60)],npix=600/25, xcharsize=1.2, $
          ycharsize=1.2,xminor=1, yminor=1 ;, position=,/noerase
        hogg_plothist, sqrt(hip[indx[0,*]].pmra^2+hip[indx[0,*]].pmdec^2), $
          /overplot, /totalweight, linestyle=2, xrange=[0,1000],npix=40, $
          xcharsize=1.2, ycharsize=1.2, xminor=1
        !P=bangP
        k_end_print
    ENDIF
ENDIF

RETURN, hip_gcs

END
