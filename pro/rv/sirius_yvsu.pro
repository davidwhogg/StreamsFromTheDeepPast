;+
;   NAME:
;      sirius_yvsu
;   PURPOSE:
;      plot the y vs U diagram for the sirius clump from gcs data
;   CALLING SEQUENCE:
;   INPUT:
;      k   - number of Gaussians
;      w   - regularization parameter
;   OPTIONAL INPUT:
;      sample- sample to use (default ALMS=6)
;      clumpcut - hard cut for clump_membership (default .5)
;   KEYWORDS:
;      pleiades - do this for pleiades
;      hyades   - do this for hyades
;      ngc1901  - do this for ngc1901
;      cmd      - also plot color magnitude diagram
;   OUPUT:
;   REVISION HISTORY:
;      2008-12-16 - written Bovy;   
;-
PRO SIRIUS_YVSU, k, w, sample=sample, clumpcut=clumpcut, $
                 plotfilename=plotfilename, pleiades=pleiades, $
                 hyades=hyades, ngc1901=ngc1901, cmd=cmd

IF ~keyword_set(sample) THEN sample= 6
IF ~keyword_set(clumpcut) THEN clumpcut= .5

;;Restore reconstruction of the velocity distribution
savefilename= './ALMSv2/fitv'+strtrim(string(K),2)+$
  '_bestV_V'+strtrim(string(w,format='(f4.1)'),2)+$
  '_sample'+strtrim(string(sample),2)+$
  '.sav'

IF file_test(savefilename) THEN BEGIN
    restore, savefilename
ENDIF ELSE BEGIN
    splog, 'Error: savefile of reconstruction does not exist in ', $
      savefilename
    splog, 'savefile must exist, returning...'
    RETURN
ENDELSE

;;restore gcs
gcsfilename= '/global/data/rv/gcs/gcs.sav'
IF file_test(gcsfilename) THEN BEGIN
    restore, gcsfilename
ENDIF ELSE BEGIN
    splog, 'Error: GCS-savefile does not exist in ', $
      genevafilename
    splog, 'GCS-savefile must exist, returning...'
    RETURN
ENDELSE
;;Get the gcs data in the right form
ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
ycovar= gcs.vrvlvbc
projection= gcs.sm

;;Compute posterior probabilities
assign_clump_members, logpost, ydata, ycovar, projection, mean, covar, amp


;;Define Sirius
siriusmean= [9.,3.,-6.]
title='Sirius'
IF keyword_set(pleiades) THEN BEGIN
    siriusmean=[-18.,-23.,-4.];;Pleiades
    title='Pleiades'
ENDIF
IF keyword_set(ngc1901) THEN BEGIN
    siriusmean=[-22.,-10.,-7.];;NGC 1901
    title='NGC 1901'
ENDIF
IF keyword_set(hyades) THEN BEGIN
    siriusmean=[-41.,-19.,-1.];;Hyades
    title='Hyades'
ENDIF
;;find the nearest Gaussian to Sirius
meandiffs= dblarr(3,k)
diffs= dblarr(k)
FOR kk= 0L, k-1 DO BEGIN
    meandiffs[*,kk]= mean[*,kk]-siriusmean
    diffs[kk] = norm(meandiffs[*,kk])
ENDFOR

mindiff= min(diffs,siriusindx)


siriusmembers= where(exp(logpost[siriusindx,*]) GE clumpcut)
nmembers= n_elements(siriusmembers)

title+= ', n='+strtrim(string(nmembers),2)

;;plot y vs u for these
;;first calculate y and u
degtorad= !DPI/180.
y= gcs[siriusmembers].radius*cos(gcs[siriusmembers].b*degtorad)*$
  sin(gcs[siriusmembers].l*degtorad)
u= dblarr(3,nmembers)
FOR ii=0L, nmembers-1 DO u[*,ii]= special_invert(gcs[siriusmembers[ii]].sm)##$
  ydata[*,siriusmembers[ii]]

djs_plot, u[0,*], u[1,*], psym=3

wait, 1

djs_plot, u[0,*], u[2,*], psym=3

wait, 1

u= u[0,*]

IF keyword_set(plotfilename) THEN k_print, filename=plotfilename

djs_plot, u, y/30.857E12, psym=5, ytitle='y [pc]', xtitle='v_x [km s^{-1}]', $
  title=title, yrange=[-200,200]

IF keyword_set(plotfilename) THEN k_end_print


t= regress(u,y,chisq=chisq,measure_errors=1/exp(logpost[siriusindx,siriusmembers]))

splog, 'Chisq= '+string(chisq)

print, t*220/30.857E12

djs_plot, gcs[siriusmembers].bvcolor,gcs[siriusmembers].vmag-5*alog10(gcs[siriusmembers].radius/30.857E12)+5,$
  xtitle='B-V',ytitle='M_V', psym=3, yrange=[6,0]


stop


END
