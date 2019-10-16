;+
;   NAME:
;      calc_probable_mg_members
;   PURPOSE:
;      calculate which are likely moving group members
;   INPUT:
;      group - which group?
;      maxindx - go down this far in the list of probable members
;   KEYWORDS:
;      /gcs - use gcs
;   OUTPUT:
;      prints sorted list
;   HISTORY:
;      2009-11-17 - Written - Bovy (NYU)
;-
PRO CALC_PROBABLE_MG_MEMBERS, group, gcs=gcs, maxindx=maxindx
IF ~keyword_set(maxindx) THEN maxindx= 10L
IF keyword_set(gcs) THEN BEGIN
    datafilename='gcs.sav'
ENDIF ELSE BEGIN
    datafilename='../lsr2/hip2-aumer.sav'
ENDELSE
restore, filename=datafilename
IF keyword_set(gcs) THEN hip= gcs
nhip= n_elements(hip.hip)
solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
IF keyword_set(gcs) THEN BEGIN
    ydata= transpose([[hip.vr]])
    ycovar= hip.vrc
    ycovar= reform(ycovar,1,1,nhip)
    projection= hip.sm[*,0]
    projection= reform(projection,3,1,nhip)
ENDIF ELSE BEGIN
    ydata= transpose([[hip.vl], [hip.vb]])
    ycovar= [hip.vlvbc]
    projection= [hip.nsm]
ENDELSE
;;Restore solution
restore, filename=solutionfilename
;;First compute all the probabilities
assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
IF group EQ 3 THEN BEGIN
    grouplogpost[group,*]= alog(exp(grouplogpost[group,*])+exp(grouplogpost[group+1,*]))
ENDIF

;;Sort by probability
sortindx= reverse(sort(grouplogpost[group,*]))
print, sortindx[0:maxindx], exp(grouplogpost[group,sortindx[0:maxindx]])
END
