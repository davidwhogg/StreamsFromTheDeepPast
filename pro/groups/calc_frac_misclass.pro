;+
;   NAME:
;      calc_frac_miscat
;   PURPOSE:
;      calculate the fraction of the group that is possibly misclassified
;   CALLING SEQUENCE:
;      alpha=calc_frac_miscat(group,pcut=pcut)
;   INPUT:
;      group - which moving group
;      pcut  - where is the p_ij reliable (default: 0.5)
;      savefilename - filename with the velocity field
;      datafilename - filename with the data
;   OUTPUT:
;      fraction possibly misclassified
;   REVISION HISTORY:
;      2009-09-15 - Written - Bovy (NYU)
;-
FUNCTION CALC_FRAC_MISCLASS, group, pcut=pcut, savefilename=savefilename, $
                             datafilename=datafilename
IF ~keyword_set(pcut) THEN pcut= 0.5
;;Restore data
IF ~keyword_set(savefilename) THEN $
  savefilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
;;Restore data
IF ~keyword_set(datafilename) THEN datafilename='../lsr2/hip2-aumer.sav'
restore, filename=datafilename
nhip= n_elements(hip.hip)
ydata= transpose([[hip.vl], [hip.vb]])
ycovar= hip.vlvbc
projection= hip.nsm
;;Restore solution
restore, filename=savefilename
;;First compute all the probabilities
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
IF group EQ 3 THEN BEGIN;;Pleiades is a special case
    logpost[group,*]= alog(exp(logpost[group,*])+exp(logpost[group+1,*]))
ENDIF
;;Then calculate alpha
alpha=0.
pijs= exp(logpost[group,*])
indx_sumpijs= where(pijs LE pcut)
alpha=total(pijs[indx_sumpijs])/total(pijs)
RETURN, alpha
END
