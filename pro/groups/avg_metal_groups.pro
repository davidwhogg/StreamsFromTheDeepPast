;+
;   NAME:
;      avg_metal_groups
;   PURPOSE:
;      calculate the average metallicity and spread of the moving
;      groups, by weighting by membership probablities
;   INPUT:
;      datafilename - filename that holds the savefile that has the
;                     data
;      solutionfilename - filename that holds the best fit mixture
;   KEYWORDS:
;      colorcut - make colorcut
;   OUTPUT:
;   HISTORY:
;      2009-12-02 - Written - Bovy (NYU)
;      2010-03-09 - Added spread - Bovy
;-
PRO AVG_METAL_GROUPS, solutionfilename=solutionfilename, $
                      datafilename=datafilename, $
                      colorcut=colorcut

IF ~keyword_set(solutionfilename) THEN $
  solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN BEGIN
    IF keyword_set(colorcut) THEN BEGIN
        datafilename='phot_gcs_0.35_BV_0.95.sav'
    ENDIF ELSE BEGIN
        datafilename='phot_gcs.sav'
    ENDELSE
ENDIF
restore, filename=datafilename
ngcs= n_elements(gcs.hip)

ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
ycovar= gcs.vrvlvbc
projection= gcs.sm
restore, filename=solutionfilename
;;First compute all the probabilities
assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
grouplogpost[3,*]= alog(exp(grouplogpost[3,*])+exp(grouplogpost[4,*]))

splog, 'Average metallicities'
;;NGC1901
splog, 'NGC1901'
splog, '--------'
mean= total(gcs.feh*exp(grouplogpost[0,*]))/total(exp(grouplogpost[0,*]))
splog, mean
splog, sqrt(total(gcs.feh^2.*exp(grouplogpost[0,*]))/total(exp(grouplogpost[0,*]))-mean^2.)
;;Sirius
splog, 'Sirius'
splog, '--------'
mean= total(gcs.feh*exp(grouplogpost[1,*]))/total(exp(grouplogpost[1,*]))
splog, mean
splog, sqrt(total(gcs.feh^2.*exp(grouplogpost[1,*]))/total(exp(grouplogpost[1,*]))-mean^2.)
;;Pleiades
splog, 'Pleiades'
splog, '--------'
mean= total(gcs.feh*exp(grouplogpost[3,*]))/total(exp(grouplogpost[3,*]))
splog, mean
splog, sqrt(total(gcs.feh^2.*exp(grouplogpost[3,*]))/total(exp(grouplogpost[3,*]))-mean^2.)
;;Hyades
splog, 'Hyades'
splog, '--------'
mean= total(gcs.feh*exp(grouplogpost[6,*]))/total(exp(grouplogpost[6,*]))
splog, mean
splog, sqrt(total(gcs.feh^2.*exp(grouplogpost[6,*]))/total(exp(grouplogpost[6,*]))-mean^2.)
;;Hercules
splog, 'Hercules'
splog, '--------'
mean= total(gcs.feh*exp(grouplogpost[7,*]))/total(exp(grouplogpost[7,*]))
splog, mean
splog, sqrt(total(gcs.feh^2.*exp(grouplogpost[7,*]))/total(exp(grouplogpost[7,*]))-mean^2.)
END
