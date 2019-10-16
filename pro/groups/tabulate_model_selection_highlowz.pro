;+
;   NAME:
;      tabulate_model_selection_highlowz
;   PURPOSE:
;      do the high/low metallicity model selection for the moving
;      groups using GCS
;   INPUT:
;      kernel - which kernel to use
;      kernel_width_full - kernel width for the full GCS sample
;      kernel_width_highz - kernel width for the high metallicity GCS
;                           sample
;      kernel_width_lowz - kernel width for the low metallicity GCS
;                          sample
;      zcut - metallicity at which to cut the sample in half
;      solutionfilename - filename that holds the best fit mixture
;   OUTPUT:
;      For each moving group, prints the background, highz foreground,
;      and lowz foreground likelihood
;   HISTORY:
;      2009-11-21 - Written - Bovy (NYU)
;-
PRO TABULATE_MODEL_SELECTION_HIGHLOWZ, kernel=kernel, full_width=full_width, $
                                       highz_width=highz_width, $
                                       lowz_width=lowz_width, zcut=zcut, $
                                       solutionfilename=solutionfilename
IF ~keyword_set(zcut) THEN zcut= -0.132429
IF ~keyword_set(solutionfilename) THEN $
  solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN datafilename='phot_gcs_0.35_BV_0.95.sav'
restore, filename=datafilename
ngcs= n_elements(gcs.hip)

;;First calculate the background probability for each star
backloglike= dblarr(ngcs)
data= dblarr(ngcs,2)
data[*,0]= gcs.bvcolor
data[*,1]= gcs.vmag-5.*alog10(100./gcs.plx)

highz_indx= where(gcs.feh GT zcut)
highz_gcs= gcs[highz_indx]
nhighz= n_elements(highz_gcs.hip)
highz_data= dblarr(nhighz,2)
highz_data[*,0]= highz_gcs.bvcolor
highz_data[*,1]= highz_gcs.vmag-5.*alog10(100./highz_gcs.plx)

lowz_indx= where(gcs.feh LT zcut)
lowz_gcs= gcs[lowz_indx]
nlowz= n_elements(lowz_gcs.hip)
lowz_data= dblarr(nlowz,2)
lowz_data[*,0]= lowz_gcs.bvcolor
lowz_data[*,1]= lowz_gcs.vmag-5.*alog10(100./lowz_gcs.plx)


splog, "Calculating background log likelihood..."
FOR ii=0L, ngcs-1 DO BEGIN
    backloglike[ii]= alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                                          gcs[ii].vmag,$
                                          gcs[ii].e_plx,data,$
                                          width=full_width,$
                                          kernel=kernel,silent=silent))
ENDFOR
totalbackloglike= total(backloglike)

;;Then, for each moving group calculate the high and low metallicity
;;likelihood
ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
ycovar= gcs.vrvlvbc
projection= gcs.sm
restore, filename=solutionfilename
;;First compute all the probabilities
assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
grouplogpost[3,*]= alog(exp(grouplogpost[3,*])+exp(grouplogpost[4,*]))


;;NGC1901
ngc1901loglike_highz= dblarr(ngcs)
ngc1901loglike_lowz= dblarr(ngcs)
splog, "Calculating NGC1901 high and low metallicity log likelihoods..."
FOR ii=0L, ngcs-1 DO BEGIN
    ngc1901loglike_highz[ii]= exp(grouplogpost[0,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,highz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[0,ii]))*backloglike[ii]
    ngc1901loglike_lowz[ii]= exp(grouplogpost[0,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,lowz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[0,ii]))*backloglike[ii]
ENDFOR


;;Sirius
siriusloglike_highz= dblarr(ngcs)
siriusloglike_lowz= dblarr(ngcs)
splog, "Calculating SIRIUS high and low metallicity log likelihoods..."
FOR ii=0L, ngcs-1 DO BEGIN
    siriusloglike_highz[ii]= exp(grouplogpost[1,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,highz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[1,ii]))*backloglike[ii]
    siriusloglike_lowz[ii]= exp(grouplogpost[1,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,lowz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[1,ii]))*backloglike[ii]
ENDFOR


;;Pleiades
pleiadesloglike_highz= dblarr(ngcs)
pleiadesloglike_lowz= dblarr(ngcs)
splog, "Calculating PLEIADES high and low metallicity log likelihoods..."
FOR ii=0L, ngcs-1 DO BEGIN
    pleiadesloglike_highz[ii]= exp(grouplogpost[3,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,highz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[3,ii]))*backloglike[ii]
    pleiadesloglike_lowz[ii]= exp(grouplogpost[3,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,lowz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[3,ii]))*backloglike[ii]
ENDFOR


;;Hyades
hyadesloglike_highz= dblarr(ngcs)
hyadesloglike_lowz= dblarr(ngcs)
splog, "Calculating HYADES high and low metallicity log likelihoods..."
FOR ii=0L, ngcs-1 DO BEGIN
    hyadesloglike_highz[ii]= exp(grouplogpost[6,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,highz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[6,ii]))*backloglike[ii]
    hyadesloglike_lowz[ii]= exp(grouplogpost[6,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,lowz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[6,ii]))*backloglike[ii]
ENDFOR


;;Hercules
herculesloglike_highz= dblarr(ngcs)
herculesloglike_lowz= dblarr(ngcs)
splog, "Calculating HERCULES high and low metallicity log likelihoods..."
FOR ii=0L, ngcs-1 DO BEGIN
    herculesloglike_highz[ii]= exp(grouplogpost[7,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,highz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[7,ii]))*backloglike[ii]
    herculesloglike_lowz[ii]= exp(grouplogpost[7,ii])*$
      alog(kernel_phot_plx(gcs[ii].plx,gcs[ii].bvcolor,$
                           gcs[ii].vmag,$
                           gcs[ii].e_plx,lowz_data,$
                           width=highz_width,$
                           kernel=kernel,silent=silent))+$
      (1.-exp(grouplogpost[7,ii]))*backloglike[ii]
ENDFOR


;;Print the results
floatformat='(F8.1)'
splog, "Background GCS log probability"
splog, "------------------------------"
splog, totalbackloglike

splog, "NGC1901 GCS log probabilities"
splog, "------------------------------"
splog, "High metallicity: "+strtrim(string(total(ngc1901loglike_highz),format=floatformat),2)
splog, "Low metallicity: "+strtrim(string(total(ngc1901loglike_lowz),format=floatformat),2)

splog, "SIRIUS GCS log probabilities"
splog, "------------------------------"
splog, "High metallicity: "+strtrim(string(total(siriusloglike_highz),format=floatformat),2)
splog, "Low metallicity: "+strtrim(string(total(siriusloglike_lowz),format=floatformat),2)

splog, "PLEIADES GCS log probabilities"
splog, "------------------------------"
splog, "High metallicity: "+strtrim(string(total(pleiadesloglike_highz),format=floatformat),2)
splog, "Low metallicity: "+strtrim(string(total(pleiadesloglike_lowz),format=floatformat),2)

splog, "HYADES GCS log probabilities"
splog, "------------------------------"
splog, "High metallicity: "+strtrim(string(total(hyadesloglike_highz),format=floatformat),2)
splog, "Low metallicity: "+strtrim(string(total(hyadesloglike_lowz),format=floatformat),2)

splog, "HERCULES GCS log probabilities"
splog, "------------------------------"
splog, "High metallicity: "+strtrim(string(total(herculesloglike_highz),format=floatformat),2)
splog, "Low metallicity: "+strtrim(string(total(herculesloglike_lowz),format=floatformat),2)

END
