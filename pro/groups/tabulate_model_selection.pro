;+
;   NAME:
;      tabulate_model_selection
;   PURPOSE:
;      make a table with the model selection outcome of various model
;      selection tests
;   INPUT:
;      width - width parameter of the kernel
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      colorcut - if set, only use colors smaller than this value
;      grid - number of grid points per Z, age dimension (actual
;             number of grid+1, best if 300 is divisible by grid)
;      agrid - grid parameter for alpha
;      sspspread - spread to include around the ssp prediction for the
;                  parallax
;      group - which group to focus on if justonegroup is set
;   KEYWORDS:
;      justonegroup - if set, only calculate one group, don't make a table
;   OUTPUT:
;   HISTORY:
;      2009-11-17 - Written - Bovy (NYU)
;-
PRO TABULATE_MODEL_SELECTION, kernel=kernel,width=width,colorcut=colorcut, $
                              grid=grid,agrid=agrid,sspspread=sspspread, $
                              justonegroup=justonegroup,group=group

IF ~keyword_set(group) THEN group= 10
ntestspermg= 4

gcs_back_loglike_filename='gcs_loglike_back.sav'
IF ~keyword_set(justonegroup) THEN BEGIN
    ;;First calculate the background loglikelihood
    gcs_back_loglike= calc_gcs_prob(kernel=kernel,width=width,colorcut=colorcut,$
                                    backloglikefilename=gcs_back_loglike_filename)
ENDIF

;;Calculate fixed alpha, marginalized gcs probabilities for all of the
;;moving groups
IF ~keyword_set(justonegroup) OR (keyword_set(justonegroup) AND group EQ 0) THEN BEGIN
    ngc1901_basic_filename='gcs_loglike_ngc1901'
    ngc1901_gcs_loglikes= dblarr(ntestspermg)
    ngc1901_gcs_loglikes[0]= calc_gcs_prob(group=0,/foreground,/fixalpha,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_fixalpha_ngc1901.sav',$
                                           grid=grid,agrid=0,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=ngc1901_basic_filename+'_fixalpha.sav',$
                                           sspspread=sspspread)
    ngc1901_gcs_loglikes[1]= calc_gcs_prob(group=0,/foreground,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_ngc1901.sav',$
                                           grid=grid,agrid=agrid,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=ngc1901_basic_filename+'.sav',$
                                           sspspread=sspspread)
    ngc1901_gcs_loglikes[2]= calc_gcs_prob(group=0,/foreground,/fixalpha,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=ngc1901_basic_filename+'_fixalpha_bestfit.sav',$
                                           /bestfitonly)
    ngc1901_gcs_loglikes[3]= calc_gcs_prob(group=0,/foreground,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=ngc1901_basic_filename+'_bestfit.sav',$
                                           /bestfitonly)
ENDIF

IF ~keyword_set(justonegroup) OR (keyword_set(justonegroup) AND group EQ 1) THEN BEGIN
    sirius_basic_filename='gcs_loglike_sirius'
    sirius_gcs_loglikes= dblarr(ntestspermg)
    sirius_gcs_loglikes[0]= calc_gcs_prob(group=1,/foreground,/fixalpha,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_fixalpha_sirius.sav',$
                                           grid=grid,agrid=0,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=sirius_basic_filename+'_fixalpha.sav',$
                                           sspspread=sspspread)
    sirius_gcs_loglikes[1]= calc_gcs_prob(group=1,/foreground,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_sirius.sav',$
                                           grid=grid,agrid=agrid,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=sirius_basic_filename+'.sav',$
                                           sspspread=sspspread)
    sirius_gcs_loglikes[2]= calc_gcs_prob(group=1,/foreground,/fixalpha,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=sirius_basic_filename+'_fixalpha_bestfit.sav',$
                                           /bestfitonly)
    sirius_gcs_loglikes[3]= calc_gcs_prob(group=1,/foreground,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=sirius_basic_filename+'_bestfit.sav',$
                                           /bestfitonly)
ENDIF

IF ~keyword_set(justonegroup) OR (keyword_set(justonegroup) AND group EQ 3) THEN BEGIN
    pleiades_basic_filename='gcs_loglike_pleiades'
    pleiades_gcs_loglikes= dblarr(ntestspermg)
    pleiades_gcs_loglikes[0]= calc_gcs_prob(group=3,/foreground,/fixalpha,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_fixalpha_pleiades.sav',$
                                           grid=grid,agrid=0,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=pleiades_basic_filename+'_fixalpha.sav',$
                                           sspspread=sspspread)
    pleiades_gcs_loglikes[1]= calc_gcs_prob(group=3,/foreground,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_pleiades.sav',$
                                           grid=grid,agrid=agrid,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=pleiades_basic_filename+'.sav',$
                                           sspspread=sspspread)
    pleiades_gcs_loglikes[2]= calc_gcs_prob(group=3,/foreground,/fixalpha,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=pleiades_basic_filename+'_fixalpha_bestfit.sav',$
                                           /bestfitonly)
    pleiades_gcs_loglikes[3]= calc_gcs_prob(group=3,/foreground,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=pleiades_basic_filename+'_bestfit.sav',$
                                           /bestfitonly)
ENDIF

IF ~keyword_set(justonegroup) OR (keyword_set(justonegroup) AND group EQ 6) THEN BEGIN
    hyades_basic_filename='gcs_loglike_hyades'
    hyades_gcs_loglikes= dblarr(ntestspermg)
    hyades_gcs_loglikes[0]= calc_gcs_prob(group=6,/foreground,/fixalpha,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_fixalpha_hyades.sav',$
                                           grid=grid,agrid=0,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=hyades_basic_filename+'_fixalpha.sav',$
                                           sspspread=sspspread)
    hyades_gcs_loglikes[1]= calc_gcs_prob(group=6,/foreground,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_hyades.sav',$
                                           grid=grid,agrid=agrid,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=hyades_basic_filename+'.sav',$
                                           sspspread=sspspread)
    hyades_gcs_loglikes[2]= calc_gcs_prob(group=6,/foreground,/fixalpha,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=hyades_basic_filename+'_fixalpha_bestfit.sav',$
                                           /bestfitonly)
    hyades_gcs_loglikes[3]= calc_gcs_prob(group=6,/foreground,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=hyades_basic_filename+'_bestfit.sav',$
                                           /bestfitonly)
ENDIF

IF ~keyword_set(justonegroup) OR (keyword_set(justonegroup) AND group EQ 7) THEN BEGIN
    hercules_basic_filename='gcs_loglike_hercules'
    hercules_gcs_loglikes= dblarr(ntestspermg)
    hercules_gcs_loglikes[0]= calc_gcs_prob(group=7,/foreground,/fixalpha,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_fixalpha_hercules.sav',$
                                           grid=grid,agrid=0,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=hercules_basic_filename+'_fixalpha.sav',$
                                           sspspread=sspspread)
    hercules_gcs_loglikes[1]= calc_gcs_prob(group=7,/foreground,$
                                           colorcut=colorcut,$
                                           logpostfilename='fit_ssp_fixpij_hercules.sav',$
                                           grid=grid,agrid=agrid,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           foreloglikefilename=hercules_basic_filename+'.sav',$
                                           sspspread=sspspread)
    hercules_gcs_loglikes[2]= calc_gcs_prob(group=7,/foreground,/fixalpha,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=hercules_basic_filename+'_fixalpha_bestfit.sav',$
                                           /bestfitonly)
    hercules_gcs_loglikes[3]= calc_gcs_prob(group=7,/foreground,$
                                           backloglikefilename=gcs_back_loglike_filename,$
                                           colorcut=colorcut,sspspread=sspspread,$
                                           foreloglikefilename=hercules_basic_filename+'_bestfit.sav',$
                                           /bestfitonly)
ENDIF

IF keyword_set(justonegroup) THEN return


;;Print all results
floatformat='(F8.1)'
splog, "Background GCS log probability"
splog, "------------------------------"
splog, gcs_back_loglike

splog, "NGC1901 GCS log probabilities"
splog, "------------------------------"
splog, "Marginalized, fixed alpha: "+strtrim(string(ngc1901_gcs_loglikes[0],format=floatformat),2)
splog, "Marginalized, free alpha: "+strtrim(string(ngc1901_gcs_loglikes[1],format=floatformat),2)
splog, "Best fit, fixed alpha: "+strtrim(string(ngc1901_gcs_loglikes[2],format=floatformat),2)
splog, "Best fit, free alpha: "+strtrim(string(ngc1901_gcs_loglikes[3],format=floatformat),2)

splog, "SIRIUS GCS log probabilities"
splog, "------------------------------"
splog, "Marginalized, fixed alpha: "+strtrim(string(sirius_gcs_loglikes[0],format=floatformat),2)
splog, "Marginalized, free alpha: "+strtrim(string(sirius_gcs_loglikes[1],format=floatformat),2)
splog, "Best fit, fixed alpha: "+strtrim(string(sirius_gcs_loglikes[2],format=floatformat),2)
splog, "Best fit, free alpha: "+strtrim(string(sirius_gcs_loglikes[3],format=floatformat),2)

splog, "PLEIADES GCS log probabilities"
splog, "------------------------------"
splog, "Marginalized, fixed alpha: "+strtrim(string(pleiades_gcs_loglikes[0],format=floatformat),2)
splog, "Marginalized, free alpha: "+strtrim(string(pleiades_gcs_loglikes[1],format=floatformat),2)
splog, "Best fit, fixed alpha: "+strtrim(string(pleiades_gcs_loglikes[2],format=floatformat),2)
splog, "Best fit, free alpha: "+strtrim(string(pleiades_gcs_loglikes[3],format=floatformat),2)

splog, "HYADES GCS log probabilities"
splog, "------------------------------"
splog, "Marginalized, fixed alpha: "+strtrim(string(hyades_gcs_loglikes[0],format=floatformat),2)
splog, "Marginalized, free alpha: "+strtrim(string(hyades_gcs_loglikes[1],format=floatformat),2)
splog, "Best fit, fixed alpha: "+strtrim(string(hyades_gcs_loglikes[2],format=floatformat),2)
splog, "Best fit, free alpha: "+strtrim(string(hyades_gcs_loglikes[3],format=floatformat),2)

splog, "HERCULES GCS log probabilities"
splog, "------------------------------"
splog, "Marginalized, fixed alpha: "+strtrim(string(hercules_gcs_loglikes[0],format=floatformat),2)
splog, "Marginalized, free alpha: "+strtrim(string(hercules_gcs_loglikes[1],format=floatformat),2)
splog, "Best fit, fixed alpha: "+strtrim(string(hercules_gcs_loglikes[2],format=floatformat),2)
splog, "Best fit, free alpha: "+strtrim(string(hercules_gcs_loglikes[3],format=floatformat),2)


END
