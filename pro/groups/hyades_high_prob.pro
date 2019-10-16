;+
;   NAME:
;      hyades_high_prob
;   PURPOSE:
;      to determine whether the Hyades cluster members from the
;      Perryman catalog have high probability of being in the Hyades
;      moving group
;   INPUT:
;      savefilename - filename for Hyades members data
;      plotfilename - name for plot
;      solutionfilename - name of the file that holds the
;                         decomposition of the velocity distribution
;   KEYWORDS:
;      colorcut     - cut everything out that has B-V > 1
;   OUTPUT:
;   REVISION HISTORY:
;      2009-12-01 - Written - Bovy (NYU)
;-
PRO HYADES_HIGH_PROB, savefilename=savefilename, $
                      plotfilename=plotfilename, $
                      colorcut=colorcut, $
                      solutionfilename=solutionfilename
;;Restore Hyades members
IF ~file_test(savefilename) THEN BEGIN
    hyades= read_hyades(savefilename=savefilename,/single,/tenpc)
ENDIF ELSE BEGIN
    restore, savefilename
ENDELSE
;;select (possible) MS members
IF keyword_set(colorcut) THEN usethese= where(hyades.vmag-5.*alog10(100./hyades.plx) GE 0.5 $
                                              AND hyades.bvcolor LT 1.) $
  ELSE usethese= where(hyades.vmag-5.*alog10(100./hyades.plx) GE 0.5)
nhyades= n_elements(usethese)
splog, 'Using '+strtrim(string(nhyades),2)+' Hyades members'
hyades=hyades[usethese]

;;Restore velocity distribution
IF ~keyword_set(solutionfilename) THEN $
  solutionfilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
ydata= transpose([[hyades.vl], [hyades.vb]])
ycovar= hyades.vlvbc
projection= hyades.nsm
;;Restore solution
restore, filename=solutionfilename

;;First compute all the probabilities
assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
logpost[3,*]= alog(exp(logpost[3,*])+exp(logpost[3+1,*]))

;;Histogram the Hyades membership probabilities
k_print, filename=plotfilename
hogg_plothist, exp(logpost[6,*]),xtitle=textoidl('p_{Hyades}'),$
  ytitle='N',$
  title='Hyades open cluster members = Hyades moving group members?',$
  /totalweight
k_end_print
high= where(exp(logpost[6,*]) GT 0.4)
splog, strtrim(string(n_elements(high)),2)+' Hyades open cluster members have high probability of being moving group members'
veryhigh= where(exp(logpost[6,*]) GT 0.7)
splog, strtrim(string(n_elements(veryhigh)),2)+' Hyades open cluster members have high probability of being moving group members'

END
