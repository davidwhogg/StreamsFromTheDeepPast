;+
;   NAME:
;      hist_loglike_back
;   PURPOSE:
;      make a histogram of the log likelihoods in the background
;      model, both for all stars as well as for the groups
;   INPUT:
;      width - width parameter of the kernel
;      plotfilename - filename to plot to
;      savefilename - filename to save the loglikes in
;      xsize   - x size for k_print (can be left unset)
;      ysize   - ysize for k_print
;      datafilename - filename with the data
;      group - which group?, 20 for all stars
;      denssavefilename - filename that holds the velocity distribution
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;   KEYWORDS:
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-04 - Written - Bovy (NYU)
;-
PRO HIST_LOGLIKE_BACK, width=width, kernel=kernel, plotfilename=plotfilename, $
                       savefilename=savefilename, xsize=xsize, ysize=ysize, $
                       datafilename=datafilename, group=group, $
                       denssavefilename=denssavefilename
IF ~keyword_set(cut) THEN cut=0.
IF ~keyword_set(denssavefilename) THEN $
  denssavefilename='../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
charsize=.7
;;Restore data
restore, filename='../lsr2/hip2-aumer.sav'
nhip= n_elements(hip.hip)
data= dblarr(nhip,2)
data[*,0]= hip.bvcolor
data[*,1]= hip.vmag-5.*alog10(100./hip.plx)
IF keyword_set(savefilename) and file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    logls= dblarr(nhip)
    FOR ii=0L, nhip-1 DO BEGIN
        logls[ii]= alog(kernel_phot_plx(hip[ii].plx,hip[ii].bvcolor,$
                                        hip[ii].vmag,$
                                        hip[ii].e_plx,data,$
                                        width=width,kernel=kernel,$
                                        silent=silent))
    ENDFOR
    IF keyword_set(savefilename) THEN save, filename=savefilename, logls
ENDELSE
IF group EQ 20 THEN BEGIN
    weight = make_array(nhip,/float,Value=1/float(nhip))
    name='All stars'
ENDIF ELSE BEGIN
    ydata= transpose([[hip.vl], [hip.vb]])
    ycovar= hip.vlvbc
    projection= hip.nsm
    ;;Restore solution
    restore, filename=denssavefilename
    ;;First compute all the probabilities
    assign_clump_members,logpost, ydata, ycovar, projection, mean, covar, amp
    ;;Sort the stars
    IF group EQ 3 THEN BEGIN
        logpost[group,*]= alog(exp(logpost[group,*])+exp(logpost[group+1,*]))
    ENDIF
    basefilename='back_loglike_'
    IF group EQ 0 THEN BEGIN
        name='NGC 1901'
        plotfilename=basefilename+'ngc1901'
    ENDIF ELSE IF group EQ 1 THEN BEGIN
        name='Sirius'
        plotfilename=basefilename+'sirius'
    ENDIF ELSE IF group EQ 2 THEN BEGIN
        name='Arcturus'
        plotfilename=basefilename+'arcturus'
    ENDIF ELSE IF group EQ 3 THEN BEGIN
        name='Pleiades'
        plotfilename=basefilename+'pleiades'
    ENDIF ELSE IF group EQ 6 THEN BEGIN
        name='Hyades'
        plotfilename=basefilename+'hyades'
    ENDIF ELSE IF group EQ 7 THEN BEGIN
        name='Hercules'
        plotfilename=basefilename+'hercules'
    ENDIF
    plotfilename+= '.ps'
    weight= exp(logpost[group,*])
    weight= weight/total(weight)
ENDELSE
;;Plot
k_print, filename=plotfilename, xsize=xsize, ysize=ysize
hogg_plothist, logls, weight=weight, xtitle=textoidl('p(\pi!10#!Xbackground model)'), $
  ytitle=textoidl('fraction'), charsize=charsize, title=name
k_end_print
END
