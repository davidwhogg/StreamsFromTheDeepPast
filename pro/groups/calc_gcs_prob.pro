;+
;   NAME:
;      calc_gcs_prob
;   PURPOSE:
;      calculate the gcs probability under the background model and
;      the various foreground models for the moving groups
;   INPUT:
;      group - which group? 20 for background
;      backdatafilename - filename that holds the background data
;      gcsdatafilename - filename that holds the GCS data
;      width - width parameter of the kernel
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      colorcut - if set, only use colors smaller than this value
;      logpostfilename - filename that holds p(age,Z(,alpha)) for this
;                        group
;      grid - number of grid points per Z, age dimension (actual
;             number of grid+1, best if 300 is divisible by grid)
;      agrid - grid parameter for alpha
;      age_range - range for the age parameter (log_age)
;      z_range - range for the metallicity
;      alpha_range - range for alpha
;      backloglikefilename - filename that either holds or will hold
;                            the background loglikelihoods
;      foreloglikefilename - filename that either holds or will hold
;                            the foreground loglikelihoods
;      solutionfilename - filename that holds the best fit mixture
;      sspspread - spread to include around the ssp prediction for the
;                  parallax
;   KEYWORDS:
;      foreground - if True, we are evaluating a foreground model
;      fixalpha - if True, we are evaluating a foreground model with
;                 fixed alpha
;      bestfitonly - only use  the best fit SSP values (these are set
;                    in stone in the code)
;   OUTPUT:
;      the logprobability of the GCS sample
;   HISTORY:
;      2009-11-17 - Written - Bovy (NYU)
;-
FUNCTION CALC_GCS_PROB, group=group, foreground=foreground, $
                        fixalpha=fixalpha,$
                        backdatafilename=backdatafilename,$
                        gcsdatafilename=gcsdatafilename, $
                        kernel=kernel, width=width, $
                        colorcut=colorcut, logpostfilename=logpostfilename,$
                        grid=grid, agrid=agrid, age_range=age_range, $
                        z_range=z_range, $
                        alpha_range=alpha_range, $
                        backloglikefilename=backloglikefilename, $
                        solutionfilename=solutionfilename, $
                        sspspread=sspspread,$
                        foreloglikefilename=foreloglikefilename,$
                        bestfitonly=bestfitonly
;;Restore background model data
IF ~keyword_set(backdatafilename) THEN restore, filename='../lsr2/hip2-aumer.sav' $
  ELSE restore, filename=backdatafilename
nhip= n_elements(hip.hip)
data= dblarr(nhip,2)
data[*,0]= hip.bvcolor
data[*,1]= hip.vmag-5.*alog10(100./hip.plx)

;;Restore gcs data
IF ~keyword_set(gcsdatafilename) THEN restore, filename='gcs.sav' ELSE restore, filename=gcsdatafilename
hip=gcs
nhip= n_elements(hip.hip)
IF keyword_set(colorcut) THEN BEGIN
    indx= where(gcs.bvcolor LT colorcut)
    hip= hip[indx]
    nhip= n_elements(indx)
ENDIF
ydata= transpose([[hip.vr]])
ycovar= hip.vrc
ycovar= reform(ycovar,1,1,nhip)
projection= hip.sm[*,0]
projection= reform(projection,3,1,nhip)

logl= 0.
IF keyword_set(foreground) THEN BEGIN
    IF ~keyword_set(grid) THEN grid= 300
    IF ~keyword_set(agrid) THEN IF keyword_set(fixalpha) THEN agrid= 0 ELSE agrid= 10
    IF ~keyword_set(age_range) THEN age_range= [6.6,10.2]
    IF ~keyword_set(z_range) THEN z_range= [0.0001,0.03]
    IF keyword_set(bestfitonly) THEN BEGIN
        grid=0
        IF keyword_set(fixalpha) THEN BEGIN
            IF group EQ 0 THEN BEGIN
                age_grid= [8.26]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 1 THEN BEGIN
                age_grid= [8.54]
                z_grid= [0.026]
            ENDIF ELSE IF group EQ 2 THEN BEGIN
                age_grid= [7.83]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 3 THEN BEGIN
                age_grid= [7.83]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 6 THEN BEGIN
                age_grid= [8.69]
                z_grid= [0.029]
            ENDIF ELSE IF group EQ 7 THEN BEGIN
                age_grid= [8.26]
                z_grid= [0.030]
            ENDIF
        ENDIF ELSE BEGIN
            IF group EQ 0 THEN BEGIN
                age_grid= [7.75]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 1 THEN BEGIN
                age_grid= [8.62]
                z_grid= [0.023]
            ENDIF ELSE IF group EQ 2 THEN BEGIN
                age_grid= [7.83]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 3 THEN BEGIN
                age_grid= [7.83]
                z_grid= [0.030]
            ENDIF ELSE IF group EQ 6 THEN BEGIN
                age_grid= [8.83]
                z_grid= [0.027]
            ENDIF ELSE IF group EQ 7 THEN BEGIN
                age_grid= [8.26]
                z_grid= [0.030]
            ENDIF
        ENDELSE
    ENDIF ELSE BEGIN
        age_step= (age_range[1]-age_range[0])/grid
        age_grid= dindgen(grid+1)*age_step+age_range[0]
        z_step= (z_range[1]-z_range[0])/grid
        z_grid= dindgen(grid+1)*z_step+z_range[0]
    ENDELSE
    IF ~keyword_set(fixalpha) AND ~keyword_set(alpha_range) THEN alpha_range=[0.00001,0.9999]
    IF keyword_set(fixalpha) THEN BEGIN
        agrid=0
        IF group EQ 0 THEN BEGIN
            name='NGC 1901'
            alpha_grid= [0.41]
        ENDIF ELSE IF group EQ 1 THEN BEGIN
            name='Sirius'
            alpha_grid= [0.53]
        ENDIF ELSE IF group EQ 2 THEN BEGIN
            name='Arcturus'
            alpha_grid= [0.37]
        ENDIF ELSE IF group EQ 3 THEN BEGIN
            name='Pleiades'
            alpha_grid= [0.57]
        ENDIF ELSE IF group EQ 6 THEN BEGIN
            name='Hyades'
            alpha_grid= [0.58]
        ENDIF ELSE IF group EQ 7 THEN BEGIN
            name='Hercules'
            alpha_grid= [0.83]
        ENDIF
    ENDIF ELSE IF keyword_set(bestfitonly) THEN BEGIN
        agrid=0
        IF group EQ 0 THEN BEGIN
            alpha_grid= [0.98]
        ENDIF ELSE IF group EQ 1 THEN BEGIN
            alpha_grid= [0.90]
        ENDIF ELSE IF group EQ 2 THEN BEGIN
            alpha_grid= [0.37];;Arcturus not run yet
        ENDIF ELSE IF group EQ 3 THEN BEGIN
            alpha_grid= [0.90]
        ENDIF ELSE IF group EQ 6 THEN BEGIN
            alpha_grid= [0.86]
        ENDIF ELSE IF group EQ 7 THEN BEGIN
            alpha_grid= [1.00]
        ENDIF
    ENDIF ELSE BEGIN
        alpha_step= (alpha_range[1]-alpha_range[0])/agrid
        alpha_grid= dindgen(agrid+1)*alpha_step+alpha_range[0]
    ENDELSE
    ;;Restore solution
    IF ~keyword_set(solutionfilename) THEN $
      solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
    restore, filename=solutionfilename
    ;;First compute all the probabilities
    assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
    IF group EQ 3 THEN BEGIN
        grouplogpost[group,*]= alog(exp(grouplogpost[group,*])+exp(grouplogpost[group+1,*]))
    ENDIF
    
    ;;Restore the posterior
    IF ~keyword_set(bestfitonly) THEN restore, filename=logpostfilename
    ;;Restore the background loglikes
    IF ~file_test(backloglikefilename) THEN BEGIN
        splog, "Error: backloglikefilename needs to be set"
        splog, "Run this procedure with the background first and save the results"
        splog, "Returning...."
        RETURN, -1
    ENDIF ELSE BEGIN
        restore, backloglikefilename
    ENDELSE
    ;;First calculate the likelihood of the foreground model
    IF file_test(foreloglikefilename) THEN BEGIN
        splog, "Restoring savefile "+foreloglikefilename
        restore, foreloglikefilename
    ENDIF ELSE BEGIN
        ii=0L
        jj=0L
        kk=0L
        gcslogpost= dblarr(agrid+1,grid+1,grid+1)
    ENDELSE
    WHILE ii LE grid DO BEGIN
        WHILE jj LE grid DO BEGIN
            WHILE kk LE agrid DO BEGIN
                IF keyword_set(fixalpha) THEN splog, 'Working on '+strtrim(string(jj+ii*(grid+1)+1),2)+'/'+strtrim(string(long((grid+1)^2.)),2) $
                ELSE splog, 'Working on '+strtrim(string(kk+(jj+ii*(grid+1))*(agrid+1)+1),2)+'/'+strtrim(string(long((agrid+1)*(grid+1)^2.)),2)
                gcslogpost[kk,ii,jj]= calc_logpost_fixedpij(group,alpha_grid[kk],age_grid[ii],z_grid[jj],sspspread=sspspread,back_loglike=logls,logpost=grouplogpost,hip=hip)
                gcslogpost[kk,ii,jj]+= total((1.-exp(grouplogpost[group,*]))*logls)
                print, gcslogpost[kk,ii,jj]
                kk+= 1L
                save, filename=foreloglikefilename, gcslogpost, ii, jj, kk
            ENDWHILE
            kk= 0L
            jj+= 1L
        ENDWHILE
        jj= 0L
        ii+= 1L
    ENDWHILE
    ;;Normalize the log posterior of age, ...
    IF ~keyword_set(bestfitonly) THEN logpost-= bovy_logsum(logpost)
    ;;Calculate the total loglikelihood
    IF keyword_set(bestfitonly) THEN logl= gcslogpost ELSE logl= bovy_logsum(gcslogpost+logpost)
ENDIF ELSE BEGIN
    ;;Background model
    IF file_test(backloglikefilename) THEN BEGIN
        splog, "savefile "+backloglikefilename+" exists"
        splog, "Restoring "+backloglikefilename
        restore, backloglikefilename
    ENDIF ELSE BEGIN
        logls= dblarr(nhip)
        FOR ii=0L, nhip-1 DO BEGIN
            logls[ii]= alog(kernel_phot_plx(hip[ii].plx,hip[ii].bvcolor,$
                                            hip[ii].vmag,$
                                            hip[ii].e_plx,data,$
                                            width=width,kernel=kernel,$
                                            silent=silent))
        ENDFOR
        IF ~keyword_set(backloglikefilename) THEN BEGIN
            splog, "You fool you should have set 'backloglikefilename' such that I could save the results for future reference"
        ENDIF ELSE BEGIN
            save, filename=backloglikefilename, logls
        ENDELSE
    ENDELSE
    logl= total(logls)
ENDELSE
RETURN, logl
END
