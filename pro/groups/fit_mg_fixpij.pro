;+
;   NAME:
;      fit_mg_fixpij
;   PURPOSE:
;      fit a moving group as an SSP, fixing the probabilities that
;      stars are part of the moving group
;   INPUT:
;      savefilename - filename for savefile
;      plotfilename - filename for plot
;      plotfilename_alpha - filename for plot of the marginal for
;                           alpha (if fixalpha is not set)
;      back_loglike_filename - filename in which the background log
;                              likelihoods can be found
;      datafilename - filename that holds the Hipparcos data
;      group - which moving group to fit?
;      sspspread - spread to include around the ssp prediction for the
;                  parallax
;      grid - number of grid points per Z, age dimension (actual
;             number of grid+1, best if 300 is divisible by grid)
;      agrid - grid parameter for alpha
;      age_range - range for the age parameter (log_age)
;      z_range - range for the metallicity
;      alpha_range - range for alpha
;      twodplotrange - plotrange for the age z density
;                      plot (number of e-cades to go down)
;      solutionfilename - filename that holds the best fit mixture
;   KEYWORDS:
;      fixalpha - if set, fix alpha at the default value (see source)
;      dontlogalpha - if set, plot the marginalized alpha
;                     probability, not the log of it
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-09 - Written - Bovy
;-
PRO FIT_MG_FIXPIJ, group, savefilename=savefilename, plotfilename=plotfilename,$
                   alphaplotfilename=alphaplotfilename, $
                   datafilename=datafilename,$
                   solutionfilename=solutionfilename, $
                   back_loglike_filename=back_loglike_filename, $
                   sspspread=sspspread, fixalpha=fixalpha, grid=grid, $
                   agrid=agrid, age_range=age_range, z_range=z_range, $
                   alpha_range=alpha_range, dontlogalpha=dontlogalpha
IF ~keyword_set(solutionfilename) THEN $
  solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
IF ~keyword_set(datafilename) THEN datafilename= '../lsr2/hip2-aumer.sav'
IF ~keyword_set(back_loglike_filename) THEN back_loglike_filename= 'back_loglike_all.sav'
IF ~keyword_set(grid) THEN grid= 300
IF ~keyword_set(agrid) THEN IF keyword_set(fixalpha) THEN agrid= 0 ELSE agrid= 10
IF ~keyword_set(age_range) THEN age_range= [6.6,10.2]
IF ~keyword_set(z_range) THEN z_range= [0.0001,0.03]
IF ~keyword_set(fixalpha) AND ~keyword_set(alpha_range) THEN alpha_range=[0.00001,0.9999]
age_step= (age_range[1]-age_range[0])/grid
age_grid= dindgen(grid+1)*age_step+age_range[0]
z_step= (z_range[1]-z_range[0])/grid
z_grid= dindgen(grid+1)*z_step+z_range[0]
IF keyword_set(fixalpha) THEN BEGIN
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
ENDIF ELSE BEGIN
    alpha_step= (alpha_range[1]-alpha_range[0])/agrid
    alpha_grid= dindgen(agrid+1)*alpha_step+alpha_range[0]
ENDELSE
;;First calculate the likelihood of the foreground model
IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    logpost= dblarr(agrid+1,grid+1,grid+1)
    ii= 0L
    jj= 0L
    kk= 0L
ENDELSE
IF file_test(back_loglike_filename) THEN BEGIN
    splog, 'Restoring background log likelihoods...'
    restore, back_loglike_filename
    back_loglike= logls
ENDIF ELSE BEGIN
    splog, 'Savefilename that holds the log likelihoods of the background model must be given'
    splog, 'Returning....'
    RETURN
ENDELSE
IF NOT (kk EQ agrid +1 AND ii EQ grid AND jj EQ grid) THEN BEGIN
    ;;Restore data
    restore, filename=datafilename
    nhip= n_elements(hip.hip)
    ydata= transpose([[hip.vl], [hip.vb]])
    ycovar= hip.vlvbc
    projection= hip.nsm
    ;;Restore solution
    restore, filename=solutionfilename
    ;;First compute all the probabilities
    assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
    ;;Sort the stars
    IF group EQ 3 THEN BEGIN
        grouplogpost[group,*]= alog(exp(grouplogpost[group,*])+exp(grouplogpost[group+1,*]))
    ENDIF
ENDIF
WHILE ii LE grid DO BEGIN
    WHILE jj LE grid DO BEGIN
        WHILE kk LE agrid DO BEGIN
            IF keyword_set(fixalpha) THEN splog, 'Working on '+strtrim(string(jj+ii*(grid+1)+1),2)+'/'+strtrim(string(long((grid+1)^2.)),2) $
              ELSE splog, 'Working on '+strtrim(string(kk+(jj+ii*(grid+1))*(agrid+1)+1),2)+'/'+strtrim(string(long((agrid+1)*(grid+1)^2.)),2)
            logpost[kk,ii,jj]= calc_logpost_fixedpij(group,alpha_grid[kk],age_grid[ii],z_grid[jj],sspspread=sspspread,back_loglike=back_loglike,logpost=grouplogpost,hip=hip)
            print, logpost[kk,ii,jj]
            kk+= 1L
            save, filename=savefilename, logpost, ii, jj, kk
        ENDWHILE
        kk= 0L
        jj+= 1L
    ENDWHILE
    jj= 0L
    ii+= 1L
ENDWHILE
;;Find and display best fit
IF keyword_set(fixalpha) THEN BEGIN
    maxloglike= max(logpost,maxindx)
    bestage= maxindx MOD (grid+1)
    bestz= maxindx/(grid+1)
    print, "Best fit age and metallicity are "+strtrim(string(10.^age_grid[bestage]/10.^6,format='(F7.2)'),2)+" and "+strtrim(string(z_grid[bestz],format='(F5.3)'),2)
ENDIF ELSE BEGIN
    maxloglike= max(logpost,maxindx)
    bestalpha= (maxindx MOD ((grid+1)*(agrid+1))) MOD (agrid+1)
    bestage= ((maxindx-bestalpha)/(agrid+1)) MOD (grid+1)
    bestz= ((maxindx-bestalpha)-bestage*(agrid+1))/(agrid+1)/(grid+1)
    print, "Best fit alpha, age and metallicity are "+strtrim(string(alpha_grid[bestalpha],format='(F4.2)'),2)+","+strtrim(string(10.^age_grid[bestage]/10.^6,format='(F7.2)'),2)+" and "+strtrim(string(z_grid[bestz],format='(F5.3)'),2)
ENDELSE

;;Now plot the 2D surface, first marginalize over alpha if necessary
IF group EQ 0 THEN BEGIN
    name='NGC 1901'
    fixed_alpha= 0.41
ENDIF ELSE IF group EQ 1 THEN BEGIN
    name='Sirius'
    fixed_alpha= 0.53
ENDIF ELSE IF group EQ 2 THEN BEGIN
    name='Arcturus'
    fixed_alpha= 0.37
ENDIF ELSE IF group EQ 3 THEN BEGIN
    name='Pleiades'
    fixed_alpha= 0.57
ENDIF ELSE IF group EQ 6 THEN BEGIN
    name='Hyades'
    fixed_alpha= 0.58
ENDIF ELSE IF group EQ 7 THEN BEGIN
    name='Hercules'
    fixed_alpha= 0.83
ENDIF
if ~keyword_set(fixalpha) THEN BEGIN
    plot_agez= dblarr(grid+1,grid+1)
    FOR ii=0L, grid DO FOR jj=0L, grid DO plot_agez[ii,jj]= bovy_logsum(logpost[*,ii,jj])
ENDIF ELSE BEGIN
    plot_agez= logpost[0,*,*]
ENDELSE
lognorm= max(plot_agez,max_indx)
;plot_agez-= lognorm
IF ~keyword_set(twodplotrange) THEN twodplotrange=[lognorm-5.,lognorm]
;;Calculate the right ranges
xrange=[age_grid[0]-age_step/2., age_grid[grid]+age_step/2.]
yrange=[z_grid[0]-z_step/2., z_grid[grid]+z_step/2.]
charsize=1.2
k_print, filename=plotfilename
bovy_density, plot_agez, xrange,yrange, $
  grid=[grid+1,grid+1], /keepbangX,title=name,charsize=charsize,charthick=2.,$
  range=range, xtitle='log Age [yr]', ytitle='Z',/flip, xcharsize=charsize, ycharsize=charsize
djs_oplot, [age_grid[max_indx MOD (grid+1)]], [z_grid[max_indx/(grid+1)]],color=djs_icolor('white'), psym=7, symsize=3, thick=4.
k_end_print

;;Plot marginal alpha distribution, if necessary
IF keyword_set(fixalpha) THEN RETURN
plot_alpha= dblarr(agrid+1)
FOR ii=0L, agrid DO BEGIN
    plot_alpha[ii]= bovy_logsum(logpost[ii,*,*])
ENDFOR
;;Normalize
lognorm= bovy_logsum(plot_alpha)
plot_alpha-= lognorm
IF keyword_set(dontlogalpha) THEN BEGIN
    ytitle= 'p(\alpha)'
    plot_alpha= exp(plot_alpha)
ENDIF ELSE BEGIN
    ytitle='log p(\alpha)'
ENDELSE
k_print, filename=alphaplotfilename
djs_plot, alpha_grid, plot_alpha, xtitle='\alpha', ytitle=ytitle, title=name,charsize=charsize,xcharsize=charsize,ycharsize=charsize,charthick=2
oplotbarx, fixed_alpha,thick=4., color=djs_icolor('gray')
k_end_print
END
