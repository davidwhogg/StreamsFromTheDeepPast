;+
;   NAME:
;      sellwood_model
;   PURPOSE:
;      investigate the space of the integrals of the motion of the
;      local disk to look for clues of Sellwood's recurrent spiral
;      structure
;   INPUT:
;      sample - 'GCS' or 'Hipparcos'
;      disk - 'Mestel' or implement something yourself
;      vc - circular velocity at the Sun
;      ro - distance from the Sun to the GC
;      plotfilename - filename for plot
;      vsunR - U velocity of the Sun
;      vsunphi - V velocity of the Sun wrt the lsr)
;      output - kind of output:
;                 1: GCS stars in (E,J) plane
;                 2: GCS stars with moving groups color-coded
;                 3: 6 panel plot, with all stars and the various
;                    moving groups
;      hardcut - hard cut to use when deciding whether a star is part
;                of a moving group
;      solutionfilename - filename with the velocity distribution
;   KEYWORDS:
;      colour - plot using color (if relevant) [there is a member
;               color in the velocity distribution file]
;   OUTPUT:
;   HISTORY:
;      2009-11-30 - Written - Bovy (NYU)
;-
FUNCTION potential_mestel, R,vc,ro
RETURN, vc^2.*alog(R/ro)
END
PRO SELLWOOD_MODEL, sample=sample, disk=disk, vc=vc, ro= ro, $
                    xsize=xsize, ysize=ysize, plotfilename=plotfilename, $
                    vsunR=vsunR, vsunphi=vsunphi, hardcut=hardcut, $
                    solutionfilename=solutionfilename,colour=colour, $
                    output=output
IF ~keyword_set(output) THEN output=1
IF ~keyword_set(vc) THEN vc= 235.
IF ~keyword_set(sample) THEN sample='GCS'
IF ~keyword_set(disk) THEN disk='mestel'
IF ~keyword_set(ro) THEN ro=8.2
IF ~keyword_set(vsunR) THEN vsunR= 10.1
IF ~keyword_set(vsunphi) THEN vsunphi= 4.0
IF ~keyword_set(hardcut) THEN hardcut= 0.4

;;Restore the sample
IF sample EQ 'GCS' THEN BEGIN
    gcssavefilename= 'phot_gcs.sav'
    splog, 'Restoring '+gcssavefilename
    restore, filename=gcssavefilename
    ngcs= n_elements(gcs.hip)
ENDIF ELSE BEGIN
    hipsavefilename= 'hip2-aumer-giants.sav'
    splog, 'Restoring '+hipsavefilename
    restore, filename=hipsavefilename
ENDELSE

IF output EQ 1 THEN BEGIN
    ;;Calculate E and J for all of the stars
    EJ= dblarr(ngcs,2)
    FOR ii=0L, ngcs-1 DO BEGIN
        l= gcs[ii].l/180.*!DPI
        dist= 1./gcs[ii].plx
        R= ro-dist*cos(l)
        phi= acos((ro^2.+R^2-dist^2.)/(2.*ro*R))
        vx= (gcs[ii].vr*sin(l)+gcs[ii].vl*cos(l))+vsunphi+vc
        vy= (-gcs[ii].vr*cos(l) + gcs[ii].vl*sin(l))-vsunR
        EJ[ii,0]= potential_mestel(R,vc,ro)+0.5*(vx^2.+vy^2.)
        EJ[ii,1]= abs(R*( vx * cos(phi)-$
                      vy * sin(phi)))
        ;;Convert into E_ran and L/L_circ
        EJ[ii,0]= (EJ[ii,0]-potential_mestel(EJ[ii,1]/vc,vc,ro)-0.5*vc^2.)/vc^2.
        EJ[ii,1]/= (ro*vc)
    ENDFOR
    charsize=.7
    k_print, filename=plotfilename, xsize=xsize, ysize=ysize
    djs_plot, EJ[*,1], EJ[*,0],xtitle='L/L_c(R_0)',$
      ytitle='(E-E_c) V_c^{-2}', title='Geneva-Copenhagen sample',$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize
    k_end_print
ENDIF ELSE IF output EQ 2 THEN BEGIN
    ;;Calculate E and J for all of the stars
    EJ= dblarr(ngcs,2)
    FOR ii=0L, ngcs-1 DO BEGIN
        l= gcs[ii].l/180.*!DPI
        dist= 1./gcs[ii].plx
        R= ro-dist*cos(l)
        phi= acos((ro^2.+R^2-dist^2.)/(2.*ro*R))
        vx= (gcs[ii].vr*sin(l)+gcs[ii].vl*cos(l))+vsunphi+vc
        vy= (-gcs[ii].vr*cos(l) + gcs[ii].vl*sin(l))-vsunR
        EJ[ii,0]= potential_mestel(R,vc,ro)+0.5*(vx^2.+vy^2.)
        EJ[ii,1]= abs(R*( vx * cos(phi)-$
                      vy * sin(phi)))
        ;;Convert into E_ran and L/L_circ
        EJ[ii,0]= (EJ[ii,0]-potential_mestel(EJ[ii,1]/vc,vc,ro)-0.5*vc^2.)/vc^2.
        EJ[ii,1]/= (ro*vc)
    ENDFOR
    ;;Then, for each moving group calculate the high and low metallicity
    ;;likelihood
    ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
    ycovar= gcs.vrvlvbc
    projection= gcs.sm
    ;;Restore solution
    IF ~keyword_set(solutionfilename) THEN $
      solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
    restore, filename=solutionfilename
    ;;First compute all the probabilities
    assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
    grouplogpost[3,*]= alog(exp(grouplogpost[3,*])+exp(grouplogpost[3+1,*]))
    ngc1901_members= where(exp(grouplogpost[0,*]) GT hardcut)
    sirius_members= where(exp(grouplogpost[1,*]) GT hardcut)
    pleiades_members= where(exp(grouplogpost[3,*]) GT hardcut)
    hyades_members= where(exp(grouplogpost[6,*]) GT hardcut)
    hercules_members= where(exp(grouplogpost[7,*]) GT hardcut)
    background_members = where( (exp(grouplogpost[0,*]) LE hardcut) AND $
      (exp(grouplogpost[1,*]) LE hardcut) AND $
      (exp(grouplogpost[3,*]) LE hardcut) AND $
      (exp(grouplogpost[6,*]) LE hardcut) AND $
      (exp(grouplogpost[7,*]) LE hardcut))
    print, n_elements(background_members)/double(ngcs) 
    ;;Plot
    charsize=.7
    k_print, filename=plotfilename, xsize=xsize, ysize=ysize
    ;;Background first
    djs_plot, EJ[background_members,1], EJ[background_members,0],$
      xtitle='L/L_c(R_0)',$
      ytitle='(E-E_c) V_c^{-2}', title='Geneva-Copenhagen sample',$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize
    ;;Moving groups next
    IF keyword_set(colour) THEN BEGIN
        djs_oplot, EJ[ngc1901_members,1], EJ[ngc1901_members,0],$
          color='red', psym=3
        legend, ['NGC 1901'], box=0.,charsize=charsize,textcolors=djs_icolor('red'),pos=[1.,.47]
        djs_oplot, EJ[sirius_members,1], EJ[sirius_members,0],$
          color='blue',psym=3
        legend, ['Sirius'], box=0.,charsize=charsize,textcolors=djs_icolor('blue'),pos=[1.,.44]
        djs_oplot, EJ[pleiades_members,1], EJ[pleiades_members,0],$
          color='green',psym=3
        legend, ['Pleiades'], box=0.,charsize=charsize,textcolors=djs_icolor('green'),pos=[1.,.41]
        djs_oplot, EJ[hyades_members,1], EJ[hyades_members,0],$
          color='yellow',psym=3
        legend, ['Hyades'], box=0.,charsize=charsize,textcolors=djs_icolor('yellow'),pos=[1.,.38]
        djs_oplot, EJ[hercules_members,1], EJ[hercules_members,0],$
          color='cyan',psym=3
        legend, ['Hercules'], box=0.,charsize=charsize,textcolors=djs_icolor('cyan'),pos=[1.,.35]
    ENDIF ELSE BEGIN
        djs_oplot, EJ[ngc1901_members,1], EJ[ngc1901_members,0],$
          psym=1,color='gray',symsize=.3
        legend, ['NGC 1901'], box=0.,charsize=charsize,pos=[1.,.47],psym=1
        djs_oplot, EJ[sirius_members,1], EJ[sirius_members,0],$
          psym=2,color='gray',symsize=.3
        legend, ['Sirius'], box=0.,charsize=charsize,pos=[1.,.44],psym=2
        djs_oplot, EJ[pleiades_members,1], EJ[pleiades_members,0],$
          psym=4,color='gray',symsize=.3
        legend, ['Pleiades'], box=0.,charsize=charsize,pos=[1.,.41],psym=4
        djs_oplot, EJ[hyades_members,1], EJ[hyades_members,0],$
          psym=5,color='gray',symsize=.3
        legend, ['Hyades'], box=0.,charsize=charsize,pos=[1.,.38],psym=5
        djs_oplot, EJ[hercules_members,1], EJ[hercules_members,0],$
          psym=7,color='gray',symsize=.3
        legend, ['Hercules'], box=0.,charsize=charsize,pos=[1.,.35],psym=7

    ENDELSE
    k_end_print
ENDIF ELSE IF output EQ 3 THEN BEGIN
    ;;Calculate E and J for all of the stars
    EJ= dblarr(ngcs,2)
    FOR ii=0L, ngcs-1 DO BEGIN
        l= gcs[ii].l/180.*!DPI
        dist= 1./gcs[ii].plx
        R= ro-dist*cos(l)
        phi= acos((ro^2.+R^2-dist^2.)/(2.*ro*R))
        vx= (gcs[ii].vr*sin(l)+gcs[ii].vl*cos(l))+vsunphi+vc
        vy= (-gcs[ii].vr*cos(l) + gcs[ii].vl*sin(l))-vsunR
        EJ[ii,0]= potential_mestel(R,vc,ro)+0.5*(vx^2.+vy^2.)
        EJ[ii,1]= abs(R*( vx * cos(phi)-$
                      vy * sin(phi)))
        ;;Convert into E_ran and L/L_circ
        EJ[ii,0]= (EJ[ii,0]-potential_mestel(EJ[ii,1]/vc,vc,ro)-0.5*vc^2.)/vc^2.
        EJ[ii,1]/= (ro*vc)
    ENDFOR
    ;;Then, for each moving group calculate the high and low metallicity
    ;;likelihood
    ydata= transpose([[gcs.vr],[gcs.vl], [gcs.vb]])
    ycovar= gcs.vrvlvbc
    projection= gcs.sm
    ;;Restore solution
    IF ~keyword_set(solutionfilename) THEN $
      solutionfilename= '../rv/ALMSv2_reconvergev2/fitv10_bestV_V4.0_sample6.sav'
    restore, filename=solutionfilename
    ;;First compute all the probabilities
    assign_clump_members,grouplogpost, ydata, ycovar, projection, mean, covar, amp
    grouplogpost[3,*]= alog(exp(grouplogpost[3,*])+exp(grouplogpost[3+1,*]))
    ngc1901_members= where(exp(grouplogpost[0,*]) GT hardcut)
    sirius_members= where(exp(grouplogpost[1,*]) GT hardcut)
    pleiades_members= where(exp(grouplogpost[3,*]) GT hardcut)
    hyades_members= where(exp(grouplogpost[6,*]) GT hardcut)
    hercules_members= where(exp(grouplogpost[7,*]) GT hardcut)
    background_members = where( (exp(grouplogpost[0,*]) LE hardcut) AND $
      (exp(grouplogpost[1,*]) LE hardcut) AND $
      (exp(grouplogpost[3,*]) LE hardcut) AND $
      (exp(grouplogpost[6,*]) LE hardcut) AND $
      (exp(grouplogpost[7,*]) LE hardcut))
    print, n_elements(background_members)/double(ngcs) 
    ;;Plot
    charsize=.6
    legendcharsize=.8
    k_print, filename=plotfilename, xsize=xsize, ysize=ysize
    ;;Set up positions
    xmargin=0.1
    ymargin= .1
    xwidth=.3
    ywidth=.3
    allposition= [xmargin,ymargin+2*ywidth,xmargin+xwidth,ywidth*3+ymargin]
    ngc1901position= [xmargin+xwidth,ymargin+2*ywidth,xmargin+2*xwidth,ywidth*3+ymargin]
    siriusposition= [xmargin,ymargin+ywidth,xmargin+xwidth,ywidth*2+ymargin]
    pleiadesposition= [xmargin+xwidth,ymargin+ywidth,xmargin+2*xwidth,ywidth*2+ymargin]
    hyadesposition= [xmargin,ymargin,xmargin+xwidth,ywidth+ymargin]
    herculesposition= [xmargin+xwidth,ymargin,xmargin+2*xwidth,ywidth+ymargin]
    ;;All first
    djs_plot, EJ[*,1], EJ[*,0],position=allposition,$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize, xtickname=strarr(30)+' ', $
      ytitle='(E-E_c) V_c^{-2}'
    legend, ['All stars'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    djs_plot, [0,1],[0,1],/nodata, $
      position=[xmargin,ymargin+2*ywidth,xmargin+2*xwidth,ywidth*3+ymargin], $
      title='Geneva-Copenhagen sample', /NOERASE,$
      xtickname=strarr(30)+' ', ytickname=strarr(30)+' ', xticks=1, yticks=1, charsize=charsize
    ;;Moving groups next
    djs_plot, EJ[ngc1901_members,1], EJ[ngc1901_members,0],$
      position=ngc1901position,/NOERASE,$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize, $
      xtickname=strarr(30)+' ', ytickname=strarr(30)+' '
    legend, ['NGC 1901'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    djs_plot, EJ[sirius_members,1], EJ[sirius_members,0],$
      position=siriusposition,/NOERASE,$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize, $
      xtickname=strarr(30)+' ', ytitle='(E-E_c) V_c^{-2}'
    legend, ['Sirius'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    djs_plot, EJ[pleiades_members,1], EJ[pleiades_members,0],$
      position=pleiadesposition,/NOERASE,$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize, $
      xtickname=strarr(30)+' ', ytickname=strarr(30)+' '
    legend, ['Pleiades'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    djs_plot, EJ[hyades_members,1], EJ[hyades_members,0],$
      position=hyadesposition,/NOERASE,$
      psym=3, yrange=[0.,0.5], xrange=[0.4,1.3],charsize=charsize, $
      ytitle='(E-E_c) V_c^{-2}', xtitle='L/L_c(R_0)'
    legend, ['Hyades'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    djs_plot, EJ[hercules_members,1], EJ[hercules_members,0],$
      position=herculesposition,/NOERASE,$
      psym=3, yrange=[0.,0.5], xrange=[0.4000001,1.3],charsize=charsize, $
      ytickname=strarr(30)+' ', xtitle='L/L_c(R_0)'
    legend, ['Hercules'], box=0.,charsize=legendcharsize,pos=[.8,.47]
    k_end_print
ENDIF
END
