;+
;   NAME:
;      tile_predictplx_iso
;   PURPOSE:
;      tile predictions for the plx, to make a pretty picture of the
;      foreground vs. background model
;   INPUT:
;      group - which group?
;      alpha - alpha parameter of the foreground model
;      age - age parameter of the foreground model (log Age)
;      z - metallicity parameter of the foreground model
;      fixed_age - best fit age while holding alpha fixed
;      fixed_z - best fit metallicity holding alpha fixed
;      width - width parameter of the kernel
;      mplot - array for multiplot (e.g. [0,4,4])
;      indices - array which holds the HIP numbers of a set of stars
;                to be plotted
;      plotfilename - filename for plot
;      sspspread - spread to add to the SSP prediction
;      seed    - seed for random
;      xsize - xsize for k_print
;      ysize - ysize for k_print
;      datafilename - filename with the data
;      nplx - discretization parameter for the function to be plotted
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      /random - pick random Hipparcos observations
;      /usegcs - use the GCS catalog, not the Hipparcos catalog
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-10 - Written - Bovy (NYU)
;-
PRO TILE_PREDICTPLX_ISO, group, alpha, age, z, fixed_age, fixed_z, $
                         sspspread=sspspread,$
                         mplot=mplot, indices=indices, $
                         xsize=xsize, ysize=ysize,$
                         plotfilename=plotfilename, seed=seed, $
                         random=random, $
                         datafilename=datafilename, width=width, $
                         kernel=kernel, $
                         nplx=nplx, usegcs=usegcs
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(mplot) THEN mplot= [0,3,2]
;;check inputs
IF ~mplot[1] EQ 3 THEN BEGIN
    splog, 'Error: mplot[1] has to be equal to 3'
    splog, 'Returning...'
    RETURN
ENDIF
ntile= mplot[1] * mplot[2]
IF ~keyword_set(random) THEN IF ~keyword_set(indices) THEN BEGIN
    splog, 'Either "/random" or "indices" must be set'
    splog, 'Returning...'
    RETURN
ENDIF ELSE BEGIN
    IF (mplot[2] NE n_elements(indices)) THEN BEGIN
        splog, 'Number of indices does not corresponds to the number of plots requested in mplot'
        splog, 'Returning...'
        RETURN
    ENDIF
ENDELSE
IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize
;!p.multi=mplot
;;If random is set, create a random set of "indices"
IF keyword_set(random) THEN BEGIN
    ;;Restore the Hipparcos data
    IF ~keyword_set(datafilename) THEN BEGIN
        IF keyword_set(usegcs) THEN BEGIN
            thisdatafilename='gcs.sav'
        ENDIF ELSE BEGIN
            thisdatafilename='../lsr2/hip2-aumer.sav'
        ENDELSE
    ENDIF ELSE BEGIN
        thisdatafilename=datafilename
    ENDELSE
    restore, filename=thisdatafilename
    IF keyword_set(usegcs) THEN hip= gcs
    indices= lonarr(ntile)
    FOR ii=0L, ntile-1 DO BEGIN
        indices[ii]= floor(n_elements(hip.hip) * randomu(seed))
    ENDFOR
ENDIF

;;Groups
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

;;Calculate some stuff related to the positioning of the subplots
xmargin=.1
ymargin=.1
xsep= 0.08
ysep= 0.05
width= (1-2*xmargin-(mplot[1]-1)*xsep)/mplot[1]
height= (1-.3-2*ymargin-(mplot[2]-1)*ysep)/mplot[2]
charsize=0.7
charthick=1.2
xcharsize=1
ycharsize=1
;;Run predictplx and predictplx_iso
indexindex= 0L
FOR ii=0L, mplot[2]-1 DO BEGIN
    ;;Left plot
    jj=0
    leftposition=[xmargin+jj*(width+xsep),$
                  1-ymargin-ii*(height+ysep)-height,$
                  xmargin+jj*(width+xsep)+width,$
                  1-ymargin-ii*(height+ysep)]
    IF ii EQ mplot[2]-1 THEN BEGIN
        predictplx, indices[indexindex], width=width, kernel=kernel, $
          datafilename=datafilename, nplx=nplx, /byindx,/hiplabel, /labelposition,$
          position=leftposition, /NOERASE, usegcs=usegcs,$
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE IF ii EQ 0 THEN BEGIN
        predictplx, indices[indexindex], width=width, kernel=kernel, $
          datafilename=datafilename, nplx=nplx, /byindx,/hiplabel,/labelposition,$
          /noxlabels, title='background', position=leftposition, usegcs=usegcs,$
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE BEGIN
        predictplx, indices[indexindex], width=width, kernel=kernel, $
          datafilename=datafilename, nplx=nplx, /byindx,/hiplabel,/labelposition,$
          /noxlabels, position=leftposition, /NOERASE, usegcs=usegcs,$
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDELSE
;;Middle plot
    jj=1
    middleposition=[xmargin+jj*(width+xsep),$
                  1-ymargin-ii*(height+ysep)-height,$
                  xmargin+jj*(width+xsep)+width,$
                  1-ymargin-ii*(height+ysep)]
    IF ii EQ mplot[2]-1 THEN BEGIN
        predictplx_iso, indices[indexindex], group, fixed_alpha, fixed_age, fixed_z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          position=middleposition, /NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE IF ii EQ 0 THEN BEGIN
        predictplx_iso, indices[indexindex], group, fixed_alpha, fixed_age, fixed_z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          /noxlabels, title=textoidl('foreground, fixed \alpha'), $
          position=middleposition, /NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE BEGIN
        predictplx_iso, indices[indexindex], group, fixed_alpha, fixed_age, fixed_z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          /noxlabels, position=middleposition, /NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDELSE
;;Right plot
    jj=2
    rightposition=[xmargin+jj*(width+xsep),$
                  1-ymargin-ii*(height+ysep)-height,$
                  xmargin+jj*(width+xsep)+width,$
                  1-ymargin-ii*(height+ysep)]
    IF ii EQ mplot[2]-1 THEN BEGIN
        predictplx_iso, indices[indexindex], group, alpha, age, z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          position=rightposition, /NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE IF ii EQ 0 THEN BEGIN
        predictplx_iso, indices[indexindex], group, alpha, age, z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          /noxlabels, title=textoidl('foreground, free \alpha'), $
          position=rightposition, /NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDIF ELSE BEGIN
        predictplx_iso, indices[indexindex], group, alpha, age, z, $
          width=width, kernel=kernel, usegcs=usegcs,$
          datafilename=datafilename, nplx=nplx, /byindx, $
          /noxlabels, position=rightposition,/NOERASE, $
          charsize=charsize, charthick=charthick, xcharsize=xcharsize, ycharsize=ycharsize
    ENDELSE
    indexindex += 1
ENDFOR
IF keyword_set(plotfilename) THEN k_end_print
END
