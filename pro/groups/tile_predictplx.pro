;+
;   NAME:
;      tile_predictplx
;   PURPOSE:
;      tile predictions for the plx, to make a pretty picture of the
;      background model
;   INPUT:
;      width - width parameter of the kernel
;      mplot - array for multiplot (e.g. [0,4,4])
;      indices - array which holds the HIP numbers of a set of stars
;                to be plotted
;      plotfilename - filename for plot
;      seed    - seed for random
;      xsize - xsize for k_print
;      ysize - ysize for k_print
;      datafilename - filename with the data
;      nplx - discretization parameter for the function to be plotted
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;      /random - pick random Hipparcos observations
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-04 - Written - Bovy (NYU)
;-
PRO TILE_PREDICTPLX, mplot=mplot, indices=indices, xsize=xsize, ysize=ysize,$
                     plotfilename=plotfilename, seed=seed, random=random, $
                     datafilename=datafilename, width=width, kernel=kernel, $
                     nplx=nplx
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN datafilename='../lsr2/hip2-aumer.sav'
restore, filename=datafilename
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(mplot) THEN mplot= [0,4,4]
;;check inputs
ntile= mplot[1] * mplot[2]
IF ~keyword_set(random) THEN IF ~keyword_set(indices) THEN BEGIN
    splog, 'Either "/random" or "indices" must be set'
    splog, 'Returning...'
    RETURN
ENDIF ELSE BEGIN
    IF (ntile NE n_elements(indices)) THEN BEGIN
        splog, 'Number of indices does not corresponds to the number of plots requested in mplot'
        splog, 'Returning...'
        RETURN
    ENDIF
ENDELSE
IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize
!p.multi=mplot
;;If random is set, create a random set of "indices"
IF keyword_set(random) THEN BEGIN
    indices= lonarr(ntile)
    FOR ii=0L, ntile-1 DO BEGIN
        indices[ii]= floor(n_elements(hip.hip) * randomu(seed))
    ENDFOR
ENDIF
;;Run predictplx
indexindex= 0L
FOR ii=0L, mplot[2]-1 DO BEGIN
    FOR jj=0L, mplot[1]-1 DO BEGIN
        IF jj EQ 0 THEN BEGIN
            IF ii EQ mplot[2]-1 THEN BEGIN
                predictplx, indices[indexindex], width=width, kernel=kernel, $
                  datafilename=datafilename, nplx=nplx, /byindx,/hiplabel
            ENDIF ELSE BEGIN
                predictplx, indices[indexindex], width=width, kernel=kernel, $
                  datafilename=datafilename, nplx=nplx, /byindx,/hiplabel, $
                  /noxlabels
            ENDELSE
        ENDIF ELSE BEGIN
            IF ii EQ mplot[2]-1 THEN BEGIN
                predictplx, indices[indexindex], width=width, kernel=kernel, $
                  datafilename=datafilename, nplx=nplx, /byindx,/hiplabel, $
                  /noylabels
            ENDIF ELSE BEGIN
                predictplx, indices[indexindex], width=width, kernel=kernel, $
                  datafilename=datafilename, nplx=nplx, /byindx,/hiplabel, $
                  /noxlabels, /noylabels
            ENDELSE
        ENDELSE
        indexindex += 1
    ENDFOR
ENDFOR
IF keyword_set(plotfilename) THEN k_end_print
END
