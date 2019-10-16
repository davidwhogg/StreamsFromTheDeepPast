;+
;   NAME:
;      tile_kw_onepred.pro
;   PURPOSE:
;      plot one prediction of a gcs star as a function of k, w
;   CALLING SEQUENCE:
;   INPUT:
;      sample     - sample to use (default ALMS=6, see fitv for
;                   different samples)
;      basedir    - base directory for density estimate filenames
;                   (ends in '/')
;      filename   - filename for plot
;      xsize      - xsize for k_print
;      ysize      - ysize for k_print
;      seed       - seed for randomu
;      rrange     - v_r range
;      yrange     - prob range
;      ks         - ks to plot
;      ws         - ws to plot
;   KEYWORDS:
;      colour     - plot in color (in the density savefile there is a
;                   member color)
;   OUTPUT:
;   REVISION HISTORY:
;      2009-01-12 - Written Bovy (NYU)
;-
PRO TILE_KW_ONEPRED, sample=sample, filename=filename, xsize=xsize, $
                     ysize=ysize, seed=seed, colour=colour, rrange=rrange, $
                     yrange=yrange, basedir=basedir, ks=ks, ws=ws

d=3

IF ~keyword_set(sample) THEN sample= 6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(yrange) THEN yrange=[0,0.037]
IF ~keyword_set(rrange) THEN rrange=[-90,90]
IF ~keyword_set(ks) THEN ks= [3,5,7,10,13,15];ks= [3,5,10,15,20,25]
nks= n_elements(ks)
IF ~keyword_set(ws) THEN ws= [1.,4.,9.]
nws= n_elements(ws)


;;Restore gcs
genevafilename= '/global/data/rv/gcs/gcs.sav'
IF file_test(genevafilename) THEN BEGIN
    restore, genevafilename
ENDIF ELSE BEGIN
    splog, 'Error: GCS-savefile does not exist in ', $
      genevafilename
    splog, 'GCS-savefile must exist, returning...'
    RETURN
ENDELSE
ngcs= n_elements(gcs.hip)

star= floor(randomu(seed)*ngcs)

IF keyword_set(filename) THEN k_print, filename=filename, $
  xsize=xsize, ysize=ysize

xlab= textoidl('v_r [km s^{-1}]')
ylab= textoidl('P(v_r)')

xmargin= .1
ymargin= .1
ywidth= .13
xwidth= .13

charsize= .45
charthick=2

FOR kk= 0L, nks-1 DO FOR ww= 0L, nws-1 DO BEGIN
    ;;Restore the right file
    densityfilename= basedir+'fitv'+strtrim(string(ks[kk]),2)+$
      '_bestV_V'+strtrim(string(ws[ww],format='(f4.1)'),2)+$
      '_sample'+strtrim(string(sample),2)+$
      '.sav'
    restore, densityfilename
    IF ww EQ 0 THEN BEGIN
        noylabels=0
    ENDIF ELSE BEGIN
        noylabels=1
    ENDELSE
    IF kk EQ nks-1 THEN BEGIN
        noxlabels=0
    ENDIF ELSE BEGIN
        noxlabels=1
    ENDELSE
    predictrv, amp, mean, covar, gcs[star].sm,$
      vr=gcs[star].vr, $
      vlvb=[[gcs[star].vl],$
            [gcs[star].vb]], $
      cvlvb=gcs[star].vlvbc,$
      cvr=gcs[star].vrc,$
      loglike=loglike,altloglike=altloglike, $
      /plot, rrange=rrange,  yrange=yrange, color=colour, $
      position=[xmargin+ww*xwidth,ymargin+(nks-kk-1)*ywidth,$
      xmargin+(ww+1)*xwidth,ymargin+(nks-kk)*ywidth], /NOERASE, $
      charsize=charsize, noxlabels=noxlabels, $
      noylabels=noylabels, charthick=charthick
         
    
    IF kk EQ 0 THEN BEGIN
        legend, ["w="+strtrim(string(ws[ww],format='(f4.1)'),2)+textoidl(' km s^{-1}')],$
          pos=[xmargin+ww*xwidth+xwidth/2-.08,ymargin+(nks-kk)*ywidth+.03], $
          box=0, charsize=1.5*charsize,/normal
    ENDIF
    IF ww EQ nws-1 THEN BEGIN
        legend, ["K="+strtrim(string(ks[kk]),2)], $
          pos=[xmargin+(ww+1)*xwidth-.01,ymargin+(nks-kk)*ywidth-ywidth/2.+.02], $
          box=0, charsize=1.5*charsize,/normal
    ENDIF
ENDFOR

splog, "Success!"
splog, 'Plot was for HIP'+strtrim(string(gcs[star].hip),2)


IF keyword_set(filename) THEN k_end_print

END
