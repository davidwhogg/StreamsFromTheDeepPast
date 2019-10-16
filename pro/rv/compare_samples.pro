;+
;   NAME:
;      compare_samples
;   PURPOSE:
;      compare the original Hipparcos sample with the new reduction
;      and the sample that combines Hipparcos with Tycho proper motions
;   CALLING SEQUENCE:
;   INPUT:
;      filename   - If set, output the plot to this file
;   KEYWORDS:
;      tyc    - If set, compare Hipparcos with Tycho
;   OUTPUT;
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO COMPARE_SAMPLES, tyc=tyc, filename=filename

IF keyword_set(filename) THEN k_print, filename=filename, $
  xsize=6.9

xmargin= .1
ymargin= .1
xwidth= .16
ywidth= .16

charsize= .4

restore, '~/streams/pro/rv/hip-dehnen.sav'
hip_tmp= hip

IF keyword_set(tyc) THEN BEGIN
    splog, 'Comparing Hipparcos with Tycho...'
    restore, '~/streams/pro/rv/tyc-dehnen.sav',/RELAXED_STRUCTURE_ASSIGNMENT
    hip2=tyc
    hip=hip_tmp
ENDIF ELSE BEGIN
    splog, 'Comparing old Hipparcos with the new reduction'
    ;;Restore both
    restore, '~/streams/pro/rv/hip2-dehnen.sav'
    hip2= hip
    hip= hip_tmp
ENDELSE

IF keyword_set(tyc) THEN xrange1=[0,19] ELSE xrange1=[0,1.45]
IF keyword_set(tyc) THEN xrange2=[0,4.4] ELSE xrange2=[0,1.45]
IF keyword_set(tyc) THEN xtickinterval1= 5 ELSE xtickinterval1=.5
IF keyword_set(tyc) THEN xtickinterval2= 1.5 ELSE xtickinterval2=.5

;    sample= where(hip.e_ra LE 1.3 AND hip.e_dec LE 1.3 AND hip.e_plx LE 2.4 AND hip.e_pmra LE 1.6 AND hip.e_pmdec LE 1.6)
;    e_ra= hip[sample].e_ra
;    e_dec= hip[sample].e_dec
;    e_plx= hip[sample].e_plx
;    e_pmra= hip[sample].e_pmra
;    e_pmdec= hip[sample].e_pmdec

e_ra= hip2.e_ra/hip.e_ra
e_dec= hip2.e_dec/hip.e_dec
e_plx= hip2.e_plx/hip.e_plx
e_pmra= hip2.e_pmra/hip.e_pmra
e_pmdec= hip2.e_pmdec/hip.e_pmdec

;;e_RA vs e_RA
IF keyword_set(tyc) THEN xtitle= textoidl('\Delta^{tyc}\alpha/\Delta^{hip}\alpha') $
  ELSE xtitle= textoidl('\Delta^{new}\alpha/\Delta^{old}\alpha')
hogg_plothist, e_ra, position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin], ytickname=REPLICATE(' ',30), $
  xtitle=xtitle, charsize=charsize, xtickinterval=xtickinterval1, xrange=xrange1

;;e_DEC vs e_DEC
IF keyword_set(tyc) THEN xtitle= textoidl('\Delta^{tyc}\delta/\Delta^{hip}\delta') $
  ELSE xtitle= textoidl('\Delta^{new}\delta/\Delta^{old}\delta')
hogg_plothist, e_dec, position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtitle=xtitle, charsize=charsize, xtickinterval=xtickinterval1, xrange=xrange1

;;e_PLX vs e_PLX
IF keyword_set(tyc) THEN xtitle= textoidl('\Delta^{tyc}\pi/\Delta^{hip}\pi') $
  ELSE xtitle= textoidl('\Delta^{new}\pi/\Delta^{old}\pi')
hogg_plothist, e_plx, position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtitle=xtitle, charsize=charsize, xtickinterval=xtickinterval2, xrange=xrange2

;;e_PMRA vs e_PMRA
IF keyword_set(tyc) THEN xtitle= textoidl('\Delta^{tyc}\mu_{\alpha}/\Delta^{hip}\mu_{\alpha}') $
  ELSE xtitle= textoidl('\Delta^{new}\mu_{\alpha}/\Delta^{old}\mu_{\alpha}')
hogg_plothist, e_pmra, position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtitle= xtitle, charsize=charsize, xtickinterval=xtickinterval2, xrange=xrange2

;;e_PMDEC vs e_PMDEC
IF keyword_set(tyc) THEN xtitle= textoidl('\Delta^{tyc}\mu_{\delta}/\Delta^{hip}\mu_{\delta}') $
  ELSE xtitle= textoidl('\Delta^{new}\mu_{\delta}/\Delta^{old}\mu_{\delta}')
hogg_plothist, e_pmdec, position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin], /NOERASE, $
  charsize=charsize, xtitle=xtitle, ytickname=REPLICATE(' ',30), xtickinterval=xtickinterval2, xrange=xrange2

IF keyword_set(filename) THEN k_end_print

END
