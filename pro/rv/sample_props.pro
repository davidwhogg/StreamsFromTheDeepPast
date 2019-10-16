;+
;   NAME:
;      sample_props
;   PURPOSE:
;      make plots of the sample properties of the hipparcos sample
;   CALLING SEQUENCE:
;   INPUT:
;      filename - if set, output ps plot to this file
;   KEYWORDS:
;      hip2    - If set, plot properties of new reduction
;      tyc     - If set, plot properties of combined Tycho and
;                Hipparcos sample
;   OUTPUT:
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO SAMPLE_PROPS, filename=filename, hip2=hip2, tyc=tyc


;;Restore
IF keyword_set(hip2) THEN BEGIN
    splog, 'Restoring hip2-dehnen.sav'
    restore, '/global/data/hipparcos/hipnewredux/hip2-dehnen.sav'
ENDIF ELSE IF keyword_set(tyc) THEN BEGIN
    splog, 'Restoring tyc-dehnen.sav'
    restore, 'tyc-dehnen.sav'
    hip= tyc
ENDIF ELSE BEGIN
    splog, 'Restoring hip-dehnen.sav'
    restore, '/global/data/hipparcos/hip-dehnen.sav'
ENDELSE


;;Make plots of RA, DEC, PLX, PMRA, PMDEC

IF keyword_set(filename) THEN k_print, filename=filename

xmargin= .1
ymargin= .1
xwidth= .16
ywidth= .16

charsize= .6


;;Make sample cuts
sample= where(hip.plx LE 40 AND abs(hip.pmra) LE 200 AND abs(hip.pmdec) LE 200)
ra= hip[sample].ra
dec= hip[sample].dec
plx= hip[sample].plx
pmra= hip[sample].pmra
pmdec= hip[sample].pmdec


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       HISTOGRAMS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;RA vs RA
hogg_plothist, ra, position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin], ytickname=REPLICATE(' ',30)
topxtitle, xrange= minmax(ra), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\alpha [deg]'), charsize=charsize, $
  position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin]

;;DEC vs DEC
hogg_plothist, dec, position=[xmargin+xwidth,1-ymargin-2*ywidth,xmargin+2*xwidth,1-ymargin-ywidth], /NOERASE

;;PLX vs PLX
hogg_plothist, plx, position=[xmargin+2*xwidth,1-ymargin-3*ywidth,xmargin+3*xwidth,1-ymargin-2*ywidth], /NOERASE

;;PMRA vs PMRA
hogg_plothist, pmra, position=[xmargin+3*xwidth,1-ymargin-4*ywidth,xmargin+4*xwidth,1-ymargin-3*ywidth], /NOERASE

;;PMDEC vs PMDEC
hogg_plothist, pmdec, position=[xmargin+4*xwidth,1-ymargin-5*ywidth,xmargin+5*xwidth,1-ymargin-4*ywidth], /NOERASE, $
  charsize=charsize, xtitle=textoidl('\mu_{\delta} [mas yr^{-1}]')



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       SCATTERPLOTS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;DEC VS RA
djs_plot, dec, ra, psym=3, position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
topxtitle, xrange= minmax(dec), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\delta [deg]'), charsize=charsize, $
  position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin]


;;PLX vs RA
djs_plot, plx, ra, psym=3, position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
topxtitle, xrange= minmax(plx), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\pi [mas]'), charsize=charsize, $
  position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin]

;;PMRA vs RA
djs_plot, pmra, ra, psym=3, position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
topxtitle, xrange= minmax(pmra), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\mu_{\alpha*} [mas yr^{-1}]'), charsize=charsize, $
  position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin]

;;PMDEC vs RA
djs_plot, pmdec, ra, psym=3, position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
rightytitle, xrange= minmax(pmdec), yrange=minmax(ra), xtickname=REPLICATE(' ',30), $
  ytitle=textoidl('\alpha [deg]'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin]
topxtitle, xrange= minmax(pmra), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\mu_{\delta} [mas yr^{-1}]'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin]

;;PLX vs DEC
djs_plot, plx, dec, psym=3, position=[xmargin+2*xwidth,1-ymargin-2*ywidth,xmargin+3*xwidth,1-ymargin-ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMRA vs DEC
djs_plot, pmra, dec, psym=3, position=[xmargin+3*xwidth,1-ymargin-2*ywidth,xmargin+4*xwidth,1-ymargin-ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs DEC
djs_plot, pmdec, dec, psym=3, position=[xmargin+4*xwidth,1-ymargin-2*ywidth,xmargin+5*xwidth,1-ymargin-ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
rightytitle, xrange= minmax(pmdec), yrange=minmax(dec), xtickname=REPLICATE(' ',30), $
  ytitle=textoidl('\delta [deg]'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-2*ywidth,xmargin+5*xwidth,1-ymargin-ywidth]


;;PMRA vs PLX
djs_plot, pmra, plx, psym=3, position=[xmargin+3*xwidth,1-ymargin-3*ywidth,xmargin+4*xwidth,1-ymargin-2*ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs PLX
djs_plot, pmdec, plx, psym=3, position=[xmargin+4*xwidth,1-ymargin-3*ywidth,xmargin+5*xwidth,1-ymargin-2*ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
rightytitle, xrange= minmax(pmdec), yrange=minmax(plx), xtickname=REPLICATE(' ',30), $
  ytitle=textoidl('\pi [mas]'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-3*ywidth,xmargin+5*xwidth,1-ymargin-2*ywidth]


;;PMDEC vs PMRA
djs_plot, pmdec, pmra, psym=3, position=[xmargin+4*xwidth,1-ymargin-4*ywidth,xmargin+5*xwidth,1-ymargin-3*ywidth], /NOERASE, $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)
rightytitle, xrange= minmax(pmdec), yrange=minmax(pmra), xtickname=REPLICATE(' ',30), $
  ytitle=textoidl('\mu_{\alpha*} [mas yr^{-1}]'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-4*ywidth,xmargin+5*xwidth,1-ymargin-3*ywidth]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       2D HISTOGRAMS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;DEC vs RA
bovy_plothist2d, hip.ra,hip.dec,/noerase,position=[xmargin,1-ymargin-2*ywidth,xmargin+1*xwidth,1-ymargin-1*ywidth], $
  charsize=charsize, ytitle=textoidl('\delta [deg]')

;;PLX vs RA
bovy_plothist2d, ra, plx,/noerase,position=[xmargin,1-ymargin-3*ywidth,xmargin+1*xwidth,1-ymargin-2*ywidth], $
  charsize=charsize, ytitle=textoidl('\pi [mas]')

;;PMRA vs RA
bovy_plothist2d, ra,pmra,/noerase,position=[xmargin,1-ymargin-4*ywidth,xmargin+1*xwidth,1-ymargin-3*ywidth], $
  charsize=charsize, ytitle=textoidl('\mu_{\alpha*} [mas yr^{-1}]')

;;PMDEC vs RA
bovy_plothist2d, ra,pmdec,/noerase,position=[xmargin,1-ymargin-5*ywidth,xmargin+1*xwidth,1-ymargin-4*ywidth], $
  charsize=charsize, ytitle=textoidl('\mu_{\delta} [mas yr^{-1}]'), xtitle=textoidl('\alpha [deg]')


;;PLX vs DEC
bovy_plothist2d, dec,plx,/noerase,position=[xmargin+xwidth,1-ymargin-3*ywidth,xmargin+2*xwidth,1-ymargin-2*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMRA vs DEC
bovy_plothist2d, dec,pmra,/noerase,position=[xmargin+xwidth,1-ymargin-4*ywidth,xmargin+2*xwidth,1-ymargin-3*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs DEC
bovy_plothist2d, dec,pmdec,/noerase,position=[xmargin+xwidth,1-ymargin-5*ywidth,xmargin+2*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\delta [deg]')


;;PMRA vs PLX
bovy_plothist2d, plx,pmra,/noerase,position=[xmargin+2*xwidth,1-ymargin-4*ywidth,xmargin+3*xwidth,1-ymargin-3*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs PLX
bovy_plothist2d, plx,pmdec,/noerase,position=[xmargin+2*xwidth,1-ymargin-5*ywidth,xmargin+3*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\pi [mas]')


;;PMDEC vs PMRA
bovy_plothist2d, pmra,pmdec,/noerase,position=[xmargin+3*xwidth,1-ymargin-5*ywidth,xmargin+4*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\mu_{\alpha*} [mas yr^{-1}]')



IF keyword_set(filename) THEN k_end_print

END
