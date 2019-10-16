;+
;   NAME:
;      sample_err_props
;   PURPOSE:
;      make plots of the error properties of the hipparcos sample
;   CALLING SEQUENCE:
;   INPUT:
;      hip2     - if set, output error properties of new reduction
;      tyc      - if set, output error properties of Tycho sample
;      filename - if set, output ps plot to this file
;   KEYWORDS:
;   OUTPUT:
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO SAMPLE_ERR_PROPS, filename=filename, hip2=hip2, tyc=tyc


;;Make plots of e_RA, e_DEC, e_PLX, e_PMRA, e_PMDEC

IF keyword_set(filename) THEN k_print, filename=filename

xmargin= .1
ymargin= .1
xwidth= .16
ywidth= .16

charsize= .4

IF keyword_set(hip2) THEN BEGIN
    splog, 'Analysing error properties of the new reduction of the Hipparcos data...'
    restore, '~/streams/pro/rv/hip2-dehnen.sav'

;;Make sample cuts
    sample= where(hip.e_ra LE 1.3 AND hip.e_dec LE 1.3 AND hip.e_plx LE 2.4 AND hip.e_pmra LE 1.6 AND hip.e_pmdec LE 1.6)
    e_ra= hip[sample].e_ra
    e_dec= hip[sample].e_dec
    e_plx= hip[sample].e_plx
    e_pmra= hip[sample].e_pmra
    e_pmdec= hip[sample].e_pmdec

;    e_ra= hip.e_ra
;    e_dec= hip.e_dec
;    e_plx= hip.e_plx
;    e_pmra= hip.e_pmra
;    e_pmdec= hip.e_pmdec

;;e_RA vs e_RA
    hogg_plothist, e_ra, position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin], ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\alpha [mas]'), charsize=charsize, xtickinterval=.5

;;e_DEC vs e_DEC
    hogg_plothist, e_dec, position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\delta [mas]'), charsize=charsize, xtickinterval=.5

;;e_PLX vs e_PLX
hogg_plothist, e_plx, position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\pi [mas]'), charsize=charsize, xtickinterval=.5

;;e_PMRA vs e_PMRA
hogg_plothist, e_pmra, position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\mu_{\alpha} [mas yr^{-1}]'), charsize=charsize, xtickinterval=.5

;;e_PMDEC vs e_PMDEC
hogg_plothist, e_pmdec, position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin], /NOERASE, $
  charsize=charsize, xtitle=textoidl('\Delta\mu_{\delta} [mas yr^{-1}]'), ytickname=REPLICATE(' ',30), xtickinterval=.5
    
ENDIF ELSE IF keyword_set(tyc) THEN BEGIN
    splog, 'Analysing the error properties of the Tycho + Hipparcos data..'
    restore, 'tyc-dehnen.sav'

;;Make sample cuts
    sample= where(tyc.e_ra LE 19 AND tyc.e_dec LE 19 AND tyc.e_plx LE 2.4 AND tyc.e_pmra LE 2.4 AND tyc.e_pmdec LE 2.4)
    splog, n_elements(sample), n_elements(tyc.hip)
    e_ra= tyc[sample].e_ra
    e_dec= tyc[sample].e_dec
    e_plx= tyc[sample].e_plx
    e_pmra= tyc[sample].e_pmra
    e_pmdec= tyc[sample].e_pmdec

;    e_ra= tyc.e_ra
;    e_dec= tyc.e_dec
;    e_plx= tyc.e_plx
;    e_pmra= tyc.e_pmra
;    e_pmdec= tyc.e_pmdec

;;e_RA vs e_RA
    hogg_plothist, e_ra, position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin], ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\alpha [mas]'), charsize=charsize, xtickinterval=5

;;e_DEC vs e_DEC
    hogg_plothist, e_dec, position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\delta [mas]'), charsize=charsize, xtickinterval=5

;;e_PLX vs e_PLX
hogg_plothist, e_plx, position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\pi [mas]'), charsize=charsize, xtickinterval=.5

;;e_PMRA vs e_PMRA
hogg_plothist, e_pmra, position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin], /NOERASE, ytickname=REPLICATE(' ',30), $
      xtitle=textoidl('\Delta\mu_{\alpha} [mas yr^{-1}]'), charsize=charsize, xtickinterval=.5

;;e_PMDEC vs e_PMDEC
hogg_plothist, e_pmdec, position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin], /NOERASE, $
  charsize=charsize, xtitle=textoidl('\Delta\mu_{\delta} [mas yr^{-1}]'), ytickname=REPLICATE(' ',30), xtickinterval=.5

ENDIF ELSE BEGIN
splog, 'Restoring hip-dehnen.sav'
restore, '~/streams/pro/rv/hip-dehnen.sav'

corr_range=[-.95,.95]

;;Make sample cuts
sample= where(hip.e_ra LE 1.3 AND hip.e_dec LE 1.3 AND hip.e_plx LE 2.4 AND hip.e_pmra LE 1.6 AND hip.e_pmdec LE 1.6)
e_ra= hip[sample].e_ra
e_dec= hip[sample].e_dec
e_plx= hip[sample].e_plx
e_pmra= hip[sample].e_pmra
e_pmdec= hip[sample].e_pmdec

c_decra= hip.c_decra
c_plxra= hip.c_plxra
c_plxdec= hip.c_plxdec
c_pmrara= hip.c_pmrara
c_pmradec= hip.c_pmradec
c_pmraplx= hip.c_pmraplx
c_pmdecra= hip.c_pmdecra
c_pmdecdec= hip.c_pmdecdec
c_pmdecplx= hip.c_pmdecplx
c_pmdecpmra= hip.c_pmdecpmra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       HISTOGRAMS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;e_RA vs e_RA
hogg_plothist, e_ra, position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin], ytickname=REPLICATE(' ',30), /NOERASE, xtickname=REPLICATE(' ',30)
topxtitle, xrange= minmax(e_ra), ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\Delta\alpha [mas]'), charsize=charsize, $
  position=[xmargin,1-ymargin-ywidth,xmargin+xwidth,1-ymargin]

;;e_DEC vs e_DEC
hogg_plothist, e_dec, position=[xmargin+xwidth,1-ymargin-2*ywidth,xmargin+2*xwidth,1-ymargin-ywidth], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtickname=REPLICATE(' ',30)

;;e_PLX vs e_PLX
hogg_plothist, e_plx, position=[xmargin+2*xwidth,1-ymargin-3*ywidth,xmargin+3*xwidth,1-ymargin-2*ywidth], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtickname=REPLICATE(' ',30)

;;e_PMRA vs e_PMRA
hogg_plothist, e_pmra, position=[xmargin+3*xwidth,1-ymargin-4*ywidth,xmargin+4*xwidth,1-ymargin-3*ywidth], /NOERASE, ytickname=REPLICATE(' ',30), $
  xtickname=REPLICATE(' ',30)

;;e_PMDEC vs e_PMDEC
hogg_plothist, e_pmdec, position=[xmargin+4*xwidth,1-ymargin-5*ywidth,xmargin+5*xwidth,1-ymargin-4*ywidth], /NOERASE, $
  charsize=charsize, xtitle=textoidl('\Delta\mu_{\delta} [mas yr^{-1}]'), ytickname=REPLICATE(' ',30), xtickinterval=.5



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       CORRELATION HISTOGRAMS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;corr DEC VS RA
hogg_plothist, c_decra, position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range
topxtitle, xrange= corr_range, ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\rho_{\delta\alpha}'), charsize=charsize, $
  position=[xmargin+xwidth,1-ymargin-ywidth,xmargin+2*xwidth,1-ymargin]

;;corr PLX vs RA
hogg_plothist, c_plxra, position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range
topxtitle, xrange= corr_range, ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\rho_{\pi\alpha}'), charsize=charsize, $
  position=[xmargin+2*xwidth,1-ymargin-ywidth,xmargin+3*xwidth,1-ymargin]

;;PMRA vs RA
hogg_plothist, c_pmrara, position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range
topxtitle, xrange= corr_range, ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\rho_{\alpha\mu_{\alpha}}'), charsize=charsize, $
  position=[xmargin+3*xwidth,1-ymargin-ywidth,xmargin+4*xwidth,1-ymargin]

;;PMDEC vs RA
hogg_plothist, c_pmdecra, position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range
topxtitle, xrange= corr_range, ytickname=REPLICATE(' ',30), $
  xtitle=textoidl('\rho_{\alpha\mu_{\delta}}'), charsize=charsize, $
  position=[xmargin+4*xwidth,1-ymargin-ywidth,xmargin+5*xwidth,1-ymargin]

;;PLX vs DEC
hogg_plothist, c_plxdec, position=[xmargin+2*xwidth,1-ymargin-2*ywidth,xmargin+3*xwidth,1-ymargin-ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range

;;PMRA vs DEC
hogg_plothist, c_pmradec, position=[xmargin+3*xwidth,1-ymargin-2*ywidth,xmargin+4*xwidth,1-ymargin-ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range

;;PMDEC vs DEC
hogg_plothist, c_pmdecdec, position=[xmargin+4*xwidth,1-ymargin-2*ywidth,xmargin+5*xwidth,1-ymargin-ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range


;;PMRA vs PLX
hogg_plothist, c_pmraplx, position=[xmargin+3*xwidth,1-ymargin-3*ywidth,xmargin+4*xwidth,1-ymargin-2*ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range

;;PMDEC vs PLX
hogg_plothist, c_pmdecplx, position=[xmargin+4*xwidth,1-ymargin-3*ywidth,xmargin+5*xwidth,1-ymargin-2*ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range


;;PMDEC vs PMRA
hogg_plothist, c_pmdecpmra, position=[xmargin+4*xwidth,1-ymargin-4*ywidth,xmargin+5*xwidth,1-ymargin-3*ywidth], ytickname=REPLICATE(' ',30), /NOERASE, $
  xtickname=REPLICATE(' ',30), xrange=corr_range



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;       2D HISTOGRAMS
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;e_DEC vs e_RA
bovy_plothist2d, e_ra, e_dec,/noerase,position=[xmargin,1-ymargin-2*ywidth,xmargin+1*xwidth,1-ymargin-1*ywidth], $
  charsize=charsize, ytitle=textoidl('\Delta\delta [mas]')


;;PLX vs RA
bovy_plothist2d, e_ra, e_plx,/noerase,position=[xmargin,1-ymargin-3*ywidth,xmargin+1*xwidth,1-ymargin-2*ywidth], $
  charsize=charsize, ytitle=textoidl('\Delta\pi [mas]')

;;PMRA vs RA
bovy_plothist2d, e_ra,e_pmra,/noerase,position=[xmargin,1-ymargin-4*ywidth,xmargin+1*xwidth,1-ymargin-3*ywidth], $
  charsize=charsize, ytitle=textoidl('\Delta\mu_{\alpha} [mas yr^{-1}]')

;;PMDEC vs RA
bovy_plothist2d, e_ra,e_pmdec,/noerase,position=[xmargin,1-ymargin-5*ywidth,xmargin+1*xwidth,1-ymargin-4*ywidth], $
  charsize=charsize, ytitle=textoidl('\Delta\mu_{\delta} [mas yr^{-1}]'), xtitle=textoidl('\Delta\alpha [mas]'), xtickinterval=.5


;;PLX vs DEC
bovy_plothist2d, e_dec,e_plx,/noerase,position=[xmargin+xwidth,1-ymargin-3*ywidth,xmargin+2*xwidth,1-ymargin-2*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMRA vs DEC
bovy_plothist2d, e_dec,e_pmra,/noerase,position=[xmargin+xwidth,1-ymargin-4*ywidth,xmargin+2*xwidth,1-ymargin-3*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs DEC
bovy_plothist2d, e_dec,e_pmdec,/noerase,position=[xmargin+xwidth,1-ymargin-5*ywidth,xmargin+2*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\Delta\delta [mas]'), xtickinterval=.5


;;PMRA vs PLX
bovy_plothist2d, e_plx,e_pmra,/noerase,position=[xmargin+2*xwidth,1-ymargin-4*ywidth,xmargin+3*xwidth,1-ymargin-3*ywidth], $
  ytickname=REPLICATE(' ',30), xtickname=REPLICATE(' ',30)

;;PMDEC vs PLX
bovy_plothist2d, e_plx,e_pmdec,/noerase,position=[xmargin+2*xwidth,1-ymargin-5*ywidth,xmargin+3*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\Delta\pi [mas]'), xtickinterval=.5


;;PMDEC vs PMRA
bovy_plothist2d, e_pmra,e_pmdec,/noerase,position=[xmargin+3*xwidth,1-ymargin-5*ywidth,xmargin+4*xwidth,1-ymargin-4*ywidth], $
  ytickname=REPLICATE(' ',30), charsize=charsize, xtitle=textoidl('\Delta\mu_{\alpha} [mas yr^{-1}]'), xtickinterval=.5

ENDELSE


IF keyword_set(filename) THEN k_end_print

END
