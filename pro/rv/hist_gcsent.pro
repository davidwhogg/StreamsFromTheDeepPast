;+
;   NAME:
;      hist_gcsent
;   PURPOSE:
;      plot histograms of the gcs entropy for the fiducial model
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      likesavefilename - save (or get) the likelihoods here
;      entsavefilename - save (or get) the entropies here
;      xsize    - xsize for k_print
;      ysize    - ysize for k_print
;   KEYWORDS:
;      /QA     - Print QA info
;   OUTPUT:
;      plot
;   REVISION HISTORY:
;      2009-01-19 - Written Bovy (NYU)
;-
PRO HIST_GCSENT, K=K, w=w, sample=sample, filename=filename, $
                  basedir=basedir, $
                 likesavefilename=likesavefilename, QA=QA, $
                 entsavefilename=entsavefilename, $
                  xsize=xsize, ysize=ysize

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'


;;Get entropies
IF file_test(entsavefilename) THEN BEGIN
    splog, "Restoring savefile "+entsavefilename
    restore, entsavefilename
ENDIF ELSE BEGIN
    ;;run tile_highlowpredict to construct the likelihoods
    tile_highlowpredict, K=K, w=w, sample=sample, basedir=basedir, $
      QA=QA, savefilename=entsavefilename
    restore, entsavefilename;;entropies are now in entropies
ENDELSE

;;Same for likelihoods
IF file_test(likesavefilename) THEN BEGIN
    splog, "Restoring savefile "+likesavefilename
    restore, likesavefilename
ENDIF ELSE BEGIN
    ;;run tile_goodbadpredict to construct the likelihoods
    tile_goodbadpredict, K=K, w=w, sample=sample, basedir=basedir, $
      QA=QA, savefilename=likesavefilename
    restore, likesavefilename;;likelihoods are now in loglikes
ENDELSE

;;Restore GCS
;genevafilename= '/global/data/rv/gcs/gcs.sav'
;IF file_test(genevafilename) THEN BEGIN
;    restore, genevafilename
;ENDIF ELSE BEGIN
;    splog, 'Error: GCS-savefile does not exist in ', $
;      genevafilename
;    splog, 'GCS-savefile must exist, returning...'
;    RETURN
;ENDELSE
;ngcs= n_elements(gcs.hip)


;;Now plot the histograms
IF keyword_set(filename) THEN k_print, filename=filename, $
  xsize=xsize, ysize=ysize

charsize=.6
charthick=2
!p.charsize=charsize
!p.charthick=charthick

erange=minmax(entropies)

xmargin=.1
ymargin=.1
spacing=.025


;;Calculate histogram
weight= dblarr(ngcs)+1D/ngcs
hogg_plothist, entropies, weight=weight, /dontplot, hist=hist, $
  /totalweight, xrange=erange
yrange=[0.,(1+.1)*max(hist)]
djs_plot, [0.],[0.],xrange=erange, yrange=yrange, $
  ytitle='Fraction of stars in GCS', $
  position=[xmargin,$
            ymargin+(1D0-2.*ymargin-spacing)/2.+spacing,$
            xmargin+(1.-2.*xmargin),$
            1.-ymargin], $
  xtickname=replicate(' ',30)
hogg_plothist, entropies, weight=weight, /overplot, /totalweight

;;Now do the 2d histogram
bovy_plothist2d, entropies,loglikes, xtitle='H', $
  ytitle='ln(likelihood)', $
  position=[xmargin,$
            ymargin, $
            xmargin+(1.-2.*xmargin),$
            ymargin+(1D0-2.*ymargin-spacing)/2.], $
  charsize=charsize, charthick=charthick,/log,$
  xrange=erange, nbins=[16,16], yrange=[-10.,max(loglikes)], /NOERASE

IF keyword_set(filename) THEN k_end_print

END
