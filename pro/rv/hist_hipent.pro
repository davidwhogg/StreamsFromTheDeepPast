;+
;   NAME:
;      hist_hipent
;   PURPOSE:
;      plot the histogram of the hip\gcs entropy for the fiducial model
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      savefilename - save (or get) the entropies here
;      xsize    - xsize for k_print
;      ysize    - ysize for k_print
;   KEYWORDS:
;      /QA     - Print QA info
;   OUTPUT:
;      plot
;   REVISION HISTORY:
;      2009-01-19 - Written Bovy (NYU)
;-
PRO HIST_HIPENT, K=K, w=w, sample=sample, filename=filename, $
                  basedir=basedir, $
                 savefilename=savefilename, QA=QA, $
                 xsize=xsize, ysize=ysize

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'


;;Get entropies
IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    ;;run tile_highlowpredict to construct the likelihoods
    tile_infohip, K=K, w=w, sample=sample, basedir=basedir, $
      QA=QA, savefilename=savefilename
    restore, savefilename;;entropies are now in entropies
ENDELSE


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
weight= dblarr(n_notingcs)+1D/n_notingcs
hogg_plothist, entropies, weight=weight, /dontplot, hist=hist, $
  /totalweight, xrange=erange
yrange=[0.,(1+.1)*max(hist)]
djs_plot, [0.],[0.],xrange=erange, yrange=yrange, $
  ytitle='Fraction of stars', $
  xtitle='H'
hogg_plothist, entropies, weight=weight, /overplot, /totalweight

IF keyword_set(filename) THEN k_end_print

END
