;+
;   NAME:
;      hist_gcslike
;   PURPOSE:
;      plot histograms of the gcs likelihood for the fiducial model
;   CALLING SEQUENCE:
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      filename - filename for plot
;      basedir - basedirectory for density estimate (ends in /)
;      savefilename - save (or get) the likelihoods here
;      xsize    - xsize for k_print
;      ysize    - ysize for k_print
;   KEYWORDS:
;      /QA     - Print QA info
;   OUTPUT:
;      plot
;   REVISION HISTORY:
;      2009-01-19 - Written Bovy (NYU)
;-
PRO HIST_GCSLIKE, K=K, w=w, sample=sample, filename=filename, $
                  basedir=basedir, savefilename=savefilename, QA=QA, $
                  xsize=xsize, ysize=ysize

;;Defaults
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'


IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    ;;run tile_goodbadpredict to construct the likelihoods
    tile_goodbadpredict, K=K, w=w, sample=sample, basedir=basedir, $
      QA=QA, savefilename=savefilename
    restore, savefilename;;likelihoods are now in loglikes
ENDELSE


;;Restore GCS
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


;;Now plot the histograms
IF keyword_set(filename) THEN k_print, filename=filename, $
  xsize=xsize, ysize=ysize

!p.multi=[0,1,2]
charsize=.6
charthick=2
!p.charsize=charsize
!p.charthick=charthick

lrange=[-10.,max(loglikes)]

xmargin=.1
ymargin=.1
spacing=.025


;;Calculate histogram
weight= dblarr(ngcs)+1D/ngcs
hogg_plothist, loglikes, weight=weight, /dontplot, hist=hist, $
  /totalweight, xrange=lrange
yrange=[0.,(1+.1)*max(hist)]
djs_plot, [0.],[0.],xrange=lrange, yrange=yrange, $
  ytitle='Fraction of stars in GCS', $
  position=[xmargin,$
            ymargin+(1D0-2.*ymargin-spacing)/2.+spacing,$
            xmargin+(1.-2.*xmargin),$
            1.-ymargin], $
  xtickname=replicate(' ',30)
hogg_plothist, loglikes, weight=weight, /overplot, /totalweight

;;Now do the 2d histogram
bovy_plothist2d, loglikes, gcs.vr, xtitle='ln(likelihood)', $
  ytitle=textoidl('v_r [km s^{-1}]'), $
  position=[xmargin,$
            ymargin, $
            xmargin+(1.-2.*xmargin),$
            ymargin+(1D0-2.*ymargin-spacing)/2.], $
  charsize=charsize, charthick=charthick, /log,$
  xrange=lrange, nbins=[16,8], yrange=[-120,120]

IF keyword_set(filename) THEN k_end_print

END
