;+
;   NAME:
;      check_quant
;   PURPOSE:
;      check the quantiles of the posterior distribution of the radial
;      velocities
;   CALLING SEQUENCE:
;      check_quant, ...
;   INPUT:
;      K     - Number of Gaussians
;      w     - regularization parameter
;      sample- sample to use (6=full)
;      plotfilename - filename to plot to
;      savefilename - filename to save the quantiles in
;      basedir - basedirectory for density estimate (ends in /)
;      yrange  - prob range
;      rrange  - vr range
;      xsize   - x size for k_print (can be left unset)
;      ysize   - ysize for k_print
;      offset  - tiny offset to add to vr such that it is
;                plotted when it is zero
;   KEYWORDS:
;      /colour - plot in color (in the density estimate there is a
;                 member color)
;      /QA     - Print QA info
;   OUTPUT:
;      plot of the quantiles
;   EXAMPLE:
;
;   REVISION HISTORY:
;      2009-01-16 - Written Bovy (NYU)
;-
PRO CHECK_QUANT, K=K, w=w, sample=sample, plotfilename=plotfilename, $
                 basedir=basedir, yrange=yrange, rrange=rrange, $
                 xsize=xsize, ysize=ysize, offset=offset, $
                 colour=colour, QA=QA, savefilename=savefilename

;;Default
if ~keyword_set(sample) THEN sample=6
IF ~keyword_set(basedir) THEN basedir='~/streams/pro/rv/ALMSv2/'
IF ~keyword_set(yrange) THEN yrange=[0,0.055]
IF ~keyword_set(rrange) THEN rrange=[-190,190]
IF ~keyword_set(offset) THEN offset= 0.00000001D;;Make's sure that vr's that are zero get plotted

IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    ;;Create the string to restore the density distribution
    densityfile=basedir+'fitv'
    densityfile+= strtrim(string(K),2)
    densityfile+='_bestV_V'
    densityfile+=strtrim(string(w,format='(f4.1)'),2)
    densityfile+= '_sample'
    densityfile+= strtrim(string(sample),2)
    densityfile+= '.sav'
    IF file_test(densityfile) THEN BEGIN
        IF keyword_set(QA) THEN splog, 'reading '+densityfile
        restore, densityfile
    ENDIF ELSE BEGIN
        splog, 'Density estimate not found in file '+densityfile
        splog, 'Returning...'
        RETURN
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

    ;;Calculate quantiles
    quantiles= dblarr(ngcs)
    FOR ii=0L, ngcs-1 DO BEGIN
        predictrv,amp,mean,covar,gcs[ii].sm,vr=gcs[ii].vr+offset,$
          vlvb=[[gcs[ii].vl],[gcs[ii].vb]],cvlvb=gcs[ii].vlvbc,$
          cvr=gcs[ii].vrc,loglike=loglike,altloglike=altloglike, $
          quantile=quantile, yrange=yrange, rrange=rrange
        quantiles[ii]= quantile
    ENDFOR
    ;;Save
    IF keyword_set(savefilename) THEN BEGIN
        splog, "Saving the result in "+savefilename
        save, filename=savefilename, quantiles, ngcs
    ENDIF
ENDELSE

;;Plot
weight = make_array(ngcs,/float,Value=1/float(ngcs))


IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize

;;Plotting stuff
!p.charsize=.6
!p.charthick=2
!x.ticks= 5
!x.minor= 4
xtickname=['0.0','0.2','0.4','0.6','0.8','1.0']
ytickname=REPLICATE('',30)
xtitle=textoidl('f(<v_r)')
!y.minor= 2

djs_plot, [0.], [0.], xrange=[0,1], $
  xtickname=xtickname, ytickname=ytickname, xtitle=xtitle, $
  ytitle=textoidl('P(f)'), yrange=[0,1.5]
hogg_plothist, quantiles, weight=weight, /overplot

IF keyword_set(plotfilename) THEN k_end_print

END
