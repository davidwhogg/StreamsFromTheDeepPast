;+
;   NAME:
;      check_quant_back_plx
;   PURPOSE:
;      check the quantiles of the distribution of the parallax
;   INPUT:
;      width - width parameter of the kernel
;      plotfilename - filename to plot to
;      savefilename - filename to save the quantiles in
;      xsize   - x size for k_print (can be left unset)
;      ysize   - ysize for k_print
;      datafilename - filename with the data
;   KEYWORDS:
;      kernel - which kernel to use (tricube,epanechnikov, or
;               gaussian)
;   OUTPUT:
;      plot of the quantiles
;   REVISION HISTORY:
;      2009-11-05 - Written Bovy (NYU)
;-
PRO CHECK_QUANT_BACK_PLX, width=width, kernel=kernel, plotfilename=plotfilename, $
                 xsize=xsize, ysize=ysize, savefilename=savefilename, $
                 datafilename=datafilename
;;Restore the Hipparcos data
IF ~keyword_set(datafilename) THEN datafilename='../lsr2/hip2-aumer.sav'
restore, filename=datafilename
nhip= n_elements(hip.hip)
IF file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    ;;Calculate quantiles
    quantiles= dblarr(nhip)
    ii=0L
ENDELSE
WHILE ii LT nhip DO BEGIN
    print, ii, nhip-1
    predictplx, ii, width=width, kernel=kernel, datafilename=datafilename, $
      /byindx,quantile=quantile,/dontplot
    quantiles[ii]= quantile
    ii+= 1
    save, filename=savefilename, quantiles, nhip, ii
ENDWHILE
;;Plot
weight = make_array(nhip,/float,Value=1/float(nhip))
IF keyword_set(plotfilename) THEN k_print, filename=plotfilename, $
  xsize=xsize, ysize=ysize
;;Plotting stuff
!p.charsize=.7
!p.charthick=2
!x.ticks= 5
!x.minor= 4
xtickname=['0.0','0.2','0.4','0.6','0.8','1.0']
ytickname=REPLICATE('',30)
xtitle=textoidl('f(<\pi)')
!y.minor= 2

djs_plot, [0.], [0.], xrange=[0,1], $
  xtickname=xtickname, ytickname=ytickname, xtitle=xtitle, $
  ytitle=textoidl('P(f)'), yrange=[0,1.5]
hogg_plothist, quantiles, weight=weight, /overplot
IF keyword_set(plotfilename) THEN k_end_print
END
