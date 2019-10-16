;+
;   NAME:
;      plot_loo_kernel
;   PURPOSE:
;      plot the leave-one-out likelihood for different widths of the
;      kernel
;   INPUT:
;      wrange - width range
;      nws - number of ws in the range
;      plotfilename - filename for plot
;      xsize - xsize for k_print
;      ysize - ysize for k_print
;      savefilename - filename for saving
;   KEYWORDS:
;      rangelog - consider a logarithmic range (base 10)
;      log - plot the log of the loo likelihood
;      kernel - which kernel to use (tricube,epanechnikov, or gaussian)
;   OUTPUT:
;   REVISION HISTORY:
;      2009-11-05 - Written - Bovy (NYU)
;-
PRO PLOT_LOO_KERNEL, wrange=wrange, nws=nws,rangelog=rangelog, $
                     log=log, kernel=kernel, plotfilename=plotfilename, $
                     xsize=xsize, ysize=ysize, savefilename=savefilename
IF ~keyword_set(wrange) THEN wrange=[.001,.1]
IF ~keyword_set(nws) THEN nws= 3
IF ~keyword_set(kernel) THEN kernel='tricube'
IF ~keyword_set(xsize) THEN xsize= 3.375
IF ~keyword_set(ysize) THEN ysize= 5.
IF keyword_set(savefilename) AND file_test(savefilename) THEN BEGIN
    splog, "Restoring savefile "+savefilename
    restore, savefilename
ENDIF ELSE BEGIN
    wstep= (wrange[1]-wrange[0])/(nws-1.)
    ws= dindgen(nws)*wstep+wrange[0]
    IF keyword_set(rangelog) THEN ws= 10.^ws
    logls= dblarr(nws)
    ww= 0L
ENDELSE
WHILE ww LT nws DO BEGIN
    splog, "Working on iteration "+strtrim(string(ww),2)+"/"+strtrim(string(nws-1),2)
    logls[ww]= loo_kernel(ws[ww],kernel=kernel,/silent)
    ww+= 1L
    IF keyword_set(savefilename) THEN save, filename=savefilename, ws, logls, ww
ENDWHILE
IF kernel EQ 'tricube' THEN BEGIN
    title= 'Tricube'
ENDIF ELSE IF kernel EQ 'epanechnikov' THEN BEGIN
    title= 'Epanechnikov'
ENDIF ELSE BEGIN
    title= 'Gaussian'
ENDELSE
title+= ' kernel: '
IF keyword_set(log) THEN title+= 'max(log L) = '+strtrim(string(round(max(logls)),format='(I)'),2) ELSE title+= 'max(log L) = '+strtrim(string(round(max(logls)),format='(I)'),2) 
IF keyword_set(rangelog) THEN ws= alog10(ws)
IF keyword_set(rangelog) THEN xtitle='log_{10} \lambda [mag]' ELSE xtitle='\lambda [mag]'
IF keyword_set(log) THEN ytitle='\Delta log L(\lambda)' ELSE ytitle='L(\lambda)/L_{max}'
k_print, filename=plotfilename, xsize=xsize, ysize=ysize
charsize=.6
;;Play nice
isinfinite= where(FINITE(logls,/INFINITY) OR FINITE(logls,/NAN))
isfinite= where(FINITE(logls))
IF isinfinite[0] NE -1 THEN logls[isinfinite]= min(logls[isfinite])
if keyword_set(log) THEN BEGIN
    djs_plot, ws, logls-max(logls), xtitle=xtitle, ytitle=ytitle, charsize=charsize, title=title
ENDIF ELSE BEGIN
    djs_plot, ws, exp(logls-max(logls)), xtitle=xtitle, ytitle=ytitle, charsize=charsize, title=title
ENDELSE    
k_end_print
END
