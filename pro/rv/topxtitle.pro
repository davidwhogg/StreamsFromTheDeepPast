;+
;   NAME:
;      topxtitle
;   PURPOSE:
;      make the xtitle of a plot at the top of the plot
;   CALLING SEQUENCE:
;   INPUT:
;      yrange  -  y range
;      xrange  -  x range
;      xtitle  -  x title
;   KEYWORDS:
;   OUPUT:
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO TOPXTITLE, yrange=yrange, xtitle=xtitle, xrange=xrange, _EXTRA= _EXTRA

IF ~keyword_set(xrange) THEN xrange= !X.RANGE
IF ~keyword_set(yrange) THEN yrange= !Y.RANGE

djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, $
  xstyle=5, /noerase, _EXTRA=_EXTRA
axis, xaxis=0, xtickformat='(A1)'
axis, xaxis=1, xtitle=xtitle, _EXTRA=_EXTRA

END
