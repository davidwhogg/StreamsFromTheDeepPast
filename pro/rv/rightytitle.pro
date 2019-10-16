;+
;   NAME:
;      rightytitle
;   PURPOSE:
;      make the ytitle of a plot on the right of the plot
;   CALLING SEQUENCE:
;   INPUT:
;      yrange  -  y range
;      xrange  -  x range
;      ytitle  - y title
;   KEYWORDS:
;   OUPUT:
;   REVISION HISTORY:
;      2008-11-26 - Written Bovy (NYU)
;-
PRO RIGHTYTITLE, yrange=yrange, ytitle=ytitle, xrange=xrange, _EXTRA= _EXTRA

IF ~keyword_set(xrange) THEN xrange= !X.RANGE
IF ~keyword_set(yrange) THEN yrange= !Y.RANGE

djs_plot, [0.], [0.], xrange=xrange, yrange=yrange, $
  ystyle=5, /noerase, _EXTRA=_EXTRA
axis, yaxis=0, ytickformat='(A1)'
axis, yaxis=1, ytitle=ytitle, _EXTRA=_EXTRA

END
