;+
;   NAME:
;      plot_mscuts
;   PURPOSE:
;      plot the main-sequence cuts that we applied to define the MS
;      sample
;      overplots these on the current figure
;   INPUT:
;      color - color of the lines
;   OUTPUT:
;      plots on the current figure
;   HISTORY:
;      2009-12-14 - Written - Bovy
;-
PRO PLOT_MSCUTS, color=color
if ~keyword_set(color) THEN color=djs_icolor('gray')
;;upper limits
djs_oplot, [!X.CRANGE[0],0.5],[7.5*!X.CRANGE[0]-3.75,0],color=color
djs_oplot, [0.5,.8],[0.,4.6], color=color
djs_oplot, [0.8,!X.CRANGE[1]],[4.6,4.43*!X.CRANGE[1]+1.055],color=color
;;lower limits
djs_oplot, [!X.CRANGE[0],0.35],[4.62*!X.CRANGE[0]+2.383,4], color=color
djs_oplot, [0.35,0.65],[4,6.5],color=color
djs_oplot, [.65,1.25],[6.5,8.5],color=color
djs_oplot, [1.25,!X.CRANGE[1]],[8.5,6.5*!X.CRANGE[1]+0.375],color=color
END
