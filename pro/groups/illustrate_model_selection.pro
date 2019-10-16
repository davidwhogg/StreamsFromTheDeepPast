;+
;   NAME:
;      illustrate_model_selection
;   PURPOSE:
;      to illustrate model selection
;   CALLING SEQUENCE:
;   INPUT:
;      plotfilename - filename for plot
;   OUTPUT:
;   REVISION HISTORY:
;      2009-09-23 - Written - Bovy (NYU)
;-
PRO ILLUSTRATE_MODEL_SELECTION, plotfilename=plotfilename
;;Two Gaussians, one narrow, the other one broad
sigma_one= 1.
sigma_two= 5.
nsamples= 1001
xs= dindgen(nsamples)/(nsamples-1)*25.-12.5
gauss_one= dblarr(nsamples)
gauss_two= dblarr(nsamples)
true_one= 1.
true_two= 5.
FOR ii=0L, nsamples-1 DO BEGIN
    gauss_one[ii]= 1./sqrt(2.*!DPI)/sigma_one*exp(-0.5*xs[ii]^2/sigma_one^2)
    gauss_two[ii]= 1./sqrt(2.*!DPI)/sigma_two*exp(-0.5*xs[ii]^2/sigma_two^2)
ENDFOR
;;Now plot
IF keyword_set(plotfilename) THEN k_print, filename=plotfilename
djs_plot, xs, gauss_one, xtickname=strarr(30)+' ', ytickname=strarr(30)+' ',$
  position=[0.1,0.5,0.5,0.9], xticks=1, yticks=1,thick=2, $
  ytitle='p(x!10#!Xmodel)'
xyouts, 0, -.04,'x', charsize=2
djs_oplot, xs, gauss_two, color='gray', thick=6
oplotbarx, true_one, thick=4.
djs_plot, xs, gauss_one, xtickname=strarr(30)+' ', ytickname=strarr(30)+' ',$
  position=[0.5,0.5,0.9,0.9], xticks=1, yticks=1, /NOERASE,thick=2
xyouts, 0, -.04,'x', charsize=2
djs_oplot, xs, gauss_two, color='gray',thick=6
oplotbarx, true_two, thick=4.
IF keyword_set(plotfilename) THEN k_end_print

END
