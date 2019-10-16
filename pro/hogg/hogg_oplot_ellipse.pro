;+
; NAME:
;   hogg_oplot_ellipse
; PURPOSE:
;   overplot an ellipse given a 2-d mean and a 2-d covariance tensor
; CALLING SEQUENCE:
;   hogg_oplot_ellipse, mean,covar,[extras]
; INPUTS:
;   x,y      - [N] arrays of x,y centers of ellipses
;   covar    - [2,2,N] array covariance tensors
; OPTIONAL INPUTS:
;   color    - [N] of colors or a single color; default to 'default'
;   npts     - number of points to use per ellipse; default 24
;   [extras] - options to pass to oplot (NOT djs_oplot)
; COMMENTS:
;   Symmetrizes covar before plotting -- beware!
; BUGS:
;   Needs extensive testing.
; REVISION HISTORY:
;   2003-03-26  written - Hogg (NYU)
;-
pro hogg_oplot_ellipse, x,y,covar,color=color,npts=npts, $
                        _EXTRA=KeywordsForOplot
if n_elements(x) GT 1 then begin
    nellipse= n_elements(x)
    for ii=0L,nellipse-1 do begin
        if keyword_set(color) then if n_elements(color) GT 1 then $
          color1= color[ii] else color1=color
        hogg_oplot_ellipse, x[ii],y[ii],covar[*,*,ii], $
          color=color1,npts=npts,_EXTRA=KeywordsForOplot
    endfor
    return
endif
if NOT keyword_set(npts) then npts= 24
covar= 5d-1*(covar+transpose(covar))
eval= eigenql(covar,eigenvectors=evec,/double)
Rmatrix= (sqrt(eval)##[1d0,1d0])*evec
theta= 2d0*!DPI*dindgen(npts+1)/double(npts)
points= transpose([[cos(theta)],[sin(theta)]])
points= Rmatrix#points
if NOT keyword_set(color) then color='default'
oplot, x+points[0,*],y+points[1,*],color=djs_icolor(color), $
  _EXTRA=KeywordsForOplot,psym=0
return
end
