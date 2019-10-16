;+
; NAME:
;   oplot_1d_gaussian
; PURPOSE:
;   overplot a mixture of 1-d gaussians on a 1-d plot
; INPUTS:
;   amp      - [N] amplitudes (amp=1 means the integral will be unity)
;   mean     - [N] means
;   var      - [N] variances or sigma^2
; OPTIONAL INPUTS:
;   numpts   - number of points in total plot
;   [etc]    - optional inputs to pass to oplot
; BUGS:
; REVISION HISTORY:
;   2003-02-21  written - Hogg
;-
pro oplot_1d_gaussian,amp,mean,var,numpts=numpts, $
                      _EXTRA=KeywordsForOplot
if NOT keyword_set(numpts) then begin
    TINY= 1D-10
    numpts= 1D1*abs(!X.CRANGE[1]-!X.CRANGE[0])/min(sqrt(var)+TINY)
    numpts= (round(numpts) > 100) < 100000
endif
ngauss= n_elements(amp)
dx= (!X.CRANGE[1]-!X.CRANGE[0])/double(numpts)
x= !X.CRANGE[0]+ dx*(dindgen(numpts)+0.5)
dx= abs(dx)
y= 0.0*x
for ii=0L,ngauss-1 do begin
    if sqrt(var[ii]) LT dx then begin
        oplot, [1,1]*mean[ii],!Y.CRANGE, $
          color=djs_icolor('red'),_EXTRA=KeywordsForOplot
    endif
    y= y+amp[ii]/sqrt(2.0*!DPI*var[ii])*exp(-0.5*(x-mean[ii])^2/var[ii])
endfor
oplot, x,y,_EXTRA=KeywordsForOplot
end
