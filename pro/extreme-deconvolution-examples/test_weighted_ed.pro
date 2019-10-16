;+
;   NAME:
;      test_weighted_ed
;   PURPOSE:
;      simple testing code for the weighted extreme-deconvolution
;      code:
;      draw samples from a mixture of two Gaussians, calculate weights,
;      fit this with and without weights
;   INPUT: (optional)
;      seed - seed for the random number generator
;      n1 - number of samples to draw for the first data set (one
;           Gaussian)
;      n2 - number of samples to draw for the second data set (mixture
;           of two Gaussians, first one equal to that of data1)
;      plot - plot histograms to output
;   OUTPUT:
;   HISTORY:
;      2010-04-01 - Written - Bovy (NYU)
;-
PRO TEST_WEIGHTED_ED, seed=seed, n1=n1, n2=n2, plot=plot
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(n1) THEN n1= 100
IF ~keyword_set(n2) THEN n2= 100
data1= randomn(seed,n1)
mean2= 5.
sig2= 2.
data2= dblarr(n2)
weight2= dblarr(n2)
FOR ii=0L, n2-1 DO BEGIN
    onequestionmark= randomu(seed)
    if onequestionmark LT 0.5 then data2[ii]= randomn(seed) else data2[ii]= randomn(seed)*sig2+mean2
    weight2[ii]= exp(-0.5*data2[ii]^2.)/(exp(-0.5*data2[ii]^2.)+1./sig2*exp(-0.5*(data2[ii]-mean2)^2./sig2^2.))
ENDFOR

print, "Input mean and variance for the second Gaussian", mean2, sig2^2.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Prepare for the extreme-deconvolution fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ydata= data2
ydata= reform(ydata,1,n2)
ycovar= dblarr(n2)
ycovar= reform(ycovar,1,n2)
weight= alog(weight2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Initial conditions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ngauss= 1
amp= [1.]
xmean= [0.]
xmean= reform(xmean,1,ngauss)
xcovar= [1.]
xcovar= reform(xcovar,1,1,ngauss)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;RUN extreme-deconvolution
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
projected_gauss_mixtures_c, ngauss, ydata, ycovar, amp, xmean, $
  xcovar, avgloglikedata=loglike, /quiet, weight=weight,/logweight
print, "Loglike was"
print, loglike
xmean= reform(xmean,ngauss)
xcovar= reform(xcovar,ngauss)
print, "Output mean and variance with weights", xmean[0], xcovar[0]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Initial conditions, again
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ngauss= 1
amp= [1.]
xmean= [0.]
xmean= reform(xmean,1,ngauss)
xcovar= [1.]
xcovar= reform(xcovar,1,1,ngauss)
projected_gauss_mixtures_c, ngauss, ydata, ycovar, amp, xmean, $
  xcovar, avgloglikedata=loglike, /quiet
print, "Output mean and variance without weights", xmean[0], xcovar[0]
END
