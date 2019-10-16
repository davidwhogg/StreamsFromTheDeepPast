;+
;   NAME:
;      plot_gcs_zdist
;   PURPOSE:
;      plot the metallicity distribution of the photometric GCS sample
;   INPUT:
;      plotfilename - filename for plot if ps is set
;      ngauss - number of gaussians to decompose the distribution as
;   KEYWORDS;
;      ps  - save as ps-file
;   OUTPUT:
;      a plot
;   HISTORY:
;      2009-02-02 - Written - Bovy (NYU)
;-
PRO PLOT_GCS_ZDIST, plotfilename=plotfilename, ngauss=ngauss, ps=ps

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Load the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF ~keyword_set(plotfilename) THEN plotfilename='gcs_zdist_test.ps'
IF ~keyword_set(ngauss) THEN ngauss=2
restore, 'test.sav'
ngcs= n_elements(test.feh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Prepare for the extreme-deconvolution fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ydata= test.feh
ydata= reform(ydata,1,ngcs)
ycovar= test.e_feh
ycovar= reform(ycovar,1,ngcs)
;projection= dblarr(ngcs)+1.
;projection= reform(projection,1,1,ngcs)
weight= dblarr(ngcs)+1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Initial conditions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF ngauss EQ 2 THEN BEGIN
    amp= [.5,.5]
    xmean= [0.,0.]
    xmean= reform(xmean,1,ngauss)
    xcovar= [1.,2.]
    xcovar= reform(xcovar,1,1,ngauss)
ENDIF ELSE IF ngauss EQ 3 THEN BEGIN
    amp= [.33,.33,.34]
    xmean= [0.,0.,0.]
    xmean= reform(xmean,1,ngauss)
    xcovar= [1.,2.,4.]
    xcovar= reform(xcovar,1,1,ngauss)
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;RUN extreme-deconvolution
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
projected_gauss_mixtures_c, ngauss, ydata, ycovar, amp, xmean, $
  xcovar, avgloglikedata=loglike, /quiet, weight=weight
print, "Loglike was"
print, loglike
xmean= reform(xmean,ngauss)
xcovar= reform(xcovar,ngauss)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Overplot best fit composition: Calculate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xrange=[-2.5,0.5]
nxs= 1001
xstep= (xrange[1]-xrange[0])/nxs
xs= dindgen(1001)*xstep+xrange[0]
ys= oned_sum_gaussians(xs,xmean,xcovar,amp)
ys1= .666666*oned_sum_gaussians(xs,xmean[0],xcovar[0],amp[0]);;Scale the thin disk down somewhat
ys2= oned_sum_gaussians(xs,xmean[1],xcovar[1],amp[1])
IF ngauss EQ 3 THEN ys3= oned_sum_gaussians(xs,xmean[2],xcovar[2],amp[2])

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;PLOTTING
;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF keyword_set(ps) THEN BEGIN
    axis_char_scale= 1.75
    tiny=1.d-4
    pold=!P
    xold=!X
    yold=!Y
    bangp=!P
    bangx=!X
    bangy=!Y
    !P.FONT= -1
    set_plot, "PS"
    !P.BACKGROUND= 16777215
    !P.COLOR= 0
    xsize= 3.375
    ysize= 4
    device, file=filename,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, $
      bits_per_pixel=64
    !P.THICK= 2.0
    !P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
    !P.CHARSIZE= 1.0
    !P.PSYM= 0
    !P.LINESTYLE= 0
ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;Prepare for histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;
binsize=0.03
nbins= (xrange[1]-xrange[0])/binsize
fehs= dindgen(nbins)*(xrange[1]-xrange[0])/double(nbins)+xrange[0];;X-axis
plot, fehs, Histogram(test.feh, BINSIZE=.03,min=xrange[0],max=xrange[1])/(double(ngcs)*binsize), PSYM=10, xtitle='[Fe/H]', position=[0.1,0.1,0.9,0.8],$
  xcharsize=1,ycharsize=1, xrange=xrange
oplot, xs, ys
oplot, xs, ys1, linestyle=2
oplot, xs, ys2, linestyle=2
IF ngauss EQ 3 THEN oplot, xs, ys3, linestyle=2
END
;;Assign members
xmean= reform(xmean,1,ngauss)
xcovar= reform(xcovar,1,1,ngauss)
assign_members, logpost, ydata,ycovar, projection, xmean,xcovar, amp
sortindx= sort(ydata)
plot, ydata[sortindx],exp(logpost[0,sortindx]),psym=-3, yrange=[0,1],xrange=xrange,$
  xtitle='[Fe/H]'
oplot, ydata[sortindx],exp(logpost[1,sortindx]),psym=-3, color=255
END
