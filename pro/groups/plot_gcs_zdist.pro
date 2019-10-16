;+
;   NAME:
;      plot_gcs_zdist
;   PURPOSE:
;      plot the metallicity distribution of the photometric GCS sample
;   INPUT:
;      plotfilename - filename for plot
;      ngauss - number of gaussians to decompose the distribution as
;   OUTPUT:
;      a plot
;   HISTORY:
;      2009-11-20 - Written - Bovy (NYU)
;-
PRO PLOT_GCS_ZDIST, plotfilename=plotfilename, ngauss=ngauss

IF ~keyword_set(ngauss) THEN ngauss=2
restore, 'phot_gcs.sav'
ngcs= n_elements(gcs.hip)
err= 0.08;;Typical GCS error

;;Calculate best fit sum of 2 gaussians using extreme deconvolution
ydata= gcs.feh
ydata= reform(ydata,1,ngcs)
ycovar= dblarr(ngcs)+err^2.
ycovar= reform(ycovar,1,1,ngcs)
projection= dblarr(ngcs)+1.
projection= reform(projection,1,1,ngcs)
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


splog, xmean, ydata[0,0]
projected_gauss_mixtures_c, ngauss, ydata, ycovar, amp, xmean, $
  xcovar, fixamp=fixamp, avgloglikedata=loglike, /quiet, logfile='ed_test'
splog, xmean
splog, "Loglike was"
splog, loglike
xmean= reform(xmean,ngauss)
xcovar= reform(xcovar,ngauss)

;;Overplot best fit composition
xrange= minmax(gcs.feh)
nxs= 1001
xstep= (xrange[1]-xrange[0])/nxs
xs= dindgen(1001)*xstep+xrange[0]
ys= oned_sum_gaussians(xs,xmean,xcovar+err^2.,amp)
ys1= .716666*oned_sum_gaussians(xs,xmean[0],xcovar[0]+err^2.,amp[0])
ys2= oned_sum_gaussians(xs,xmean[1],xcovar[1]+err^2.,amp[1])
IF ngauss EQ 3 THEN ys3= oned_sum_gaussians(xs,xmean[2],xcovar[2]+err^2.,amp[2])


;;Add legend
thinlegend= textoidl('[Fe/H]_{thin} = '+strtrim(string(xmean[0],format='(F5.2)'),2)+$
                     ' !9+!X ' + strtrim(string(sqrt(xcovar[0]),format='(F5.2)'),2)+' dex')
thicklegend= textoidl('[Fe/H]_{thick} = '+strtrim(string(xmean[1],format='(F5.2)'),2)+$
                     ' !9+!X ' + strtrim(string(sqrt(xcovar[1]),format='(F5.2)'),2)+' dex')
fthinlegend= textoidl('f_{thin} = '+strtrim(string(amp[0],format='(F5.3)'),2))
IF ngauss EQ 3 THEN BEGIN
    halolegend= textoidl('[Fe/H]_{halo} = '+strtrim(string(xmean[2],format='(F5.2)'),2)+$
                     ' !9+!X ' + strtrim(string(sqrt(xcovar[2]),format='(F5.2)'),2)+' dex')
    fhalolegend= textoidl('f_{halo} = '+strtrim(string(amp[2],format='(F5.3)'),2))
ENDIF

weights= dblarr(n_elements(gcs.hip))+1./n_elements(gcs.hip)
k_print, filename=plotfilename, xsize=3.375, ysize=4 ; actual width of ApJ column!
hogg_plothist, gcs.feh, npix=100,$
  xtitle='[Fe/H]', position=[0.1,0.1,0.9,0.8],$
  xcharsize=1,ycharsize=1, weight=weights
djs_oplot, xs, ys
djs_oplot, xs, ys1, linestyle=2
djs_oplot, xs, ys2, linestyle=2
IF ngauss EQ 3 THEN djs_oplot, xs, ys3, linestyle=2
legendcharsize= .75
xposlegend= -2.35
yposlegend= 2.7
yposincr= -.2
legend, [thinlegend], box=0., position=[xposlegend,yposlegend], charsize= legendcharsize
legend, [thicklegend], box=0., position=[xposlegend,yposlegend+yposincr], charsize= legendcharsize
IF ngauss EQ 3 THEN BEGIN
    legend, [halolegend], box=0., position=[xposlegend,yposlegend+2*yposincr], charsize= legendcharsize
    legend, [fthinlegend], box=0., position=[xposlegend,yposlegend+3*yposincr], charsize= legendcharsize
    legend, [fhalolegend], box=0., position=[xposlegend,yposlegend+4*yposincr], charsize= legendcharsize
ENDIF ELSE BEGIN    
    legend, [fthinlegend], box=0., position=[xposlegend,yposlegend+2*yposincr], charsize= legendcharsize
ENDELSE
;;And some arrows
plots, [-1.5,-.8], [1.,.6],/data
legend, ['thick'], pos=[-2.,1.2],box=0, charsize=legendcharsize
plots, [-0.01,.22], [1.35,2.],/data
legend, ['thin'], pos=[0.,2.2],box=0, charsize=legendcharsize
k_end_print
END
