;+
;
;-
FUNCTION tricube_kernel, t
;;No checking whether this is in the range t<0
out= dblarr(n_elements(t))
FOR ii=0L, n_elements(t)-1 DO out[ii]= (1.-t[ii]^3.)^3.
RETURN, out
END
PRO ILLUSTRATE_NONPARAM, plotfilename=plotfilename

;;Restore the data
restore, filename='../../lsr2/hip2-aumer.sav'
nhip= n_elements(hip.hip)
;;Select random sample of 100 stars
seed=-1
indx= randomu(seed,100)*nhip

k_print, filename=plotfilename+"_1.ps"

phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill

djs_plot, hip[indx].bvcolor, hip[indx].vmag-$
  5.*alog10(100./hip[indx].plx), $
  psym=8,symsize=0.5, $
  xrange=[-0.3,1.6], yrange=[10,-6]

XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Color [mag]',charthick=10, color=djs_icolor('black'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'Absolute Magnitude [mag]', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('black'),charsize=2

;;Overplot Tricube
width=0.1
xs= dindgen(100)*2*width/100.-width
ys= -tricube_kernel(abs(xs/width))*4.+10
xs+= 0.5
djs_oplot, xs, ys,thick=10.

k_end_print

k_print, filename=plotfilename+"_2.ps"
djs_plot, hip[indx].bvcolor, hip[indx].vmag-$
  5.*alog10(100./hip[indx].plx), $
  psym=8,symsize=0.5, $
  xrange=[-0.3,1.6], yrange=[10,-6]

XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Color [mag]',charthick=10, color=djs_icolor('black'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'Absolute Magnitude [mag]', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('black'),charsize=2

;;Overplot Tricube
width=0.1
xs= dindgen(100)*2*width/100.-width
ys= -tricube_kernel(abs(xs/width))*4.+10
xs+= 0.5
djs_oplot, xs, ys,thick=10.

thiship= hip[indx]
thesestars= thiship[where(thiship.bvcolor GT 0.5-width AND thiship.bvcolor LT 0.5+width)]
nstars= n_elements(thesestars.hip)
FOR ii=0L, nstars-1 DO BEGIN
    djs_oplot, [thesestars[ii].bvcolor,thesestars[ii].bvcolor],$
      [-tricube_kernel(abs((thesestars[ii].bvcolor-0.5)/width))*4.+10.,thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx)],color='red',thick=tricube_kernel(abs((thesestars[ii].bvcolor-0.5)/width))*4.
ENDFOR

k_end_print


k_print, filename=plotfilename+"_3.ps"

ii=0
djs_plot, [thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx),thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx)],$
  [0.,tricube_kernel(abs((thesestars[ii].bvcolor-0.5)/width))*4.],xrange=[2.,6.],yrange=[0.,4.]
FOR ii=1L, nstars-1 DO BEGIN
    djs_oplot, [thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx),thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx)],$
      [0.,tricube_kernel(abs((thesestars[ii].bvcolor-0.5)/width))*4.]
ENDFOR
XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Absolute Magnitude M_V [mag]',charthick=10, color=djs_icolor('black'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'p(M_V)', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('black'),charsize=2
k_end_print



k_print, filename=plotfilename+"_4.ps"

width=0.1
thesestars= hip[where(hip.bvcolor GT 1.-width AND hip.bvcolor LT 1.+width)]
nstars= n_elements(thesestars.hip)
ii=0
djs_plot, [thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx),thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx)],$
  [0.,tricube_kernel(abs((thesestars[ii].bvcolor-1.)/width))*4.], yrange=[0.,4.],xrange=[4.,9.]
FOR ii=1L, nstars-1 DO BEGIN
    djs_oplot, [thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx),thesestars[ii].vmag- 5.*alog10(100./thesestars[ii].plx)],$
      [0.,tricube_kernel(abs((thesestars[ii].bvcolor-1.)/width))*4.]
ENDFOR
XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Absolute Magnitude M_V [mag]',charthick=10, color=djs_icolor('black'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'p(M_V)', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('black'),charsize=2
k_end_print


k_print, filename=plotfilename+"_5.ps"

width=0.1
thesestars= hip[where(hip.bvcolor GT 1.-width AND hip.bvcolor LT 1.+width)]
nstars= n_elements(thesestars.hip)
xrange=[4.,9.]
smoothed_xs= dindgen(1001)/1000.*5+4.
amp= dblarr(nstars)
for ii=0L, nstars-1 do begin
    amp[ii]= tricube_kernel(abs((thesestars[ii].bvcolor-1.)/width))*4.
endfor
amp= amp/total(amp)
mean= thesestars.vmag- 5.*alog10(100./thesestars.plx)
covar= dblarr(nstars)+0.1
smoothed= oned_sum_gaussians(smoothed_xs,mean,covar,amp)
djs_plot, smoothed_xs, smoothed;, yrange=[0.,4.],xrange=[4.,9.]
XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Absolute Magnitude M_V [mag]',charthick=10, color=djs_icolor('black'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'p(M_V)', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('black'),charsize=2
k_end_print



END
