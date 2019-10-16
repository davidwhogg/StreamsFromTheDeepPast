;+
;   NAME:
;      plot_gcs_cmd
;   PURPOSE:
;      plot the CMD of the photometric GCS sample
;   INPUT:
;      plotfilename - filename for plot
;   OUTPUT:
;      a plot
;   HISTORY:
;      2009-11-20 - Written - Bovy (NYU)
;-
PRO PLOT_GCS_CMD, plotfilename=plotfilename

restore, 'phot_gcs.sav'

k_print, filename=plotfilename, xsize=3.375, ysize=4 ; actual width of ApJ column!
djs_plot, gcs.bvcolor, gcs.vmag-5.*alog10(100./gcs.plx),$
  yrange=[8,-1], color=0,$
  xtitle='B-V [mag]', ytitle='M_{V} [mag]',xrange=[0.1,1.2], $
  xcharsize=1,ycharsize=1, position=[0.1,0.1,0.9,0.9],$
  title='Geneva-Copenhagen sample', psym=3
legend, [textoidl('\pi/\sigma_\pi \geq 10')], box=0.,/right
k_end_print
END
