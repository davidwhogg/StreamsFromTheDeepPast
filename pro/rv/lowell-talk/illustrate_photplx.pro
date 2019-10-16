;+
;   NAME:
;      illustrate_photplx
;   PURPOSE:
;     illustrate a photometric parallax
;   INPUT:
;     plotfilename - filename for plot
;   OUTPUT:
;   HISTORY:
;      2010-01-11 - Written - Bovy
;-
PRO ILLUSTRATE_PHOTPLX, plotfilename=plotfilename

datafilename='../../groups/pleiades_isochrone.dat'

k_print, filename=plotfilename

iso= read_isochrone(datafilename)
djs_plot, iso.h09-iso.h10, iso.h10,psym=-3, thick=15,$
  xrange=[-0.3,1.6], yrange=[10,-6]
XYOutS, (!X.Window[1] - !X.Window[0]) / 2 + !X.Window[0], 0.05, /Normal, $
  Alignment=0.5, 'Color [mag]',charthick=10, color=djs_icolor('red'),charsize=2
XYOutS, 0.08,(!Y.Window[1] - !Y.Window[0]) / 2 + !Y.Window[0], 'Absolute Magnitude [mag]', /Normal, $
  Alignment=0.5, orientation=90.,charthick=10, color=djs_icolor('blue'),charsize=2
djs_oplot, [0.449,0.449],[10.,3.7210000],psym=-3, color='red',thick=20.
djs_oplot, [-0.3,0.449],[3.721,3.7210000],psym=-3, color='blue',thick=20.


k_end_print
END
