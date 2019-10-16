;+
;   NAME:
;      plot_hyades_cmd
;   PURPOSE:
;      plot the cmd of the Hyades open cluster
;   CALLING SEQUENCE:
;   INPUT:
;   OUTPUT:
;   REVISION HISTORY:
;      2009-09-14 - Written - Bovy (NYU)
;-
PRO PLOT_HYADES_CMD, savefilename=savefilename,plotfilename=plotfilename

IF ~file_test(savefilename) THEN BEGIN
    hyades= read_hyades(savefilename=savefilename,/single,/tenpc)
ENDIF ELSE BEGIN
    restore, savefilename
ENDELSE

;;Plot the CMD
k_print, filename=plotfilename, xsize=3.375, ysize=4
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill
djs_plot, hyades.bvcolor, hyades.vmag-5.*alog10(100./hyades.plx), $
  xtitle='B-V [mag]', ytitle='M_{V} [mag]',xrange=[-0.3,1.6], $
  xcharsize=1,ycharsize=1, position=[0.1,0.1,0.9,0.9], yrange=[10,-6], $
  title='Hyades cluster',psym=8,symsize=0.5
;;overplot the isochrone
datafilename='hyades_isochrone2.dat'
iso= read_isochrone(datafilename)
djs_oplot, iso.h09-iso.h10, iso.h10,psym=-3


END
