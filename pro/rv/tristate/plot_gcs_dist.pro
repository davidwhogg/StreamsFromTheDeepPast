;+
;   NAME:
;      plot_gcs_dist
;   PURPOSE:
;      plot the GCS UV distribution
;   INPUT:
;   OUTPUT:
;   REVISION HISTORY:
;      2009-10-22 - Written - Bovy (NYU)
;-
PRO plot_gcs_dist
;;Read the GCS
restore, '/global/data/rv/gcs/gcs.sav'
ngcs=n_elements(gcs.hip)
;;Rotate all velocities to UVW
UVW=dblarr(ngcs,3)
FOR ii=0L, ngcs-1 DO BEGIN
    rot1=dblarr(3,3)
    rot2=dblarr(3,3)
    l= gcs[ii].l/180.*!DPI
    b= gcs[ii].b/180.*!DPI
    rot1[0,0]= cos(l)
    rot1[1,0]= sin(l)
    rot1[0,1]= -sin(l)
    rot1[1,1]= cos(l)
    rot1[2,2]= 1
    rot2[0,0]= cos(b)
    rot2[2,0]= sin(b)
    rot2[0,2]= -sin(b)
    rot2[2,2]= cos(b)
    rot2[1,1]= 1
    thisUVW= transpose(rot1)##(transpose(rot2)##transpose([gcs[ii].vr,gcs[ii].vl,gcs[ii].vb]))
    UVW[ii,*]= thisUVW
ENDFOR
k_print, filename='gcs_dist.ps'
djs_plot, UVW[*,0],UVW[*,1],psym=3,xrange=[-130,120],yrange=[-120,60],xtitle='v_x [km s^{-1}]', ytitle='v_y [km s^{-1}]',/isotropic,title='Geneva-Copenhagen Survey'
k_end_print
END
