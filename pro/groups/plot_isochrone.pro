;+
;   NAME:
;      plot_isochrone
;   PURPOSE:
;      plot one of the isochrones corresponding to a certain age and
;      metallicity (this function is for debugging the isochrones
;      only)
;   INPUT:
;      age - log age
;      z - metallicity
;   KEYWORDS:
;      ms  - color the main-sequence
;      overplot - overplot
;      gzip - if set, gunzip and gzip the isochrones
;   OUTPUT:
;      plot to output device
;   REVISION HISTORY:
;      2009-11-10 - Written - Bovy (NYU)
;-
PRO PLOT_ISOCHRONE, age, z, ms=ms, overplot=overplot, gzip=gzip
isodir='isochrones/'
;;Read the isochrones
isoname= get_isoname(z,metal=zout)
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gunzip '+isodir+isoname+'.gz'
    spawn, cmd
ENDIF
iso= read_isochrone(isodir+isoname)
IF keyword_set(gzip) THEN BEGIN
    cmd= 'gzip '+isodir+isoname
    spawn, cmd
ENDIF
ages= iso.h01[UNIQ(iso.h01, SORT(iso.h01))]
ii=0L
WHILE ii LT n_elements(ages)-1 DO BEGIN
    IF ages[ii] GE age THEN BREAK
    ii+= 1
ENDWHILE
IF ii GE n_elements(ages)-1 THEN BEGIN
    isoage= ages[ii]
ENDIF ELSE IF ii EQ 0 THEN BEGIN
    isoage= ages[0]
ENDIF ELSE IF abs(age-ages[ii]) GT abs(age-ages[ii-1]) THEN BEGIN
    isoage= ages[ii-1]
ENDIF ELSE BEGIN
    isoage= ages[ii]
ENDELSE
thisiso_indx= where(iso.h01 EQ isoage)
nthisiso= n_elements(thisiso_indx)
thisiso= dblarr(2,nthisiso)
thisiso[0,*]= iso.h09[thisiso_indx]
thisiso[1,*]= iso.h10[thisiso_indx]
print, thisiso[0,0]-thisiso[1,0]
IF keyword_set(overplot) THEN djs_oplot, thisiso[0,*]-thisiso[1,*],thisiso[1,*], xtitle='B-V [mag]',ytitle='M_V [mag]',$
  yrange=[max(thisiso[1,*]),min(thisiso[1,*])] ELSE $
  djs_plot, thisiso[0,*]-thisiso[1,*],thisiso[1,*], xtitle='B-V [mag]',ytitle='M_V [mag]',$
  yrange=[max(thisiso[1,*]),min(thisiso[1,*])]
IF keyword_set(ms) THEN BEGIN
    ii=0
    isoBV= thisiso[0,*]-thisiso[1,*]
    WHILE isoBV[ii] GE isoBV[ii+5] DO ii+= 1
;    WHILE isoBV[ii] GE isoBV[ii+1] DO ii+= 1
    print, isoBV[ii], isoBV[ii+1]
    djs_oplot, thisiso[0,0:ii+15]-thisiso[1,0:ii+15],thisiso[1,0:ii+15], color=djs_icolor('red')
ENDIF

END
